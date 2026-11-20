import json
from pathlib import Path

import streamlit as st

from biomy.ast_parser import BioASTParser
from biomy.calculations import BioCalculations
from biomy.pagerank import MetabolicPageRank
from biomy.ripple import ClinicalRippleEngine, Interaction
from biomy.catalog import load_temporal_catalog
from biomy.tpsm import DoseEvent, TemporalPKStateMachine


st.set_page_config(page_title="BioMy", layout="wide")
st.title("BioMy - Clinical & Ayurvedic Intelligence Prototype")
st.warning("Educational prototype only. Not a diagnosis or treatment recommendation.")

tab_lab, tab_prescription, tab_interaction = st.tabs(
    ["Lab calculations", "Prescription parser", "Interaction graph"]
)

with tab_lab:
    left, right = st.columns(2)
    with left:
        st.subheader("Deterministic calculations")
        creatinine = st.number_input("Creatinine (mg/dL)", min_value=0.1, value=1.0)
        age = st.number_input("Age", min_value=1, value=40)
        weight = st.number_input("Weight (kg)", min_value=1.0, value=70.0)
        female = st.checkbox("Female")

        egfr = BioCalculations.egfr_ckd_epi_2021(creatinine, int(age), female)
        clearance = BioCalculations.creatinine_cockcroft_gault(creatinine, int(age), weight, female)
        st.metric("CKD-EPI eGFR", f"{egfr:.1f} mL/min/1.73m^2")
        st.metric("Creatinine clearance", f"{clearance:.1f} mL/min")

    with right:
        st.subheader("Laboratory interpretation")
        if egfr >= 90:
            st.success("eGFR: Normal or high range")
        elif egfr >= 60:
            st.info("eGFR: Mildly reduced range")
        elif egfr >= 30:
            st.warning("eGFR: Moderately reduced range")
        else:
            st.error("eGFR: Severely reduced range")
        st.caption("Screening aid only; kidney staging requires clinical context.")

with tab_prescription:
    prescription = st.text_area(
        "Type or paste prescription text",
        "Tab. Warfarin 5 mg OD\nPiperine extract 20 mg BD",
        height=140,
    )
    analytes, exposures = BioASTParser().parse(prescription)
    col_a, col_b = st.columns(2)
    with col_a:
        st.subheader("Parsed analytes")
        st.dataframe(
            [{"name": item.name, "value": item.value, "unit": item.unit} for item in analytes],
            use_container_width=True,
        )
    with col_b:
        st.subheader("Parsed exposures")
        st.dataframe(
            [
                {
                    "name": item.name,
                    "strength": item.strength,
                    "unit": item.unit,
                    "frequency": item.route,
                }
                for item in exposures
            ],
            use_container_width=True,
        )

with tab_interaction:
    st.subheader("Herb-drug interaction screen")
    left, right = st.columns([1, 2])
    with left:
        herb = st.text_input("Herb", "Piperine")
        drug = st.text_input("Drug", "Warfarin")
        check = st.button("Check interaction", type="primary")

    with right:
        if check:
            engine = ClinicalRippleEngine()
            seed_path = Path(__file__).parent / "data" / "seed_interactions.json"
            with seed_path.open(encoding="utf-8") as handle:
                entries = json.load(handle)
            engine.add_many(Interaction(**entry) for entry in entries)

            result = engine.screen_pair(herb, drug)

            if result["status"] == "direct_evidence":
                if result["severity"].lower() == "high":
                    st.error("High-severity evidence entry - professional pharmacist/clinician review.")
                elif result["severity"].lower() == "moderate":
                    st.warning("Moderate-severity evidence entry - verify against current guidance.")
                else:
                    st.info("Unknown or limited evidence - do not use as safety clearance.")
                st.subheader("Evidence-grounded findings")
            elif result["status"] == "indirect_mechanism_only":
                st.info("Indirect mechanism only - no direct pair-specific evidence in the seed database.")
                st.subheader("Hypothetical mechanism")
            else:
                st.success("No match in the current seed database.")
                st.caption("No match means absence of evidence in this limited demo, not safety clearance.")

            if result["status"] != "no_interaction_found":
                st.dataframe([result], use_container_width=True)

            ranking = MetabolicPageRank().rank(engine.graph)
            top = sorted(ranking.items(), key=lambda item: item[1], reverse=True)[:5]
            st.subheader("Metabolic graph centrality")
            st.dataframe(
                [{"node": node, "pagerank": round(score, 4)} for node, score in top],
                use_container_width=True,
            )

            dot = "digraph G {\n"
            dot += "  rankdir=LR;\n"
            for source, target, attrs in engine.graph.edges(data=True):
                effect = attrs.get("effect", "").replace('"', "")
                dot += f'  "{source}" -> "{target}" [label="{effect}"];\n'
            dot += "}"
            st.graphviz_chart(dot)

    st.caption(
        "Seed entries are illustrative and are not a substitute for PubMed review, "
        "prescribing information, or professional clinical assessment."
    )

st.divider()
st.header("Temporal Pharmacokinetic State Machine (TPSM)")
st.caption(
    "Research-grade temporal PK prototype: models time-dependent drug concentrations, first-order absorption, "
    "dynamic clearance, Mechanism-Based Inactivation (MBI), de novo enzyme turnover hysteresis, and non-linear saturation. "
    "Does not calculate safe clinical doses or replace professional therapeutic drug monitoring."
)

catalog_drugs, catalog_herbs, catalog_meta = load_temporal_catalog()

col_cfg1, col_cfg2, col_cfg3 = st.columns(3)
with col_cfg1:
    drug_name = st.selectbox(
        "Victim Drug",
        sorted(catalog_drugs),
        format_func=lambda value: catalog_drugs[value].name,
    )
    herb_name = st.selectbox(
        "Herbal Perpetrator (optional)",
        ["None"] + sorted(catalog_herbs),
        format_func=lambda value: "No herb" if value == "None" else catalog_herbs[value].name,
    )

with col_cfg2:
    comp_choice = st.selectbox(
        "Disposition Model",
        ["1-compartment", "2-compartment"],
        index=0,
    )
    solver_choice = st.selectbox(
        "ODE Solver Method",
        ["Radau", "BDF", "RK45"],
        index=0,
        help="Radau is recommended for stiff PK/PD differential equations with disparate time constants.",
    )

with col_cfg3:
    duration = st.slider("Simulation duration (hours)", 12, 168, 72, step=12)
    step_size = st.select_slider("Discretization step (hours)", options=[0.1, 0.25, 0.5, 1.0], value=0.5)

st.subheader("Dosing Schedule")
col_dose1, col_dose2, col_dose3, col_dose4 = st.columns(4)
with col_dose1:
    dose = st.number_input("Drug dose (mg)", min_value=0.01, value=5.0)
with col_dose2:
    drug_freq = st.selectbox(
        "Drug frequency",
        ["Single dose (t=0)", "Every 6 hours (Q6H)", "Every 12 hours (BID)", "Every 24 hours (QD)"],
        index=3,
    )
with col_dose3:
    herb_dose = st.number_input("Herb dose (mg)", min_value=0.0, value=20.0 if herb_name != "None" else 0.0)
with col_dose4:
    herb_regimen = st.selectbox(
        "Herb regimen",
        ["Single dose (t=0)", "Co-administered QD", "Delayed initiation (t=24h)"],
        index=1,
    )

if st.button("Simulate Timeline", type="primary"):
    selected_drugs = [catalog_drugs[drug_name]]
    selected_herbs = [] if herb_name == "None" else [catalog_herbs[herb_name]]

    # Build events list according to regimen
    events = []
    # Drug events
    if drug_freq == "Single dose (t=0)":
        events.append(DoseEvent(name=catalog_drugs[drug_name].name, dose_mg=dose, time_h=0.0))
    elif drug_freq == "Every 6 hours (Q6H)":
        for t in range(0, int(duration), 6):
            events.append(DoseEvent(name=catalog_drugs[drug_name].name, dose_mg=dose, time_h=float(t)))
    elif drug_freq == "Every 12 hours (BID)":
        for t in range(0, int(duration), 12):
            events.append(DoseEvent(name=catalog_drugs[drug_name].name, dose_mg=dose, time_h=float(t)))
    elif drug_freq == "Every 24 hours (QD)":
        for t in range(0, int(duration), 24):
            events.append(DoseEvent(name=catalog_drugs[drug_name].name, dose_mg=dose, time_h=float(t)))

    # Herb events
    if selected_herbs and herb_dose > 0:
        herb_obj = selected_herbs[0]
        if herb_regimen == "Single dose (t=0)":
            events.append(DoseEvent(name=herb_obj.name, dose_mg=herb_dose, time_h=0.0))
        elif herb_regimen == "Co-administered QD":
            for t in range(0, int(duration), 24):
                events.append(DoseEvent(name=herb_obj.name, dose_mg=herb_dose, time_h=float(t)))
        elif herb_regimen == "Delayed initiation (t=24h)":
            for t in range(24, int(duration), 24):
                events.append(DoseEvent(name=herb_obj.name, dose_mg=herb_dose, time_h=float(t)))

    from biomy.tpsm.engine import TPSMEngine
    from biomy.tpsm.models import CompartmentMode
    from biomy.tpsm.serialization import TPSMSerializer

    comp_mode = (
        CompartmentMode.TWO_COMPARTMENT
        if comp_choice == "2-compartment"
        else CompartmentMode.ONE_COMPARTMENT
    )

    engine = TPSMEngine(
        drugs=selected_drugs,
        herbs=selected_herbs,
        compartment_mode=comp_mode,
        solver_method=solver_choice,
    )

    try:
        result = engine.simulate(events=events, duration_h=float(duration), step_h=float(step_size))
        d_key = catalog_drugs[drug_name].name

        # Metric summary row
        m_col1, m_col2, m_col3, m_col4 = st.columns(4)
        with m_col1:
            st.metric("Peak Concentration (Cmax)", f"{result.c_max.get(d_key, 0.0):.4f} mg/L")
        with m_col2:
            st.metric("Trough Concentration (Cmin)", f"{result.c_min.get(d_key, 0.0):.4f} mg/L")
        with m_col3:
            st.metric("AUC (0 -> t)", f"{result.auc.get(d_key, 0.0):.4f} mg*h/L")
        with m_col4:
            first_breach = result.first_toxicity_breach_time.get(d_key)
            if first_breach is not None:
                st.metric("Toxicity Threshold Breach", f"t = {first_breach:.1f} h")
            else:
                st.metric("Toxicity Threshold Breach", "None detected")

        # Clinical warning banner if threshold crossed
        if first_breach is not None:
            st.error(
                f"Modelled exposure crossed the configured research threshold "
                f"({selected_drugs[0].toxicity_threshold_mg_l} mg/L) at t = {first_breach:.1f}h; "
                f"professional review is required."
            )

        # Plot 1: Drug & Herb Concentration Curves
        st.subheader("Systemic Concentration Profiles")
        chart_data = {"Time (h)": result.time_points}
        for d in selected_drugs:
            if d.name in result.concentration_curves:
                chart_data[f"{d.name} (mg/L)"] = result.concentration_curves[d.name]
        for h in selected_herbs:
            if h.name in result.concentration_curves:
                chart_data[f"{h.name} (mg/L)"] = result.concentration_curves[h.name]
        st.line_chart(chart_data, x="Time (h)")

        # Plot 2: Active Enzyme Fraction (Recovery Hysteresis)
        if result.enzyme_activity_curves:
            st.subheader("Metabolic Pathway & Enzyme Activity Fraction (0.0 to 1.0)")
            st.caption("Demonstrates Mechanism-Based Inactivation (MBI) depletion and post-washout de novo synthesis recovery lag.")
            enz_chart_data = {"Time (h)": result.time_points}
            for enz, curve in result.enzyme_activity_curves.items():
                enz_chart_data[f"{enz} Activity"] = curve
            st.line_chart(enz_chart_data, x="Time (h)")

        # Expanders for Provenance, Assumptions, LLM Serializer
        with st.expander("Parameter Provenance & Evidence Base", expanded=False):
            st.json(result.parameter_provenance)

        with st.expander("Biophysical Assumptions & Solver Diagnostics", expanded=False):
            st.write("**Model Assumptions:**")
            for a in result.assumptions:
                st.markdown(f"- {a}")
            if result.warnings:
                st.write("**Clinical Considerations:**")
                for w in result.warnings:
                    st.warning(w)
            st.write("**Uncertainty Metadata:**")
            st.json(result.uncertainty_metadata)

        with st.expander("Structured BioMistral Prompt Context (Zero-Hallucination Guard)", expanded=False):
            st.caption("Verbatim deterministic summary emitted for LLM consumption. The LLM explains this output and cannot recalculate numbers.")
            st.code(TPSMSerializer.to_llm_prompt_context(result), language="markdown")

    except Exception as exc:
        st.error(f"Simulation failed: {exc}")
