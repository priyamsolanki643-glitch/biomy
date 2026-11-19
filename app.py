import json
from pathlib import Path

import streamlit as st

from biomy.ast_parser import BioASTParser
from biomy.calculations import BioCalculations
from biomy.pagerank import MetabolicPageRank
from biomy.ripple import ClinicalRippleEngine, Interaction


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

            matches = engine.ripple_frontier([herb])[herb]
            matching_records = [
                record
                for record in matches
                if drug.lower() in record["node"].lower()
                or drug.lower() in record["from"].lower()
            ]

            if matching_records:
                severities = {str(record["severity"]).lower() for record in matching_records}
                if "high" in severities:
                    st.error("High-severity evidence entry - professional pharmacist/clinician review.")
                elif "moderate" in severities:
                    st.warning("Moderate-severity evidence entry - verify against current guidance.")
                else:
                    st.info("Unknown or limited evidence - do not use as safety clearance.")

                st.subheader("Evidence-grounded findings")
                st.dataframe(matching_records, use_container_width=True)
            else:
                st.success("No match in the current seed database.")
                st.caption("No match means absence of evidence in this limited demo, not safety clearance.")

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
