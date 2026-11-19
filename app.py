import json
from pathlib import Path

import streamlit as st

from biomy.calculations import BioCalculations
from biomy.ripple import ClinicalRippleEngine, Interaction
from biomy.verification import ClearanceParityChecker


st.set_page_config(page_title="BioMy", layout="wide")
st.title("BioMy — Clinical & Ayurvedic Intelligence Prototype")
st.warning("Educational prototype only. Not a diagnosis or treatment recommendation.")

left, right = st.columns(2)
with left:
    st.subheader("Lab calculations")
    creatinine = st.number_input("Creatinine (mg/dL)", min_value=0.1, value=1.0)
    age = st.number_input("Age", min_value=1, value=40)
    weight = st.number_input("Weight (kg)", min_value=1.0, value=70.0)
    female = st.checkbox("Female")

    egfr = BioCalculations.egfr_ckd_epi_2021(creatinine, int(age), female)
    clearance = BioCalculations.creatinine_cockcroft_gault(creatinine, int(age), weight, female)
    st.metric("CKD-EPI eGFR", f"{egfr:.1f} mL/min/1.73m²")
    st.metric("Creatinine clearance", f"{clearance:.1f} mL/min")

with right:
    st.subheader("Herb-drug interaction screen")
    herb = st.text_input("Herb", "Piperine")
    drug = st.text_input("Drug", "Warfarin")
    if st.button("Check"):
        with open(Path(__file__).parent / "data" / "seed_interactions.json", encoding="utf-8") as handle:
            entries = json.load(handle)
        engine = ClinicalRippleEngine()
        engine.add_many(Interaction(**entry) for entry in entries)
        matches = engine.ripple_frontier([herb])[herb]
        st.write(matches if matches else "No verified entry in the seed database.")
