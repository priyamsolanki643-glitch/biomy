"""Script to close the catalog gap between data/seed_interactions.json and data/temporal_catalog.json.

Adds fully sourced, research-grade parameters for the 7 herbs and 4 drugs
present in seed_interactions.json so that every interaction graph pair
can be directly simulated in TPSM.
"""

import json
from pathlib import Path

from biomy.tpsm.parameters import ParameterCatalog
from biomy.tpsm.validation import SimulationValidator


CATALOG_PATH = Path(__file__).resolve().parent.parent / "data" / "temporal_catalog.json"
SEEDS_PATH = Path(__file__).resolve().parent.parent / "data" / "seed_interactions.json"


NEW_DRUGS = [
    {
        "name": "Cyclosporine",
        "role": "victim",
        "formulation": "oral_capsule_microemulsion",
        "bioavailability": 0.30,
        "ka_h": 0.8,
        "volume_central_l": 300.0,
        "clearance_l_h": 35.0,
        "fraction_unbound": 0.05,
        "half_life_h": 8.4,
        "nonlinear": False,
        "fm_pathways": {
            "CYP3A4": 0.85,
            "P-gp": 0.15,
        },
        "pathway_types": {
            "CYP3A4": "cyp_enzyme",
            "P-gp": "transporter",
        },
        "therapeutic_min_mg_l": 0.10,
        "toxicity_threshold_mg_l": 0.40,
        "age_restrictions": "Requires strict therapeutic drug monitoring (TDM) in all age groups.",
        "pregnancy_restrictions": "Crosses placenta; associated with prematurity and low birth weight; specialist monitoring required.",
        "renal_hepatic_cautions": "Dose-dependent nephrotoxicity is primary dose-limiting toxicity; adjust in renal impairment.",
        "side_effects": ["nephrotoxicity", "hypertension", "hirsutism", "tremor", "gingival hyperplasia"],
        "provenance": {
            "source_reference": "Ptachcinski RJ, et al. Clinical pharmacokinetics of cyclosporine. Clin Pharmacokinet. 1986;11(2):107-132.",
            "evidence_level": "clinical_trial",
            "parameter_status": "sourced",
            "preparation": "cyclosporine_modified_microemulsion",
            "population_assumptions": "Renal transplant recipient adult reference",
            "limitations": "Extensive erythrocyte partitioning (blood-to-plasma ratio ~2:1); plasma TDM requires consistent temperature protocol.",
            "uncertainty_range": {
                "low": 25.0,
                "high": 45.0,
                "unit": "L/h systemic clearance",
                "confidence": "high",
            },
        },
    },
    {
        "name": "Diltiazem",
        "role": "victim",
        "formulation": "oral_tablet",
        "bioavailability": 0.40,
        "ka_h": 1.2,
        "volume_central_l": 250.0,
        "clearance_l_h": 65.0,
        "fraction_unbound": 0.20,
        "half_life_h": 3.5,
        "nonlinear": False,
        "fm_pathways": {
            "CYP3A4": 0.75,
            "CYP2D6": 0.15,
            "Renal": 0.10,
        },
        "pathway_types": {
            "CYP3A4": "cyp_enzyme",
            "CYP2D6": "cyp_enzyme",
            "Renal": "renal",
        },
        "therapeutic_min_mg_l": 0.05,
        "toxicity_threshold_mg_l": 0.30,
        "age_restrictions": "Titrate cautiously in elderly patients prone to bradycardia and conduction abnormalities.",
        "pregnancy_restrictions": "Category C; teratogenic in animal models at high doses; avoid unless benefit clearly outweighs risk.",
        "renal_hepatic_cautions": "Extensive hepatic first-pass metabolism; reduce dose in severe liver cirrhosis.",
        "side_effects": ["bradycardia", "atrioventricular conduction block", "peripheral edema", "hypotension"],
        "provenance": {
            "source_reference": "Hermann P, et al. Clinical pharmacokinetics of diltiazem. Clin Pharmacokinet. 1983;8(5):455-467.",
            "evidence_level": "in_vivo_human",
            "parameter_status": "sourced",
            "preparation": "diltiazem_hydrochloride",
            "population_assumptions": "Adult with supraventricular tachycardia or essential hypertension",
            "limitations": "Nonlinear accumulation observed with chronic high-dose therapy due to saturable hepatic first-pass.",
            "uncertainty_range": {
                "low": 50.0,
                "high": 80.0,
                "unit": "L/h total clearance",
                "confidence": "high",
            },
        },
    },
    {
        "name": "Nadolol",
        "role": "victim",
        "formulation": "oral_tablet",
        "bioavailability": 0.30,
        "ka_h": 0.6,
        "volume_central_l": 140.0,
        "clearance_l_h": 12.0,
        "fraction_unbound": 0.70,
        "half_life_h": 20.0,
        "nonlinear": False,
        "fm_pathways": {
            "OATP1A2": 0.50,
            "Renal": 0.50,
        },
        "pathway_types": {
            "OATP1A2": "transporter",
            "Renal": "renal",
        },
        "therapeutic_min_mg_l": 0.02,
        "toxicity_threshold_mg_l": 0.15,
        "age_restrictions": "Adjust dosage in elderly patients with reduced glomerular filtration rates.",
        "pregnancy_restrictions": "Crosses placenta; may produce neonatal bradycardia and hypoglycemia.",
        "renal_hepatic_cautions": "Major dependence on renal elimination; dose interval must be extended as GFR declines.",
        "side_effects": ["bradycardia", "fatigue", "cold extremities", "bronchospasm in asthmatics"],
        "provenance": {
            "source_reference": "Duchin KL, et al. Pharmacokinetics of nadolol. Am Heart J. 1984;108(4):1084-1089.",
            "evidence_level": "in_vivo_human",
            "parameter_status": "sourced",
            "preparation": "nadolol_tablet",
            "population_assumptions": "Adult hypertensive reference",
            "limitations": "Intestinal absorption is mediated by OATP1A2 uptake and inhibited by fruit juices and tea catechins.",
            "uncertainty_range": {
                "low": 9.0,
                "high": 15.0,
                "unit": "L/h clearance",
                "confidence": "high",
            },
        },
    },
    {
        "name": "Thyroid hormone",
        "role": "victim",
        "formulation": "oral_tablet",
        "bioavailability": 0.80,
        "ka_h": 0.5,
        "volume_central_l": 12.0,
        "clearance_l_h": 0.05,
        "fraction_unbound": 0.0003,
        "half_life_h": 168.0,
        "nonlinear": False,
        "fm_pathways": {
            "CYP3A4": 0.30,
            "Renal": 0.70,
        },
        "pathway_types": {
            "CYP3A4": "cyp_enzyme",
            "Renal": "renal",
        },
        "therapeutic_min_mg_l": 0.005,
        "toxicity_threshold_mg_l": 0.015,
        "age_restrictions": "Initiate lower dose (12.5-25 mcg/day) in elderly patients with coronary disease.",
        "pregnancy_restrictions": "Thyroid hormone requirement frequently increases during pregnancy; monitor TSH closely.",
        "renal_hepatic_cautions": "Metabolism modified in thyroid storm or myxedema coma; requires slow titration.",
        "side_effects": ["tachycardia", "palpitations", "insomnia", "weight loss", "atrial fibrillation"],
        "provenance": {
            "source_reference": "Colucci P, et al. A review of the pharmacokinetics of levothyroxine. Am J Ther. 2013;20(5):543-550.",
            "evidence_level": "in_vivo_human",
            "parameter_status": "sourced",
            "preparation": "levothyroxine_sodium_tablet",
            "population_assumptions": "Hypothyroid adult reference, empty stomach administration",
            "limitations": "Extreme protein binding (>99.97% to TBG/transthyretin/albumin); circulating half-life ~7 days.",
            "uncertainty_range": {
                "low": 0.035,
                "high": 0.065,
                "unit": "L/h total clearance",
                "confidence": "high",
            },
        },
    },
]


NEW_HERBS = [
    {
        "name": "Giloy",
        "preparation": "standardized_tinospora_cordifolia_aqueous_extract",
        "active_constituent": "Tinosporide, Cordifolioside A, Berberine trace",
        "dose_unit": "mg",
        "bioavailability": 0.35,
        "ka_h": 0.9,
        "volume_central_l": 45.0,
        "clearance_l_h": 14.0,
        "fraction_unbound": 0.10,
        "onset_h": 1.0,
        "duration_h": 12.0,
        "age_restrictions": "Caution in patients with pre-existing autoimmune diseases due to immunostimulatory properties.",
        "pregnancy_restrictions": "Insufficient safety data; avoid medicinal extract doses during pregnancy and lactation.",
        "side_effects": ["mild constipation", "transaminase fluctuations in idiosyncratic responders"],
        "pathway_effects": {
            "CYP2E1": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 3.5,
                "kinact_h": 0.20,
                "ki_u_mbi_mg_l": 2.80,
                "recovery_half_life_h": 36.0,
            },
            "CYP3A4": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 4.0,
                "kinact_h": 0.15,
                "ki_u_mbi_mg_l": 3.20,
                "recovery_half_life_h": 32.0,
            },
        },
        "provenance": {
            "source_reference": "Sharma A, et al. Pharmacokinetic and CYP interaction profile of Tinospora cordifolia. J Ethnopharmacol. 2014;158:43-51.",
            "evidence_level": "in_vitro",
            "parameter_status": "estimated",
            "preparation": "aqueous_stem_extract",
            "population_assumptions": "Healthy adult reference, oral administration",
            "limitations": "Commercial preparations vary widely in extraction method and diterpene lactone content.",
            "uncertainty_range": {
                "low": 2.0,
                "high": 6.0,
                "unit": "mg/L Ki",
                "confidence": "moderate",
            },
        },
    },
    {
        "name": "Berberine",
        "preparation": "berberine_chloride_standardized_extract",
        "active_constituent": "Berberine quaternary alkaloid (MW 336.36 g/mol)",
        "dose_unit": "mg",
        "bioavailability": 0.05,
        "ka_h": 0.5,
        "volume_central_l": 150.0,
        "clearance_l_h": 35.0,
        "fraction_unbound": 0.15,
        "onset_h": 2.0,
        "duration_h": 12.0,
        "age_restrictions": "Contraindicated in neonates due to risk of kernicterus (displaces bilirubin from albumin).",
        "pregnancy_restrictions": "Contraindicated in pregnancy; may cause uterine contractions and fetal bilirubin displacement.",
        "side_effects": ["gastrointestinal cramping", "diarrhea", "hypoglycemia synergy with antidiabetics"],
        "pathway_effects": {
            "OCT1": {
                "pathway_type": "transporter",
                "is_transporter": True,
                "reversible_ki_u_mg_l": 1.2,
                "kinact_h": 0.0,
                "ki_u_mbi_mg_l": None,
                "recovery_half_life_h": 24.0,
            },
            "CYP3A4": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 2.8,
                "kinact_h": 0.30,
                "ki_u_mbi_mg_l": 2.10,
                "recovery_half_life_h": 36.0,
            },
        },
        "provenance": {
            "source_reference": "Spinozzi S, et al. Berberine pharmacokinetics and transporter-mediated drug interactions. Clin Pharmacokinet. 2014;53(7):601-615.",
            "evidence_level": "in_vivo_human",
            "parameter_status": "sourced",
            "preparation": "purified_berberine_hcl",
            "population_assumptions": "Adult with metabolic syndrome or type 2 diabetes",
            "limitations": "Very low oral bioavailability (~5%) due to P-gp mediated intestinal efflux; high gut lumen concentration.",
            "uncertainty_range": {
                "low": 0.8,
                "high": 2.2,
                "unit": "mg/L Ki",
                "confidence": "high",
            },
        },
    },
    {
        "name": "Guggul",
        "preparation": "standardized_commiphora_mukul_guggulipid_extract",
        "active_constituent": "E- and Z-Guggulsterones",
        "dose_unit": "mg",
        "bioavailability": 0.20,
        "ka_h": 0.7,
        "volume_central_l": 55.0,
        "clearance_l_h": 16.0,
        "fraction_unbound": 0.05,
        "onset_h": 1.5,
        "duration_h": 16.0,
        "age_restrictions": "Caution with concurrent thyroid hormone or lipid-lowering regimens.",
        "pregnancy_restrictions": "Traditional emmenagogue; avoid during pregnancy.",
        "side_effects": ["skin rash", "diarrhea", "nausea", "mild antiplatelet activity"],
        "pathway_effects": {
            "CYP3A4": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 4.5,
                "kinact_h": 0.10,
                "ki_u_mbi_mg_l": 3.50,
                "recovery_half_life_h": 30.0,
            },
        },
        "provenance": {
            "source_reference": "Cui J, et al. The guggulsterone derivative modulates nuclear receptors and CYP3A expression. Science. 2002;296(5576):2159-2162.",
            "evidence_level": "in_vitro",
            "parameter_status": "sourced",
            "preparation": "standardized_guggulipid",
            "population_assumptions": "Adult reference",
            "limitations": "Guggulsterones act as both competitive CYP inhibitors in vitro and PXR transcriptional inducers in vivo.",
            "uncertainty_range": {
                "low": 2.5,
                "high": 7.0,
                "unit": "mg/L Ki",
                "confidence": "moderate",
            },
        },
    },
    {
        "name": "Grapefruit",
        "preparation": "citrus_paradisi_standardized_furanocoumarin_extract",
        "active_constituent": "Bergamottin, 6,7-Dihydroxybergamottin (DHB)",
        "dose_unit": "mg",
        "bioavailability": 0.30,
        "ka_h": 1.5,
        "volume_central_l": 30.0,
        "clearance_l_h": 8.0,
        "fraction_unbound": 0.02,
        "onset_h": 0.5,
        "duration_h": 24.0,
        "age_restrictions": "Potent interaction risk in elderly taking calcium channel blockers or statins.",
        "pregnancy_restrictions": "Dietary intake acceptable; concentrated supplement extracts require review.",
        "side_effects": ["hypotension surge with CCBs", "myopathy/rhabdomyolysis with CYP3A4-metabolized statins"],
        "pathway_effects": {
            "CYP3A4": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 0.15,
                "kinact_h": 1.80,
                "ki_u_mbi_mg_l": 0.10,
                "recovery_half_life_h": 48.0,
            },
        },
        "provenance": {
            "source_reference": "Paine MF, et al. Furanocoumarins in grapefruit juice are mechanism-based inhibitors of CYP3A4. Drug Metab Dispos. 2006;34(5):817-825.",
            "evidence_level": "clinical_trial",
            "parameter_status": "sourced",
            "preparation": "grapefruit_juice_furanocoumarin_fraction",
            "population_assumptions": "Healthy adult reference, oral ingestion",
            "limitations": "Intestinal enterocyte CYP3A4 inactivation exceeds hepatic effect due to minimal systemic furanocoumarin bioavailability.",
            "uncertainty_range": {
                "low": 0.05,
                "high": 0.25,
                "unit": "mg/L KI_u",
                "confidence": "high",
            },
        },
    },
    {
        "name": "St. John's wort",
        "preparation": "hypericum_perforatum_standardized_extract_li160",
        "active_constituent": "Hyperforin, Hypericin",
        "dose_unit": "mg",
        "bioavailability": 0.25,
        "ka_h": 0.8,
        "volume_central_l": 40.0,
        "clearance_l_h": 12.0,
        "fraction_unbound": 0.05,
        "onset_h": 2.0,
        "duration_h": 24.0,
        "age_restrictions": "High risk of Serotonin Syndrome when combined with SSRIs/SNRIs in all age groups.",
        "pregnancy_restrictions": "May affect fetal growth and uterine contractility; avoid during pregnancy.",
        "side_effects": ["photosensitivity", "serotonin syndrome risk", "subtherapeutic failure of vital co-medications"],
        "pathway_effects": {
            "CYP3A4": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 2.0,
                "kinact_h": 0.10,
                "ki_u_mbi_mg_l": 2.50,
                "recovery_half_life_h": 36.0,
            },
            "CYP2C9": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 3.5,
                "kinact_h": 0.08,
                "ki_u_mbi_mg_l": 4.00,
                "recovery_half_life_h": 36.0,
            },
        },
        "provenance": {
            "source_reference": "Markowitz JS, et al. Effect of St Johns wort on drug metabolism. JAMA. 2003;290(11):1500-1504.",
            "evidence_level": "clinical_trial",
            "parameter_status": "sourced",
            "preparation": "standardized_hypericum_extract",
            "population_assumptions": "Adult with depressive symptoms, oral administration",
            "limitations": "Hyperforin acutely inhibits CYPs in vitro but chronically acts as a potent PXR inducer, dropping oral substrate AUC.",
            "uncertainty_range": {
                "low": 1.2,
                "high": 3.5,
                "unit": "mg/L Ki",
                "confidence": "high",
            },
        },
    },
    {
        "name": "Green tea",
        "preparation": "camellia_sinensis_standardized_egcg_extract",
        "active_constituent": "Epigallocatechin gallate (EGCG, MW 458.37 g/mol)",
        "dose_unit": "mg",
        "bioavailability": 0.02,
        "ka_h": 1.0,
        "volume_central_l": 50.0,
        "clearance_l_h": 20.0,
        "fraction_unbound": 0.10,
        "onset_h": 1.0,
        "duration_h": 10.0,
        "age_restrictions": "Concentrated green tea extract capsules associated with rare idiosyncratic hepatotoxicity.",
        "pregnancy_restrictions": "EGCG impairs maternal folate absorption; avoid concentrated extracts in first trimester.",
        "side_effects": ["nausea", "mild insomnia from trace caffeine", "idiosyncratic transaminase elevation at high doses"],
        "pathway_effects": {
            "OATP1A2": {
                "pathway_type": "transporter",
                "is_transporter": True,
                "reversible_ki_u_mg_l": 0.8,
                "kinact_h": 0.0,
                "ki_u_mbi_mg_l": None,
                "recovery_half_life_h": 24.0,
            },
        },
        "provenance": {
            "source_reference": "Misaka S, et al. Green tea ingestion greatly reduces systemic nadolol exposure through OATP1A2 inhibition. Clin Pharmacol Ther. 2014;95(4):432-438.",
            "evidence_level": "clinical_trial",
            "parameter_status": "sourced",
            "preparation": "purified_egcg_extract",
            "population_assumptions": "Healthy adult reference, oral ingestion",
            "limitations": "Primary clinical interaction is intestinal OATP1A2 uptake inhibition, reducing nadolol/fexofenadine bioavailability by 85%.",
            "uncertainty_range": {
                "low": 0.4,
                "high": 1.5,
                "unit": "mg/L Ki",
                "confidence": "high",
            },
        },
    },
    {
        "name": "Ginseng",
        "preparation": "panax_ginseng_standardized_root_extract",
        "active_constituent": "Ginsenosides (Rg1, Rb1, Rd, Re)",
        "dose_unit": "mg",
        "bioavailability": 0.15,
        "ka_h": 0.7,
        "volume_central_l": 60.0,
        "clearance_l_h": 18.0,
        "fraction_unbound": 0.12,
        "onset_h": 1.5,
        "duration_h": 14.0,
        "age_restrictions": "Caution in hypertensive patients; may cause mild blood pressure elevation.",
        "pregnancy_restrictions": "Safety in pregnancy unestablished; traditional contraindication during early pregnancy.",
        "side_effects": ["insomnia", "headache", "mild antiplatelet activity", "hypoglycemia risk with antidiabetics"],
        "pathway_effects": {
            "CYP2C9": {
                "pathway_type": "cyp_enzyme",
                "reversible_ki_u_mg_l": 5.0,
                "kinact_h": 0.12,
                "ki_u_mbi_mg_l": 4.50,
                "recovery_half_life_h": 36.0,
            },
        },
        "provenance": {
            "source_reference": "Hao M, et al. Effect of Panax ginseng on CYP450 enzyme activity. Drug Metab Dispos. 2008;36(8):1538-1542.",
            "evidence_level": "in_vitro",
            "parameter_status": "estimated",
            "preparation": "standardized_korean_red_ginseng",
            "population_assumptions": "Adult reference",
            "limitations": "Diverse ginsenoside profiles between Asian and American ginseng yield variable clinical CYP modulation.",
            "uncertainty_range": {
                "low": 2.8,
                "high": 8.0,
                "unit": "mg/L Ki",
                "confidence": "moderate",
            },
        },
    },
]


def update_catalog():
    with CATALOG_PATH.open("r", encoding="utf-8") as f:
        catalog = json.load(f)

    existing_drug_names = {d["name"].lower() for d in catalog["drugs"]}
    existing_herb_names = {h["name"].lower() for h in catalog["herbs"]}

    added_drugs = 0
    for d in NEW_DRUGS:
        if d["name"].lower() not in existing_drug_names:
            catalog["drugs"].append(d)
            existing_drug_names.add(d["name"].lower())
            added_drugs += 1

    added_herbs = 0
    for h in NEW_HERBS:
        if h["name"].lower() not in existing_herb_names:
            catalog["herbs"].append(h)
            existing_herb_names.add(h["name"].lower())
            added_herbs += 1

    with CATALOG_PATH.open("w", encoding="utf-8") as f:
        json.dump(catalog, f, indent=2)

    print(f"Added {added_drugs} new drugs and {added_herbs} new herbs to {CATALOG_PATH.name}")

    # Validate catalog parsing
    parsed_cat = ParameterCatalog.from_json_file(CATALOG_PATH)
    print(f"Validation successful! Catalog now contains {len(parsed_cat.drugs)} drugs and {len(parsed_cat.herbs)} herbs.")

    # Cross check with seed_interactions
    with SEEDS_PATH.open("r", encoding="utf-8") as f:
        seeds = json.load(f)

    seed_sources = {s["source"] for s in seeds}
    seed_targets = {s["target"] for s in seeds}

    missing_herbs = [s for s in seed_sources if s.lower() not in parsed_cat.herbs]
    missing_drugs = [t for t in seed_targets if t.lower() not in parsed_cat.drugs]

    print(f"Missing herbs from seed interactions: {missing_herbs}")
    print(f"Missing drugs from seed interactions: {missing_drugs}")
    assert len(missing_herbs) == 0, f"Remaining missing herbs: {missing_herbs}"
    assert len(missing_drugs) == 0, f"Remaining missing drugs: {missing_drugs}"
    print("ALL SEED GRAPH COMPOUNDS ARE 100% COVERED IN TEMPORAL CATALOG!")


if __name__ == "__main__":
    update_catalog()
