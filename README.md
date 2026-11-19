# BioMy

BioMy is an auditable, evidence-first clinical and Ayurvedic intelligence prototype.
It combines deterministic laboratory calculations, structured clinical text parsing,
metabolic graph traversal, PageRank ranking, model-safe verification, and an
auditable report interface.

## Core systems

1. **UniversalBioAST**: converts typed lab text and exposures into typed nodes.
2. **Continuous scalar projector**: Fourier features for continuous measurements.
3. **Deterministic math**: Python implementations of common clinical equations.
4. **Clinical ripple engine**: directed herb/drug/metabolite graph traversal.
5. **Metabolic PageRank**: ranks graph bottlenecks.
6. **BioMistral adapter**: explanation-only model interface.
7. **Invariant gate**: rejects claims outside an approved evidence set.
8. **Clearance parity checker**: dose and renal metric validation.
9. **Auditable HUD**: structured findings, evidence, and disclaimer.

## Current MVP features

- Enzymes and transporters are first-class metabolic graph nodes.
- Ripple traversal follows herb -> enzyme/transporter -> drug chains.
- PageRank identifies graph bottlenecks.
- Prescription text accepts compact doses such as `5mg` and common frequencies.
- Streamlit provides separate lab, prescription, interaction, and evidence views.
- The seed interaction database is a demo dataset, not a validated clinical knowledge base.

## Evidence and scope fixes

- Graph screens require pair-specific evidence; shared-enzyme matches are labeled indirect mechanisms.
- Prescription parsing processes each line independently so dose and frequency fields do not bleed across medicines.
- Lab values are excluded from exposure parsing.
- Adult clinical equations validate age and reject non-finite or out-of-range inputs.

## TPSM status

TPSM is intentionally not included in this implementation. Its future demo should use one
published human interaction study, separately identified assumptions, and observed-versus-modeled points.

## Run

```powershell
python -m pip install -e ".[app,test]"
python -m pytest
python -m streamlit run app.py
```

## Training status

BioMistral fine-tuning is intentionally deferred until a verified dataset is
available. The model adapter currently provides a deterministic integration
boundary for later QLoRA work without claiming an already-trained model.

## No medical dependency claims

This repository does not replace a clinician or pharmacist. The included
interaction entries are seed examples and must be validated against current
peer-reviewed and authoritative clinical databases before use.
