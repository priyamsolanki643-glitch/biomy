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
