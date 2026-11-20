# BioMy: Auditable Ayurvedic Herb-Drug Clinical Intelligence & Research-Grade TPSM

BioMy is an auditable, evidence-first clinical and Ayurvedic intelligence prototype.
It combines deterministic laboratory calculations, structured clinical text parsing,
metabolic interaction graph traversal, PageRank bottleneck analysis, model-safe invariant
gates, an auditable HUD, and a genuinely complete, research-grade **Temporal Pharmacokinetic
State Machine (TPSM)**.

---

## Core Systems

1. **UniversalBioAST (`biomy/ast_parser.py`)**: Converts clinical prescriptions and laboratory text into typed entity representations.
2. **Continuous Scalar Projector (`biomy/calculations.py`)**: Fourier features for continuous physiological measurements.
3. **Deterministic Math (`biomy/calculations.py`)**: Exact implementations of CKD-EPI (2021) eGFR and Cockcroft-Gault CrCl equations with clinical boundary gates.
4. **Clinical Ripple Engine (`biomy/ripple.py`)**: Directed herb -> enzyme/transporter -> drug causal cascade traversal.
5. **Metabolic PageRank (`biomy/pagerank.py`)**: Centrality scoring to highlight metabolic bottlenecks and high-risk convergence nodes.
6. **Biochemical Invariant Gate & Parity Checker (`biomy/verification.py`)**: Rejects model hallucinations and validates mass-balance invariants.
7. **Research-Grade TPSM (`biomy/tpsm/`)**: Stiff ODE-based multi-compartment pharmacokinetic/pharmacodynamic simulation engine.
8. **Auditable Streamlit HUD (`app.py`)**: Visualizes deterministic lab scores, interaction networks, and time-dependent PK/PD curves.

---

## Research-Grade TPSM Architecture

The Temporal Pharmacokinetic State Machine (`biomy/tpsm/`) is a research-grade biophysical simulation engine designed to compute the biological state dynamically. The LLM (e.g., BioMistral) serves strictly to explain TPSM outputs, with zero authority to alter numerical concentrations or timestamps.

### Module Structure

```text
biomy/tpsm/
    __init__.py          # Public API exports
    models.py            # Typed dataclasses, enums, ProvenanceRecord, SimulationResult
    parameters.py        # ParameterCatalog loader with uncertainty mapping
    events.py            # EventTimelineProcessor for exact impulse segmentation
    equations.py         # StateLayout & PKPDEquations (linear/nonlinear/MBI/reversible)
    solver.py            # ODESolverWrapper (scipy.integrate.solve_ivp Radau/BDF/RK45)
    engine.py            # TPSMEngine simulation coordinator & backward-compatible facade
    metrics.py           # PKMetricsCalculator (AUC, Cmax, Cmin, TMR, toxicity root finding)
    validation.py        # SimulationValidator boundary and kinetic constraint checks
    serialization.py     # TPSMSerializer for zero-hallucination BioMistral prompt context
```

### Mathematical Formulation

#### 1. Linear Disposition (1-Compartment & 2-Compartment)
$$\frac{dA_g}{dt} = -k_a \cdot A_g$$
$$\frac{dA_c}{dt} = k_a \cdot A_g - \frac{CL(t)}{V_c} A_c - Q \left(\frac{A_c}{V_c} - \frac{A_p}{V_p}\right)$$
$$\frac{dA_p}{dt} = Q \left(\frac{A_c}{V_c} - \frac{A_p}{V_p}\right)$$
$$C(t) = \frac{A_c(t)}{V_c}$$

#### 2. Michaelis-Menten Non-Linear Elimination
$$\frac{dA_c}{dt} = k_a \cdot A_g - \frac{V_{\max}(t) \cdot C(t)}{K_m + C(t)} \cdot V_c$$

#### 3. Mechanism-Based Enzyme Inactivation (MBI) & De Novo Turnover
$$\frac{dE_{\text{rel}}}{dt} = k_{\text{deg}} \cdot (1 - E_{\text{rel}}) - \frac{k_{\text{inact}} \cdot I_u}{K_{I_u} + I_u} E_{\text{rel}}$$
$$k_{\text{deg}} = \frac{\ln 2}{t_{1/2, \text{recovery}}}$$

#### 4. Dynamic Clearance Scaling (Reversible + Irreversible Inhibition)
$$CL_{\text{scale}}(t) = f_{m,0} + \sum_j \frac{f_{m,j} \cdot E_{\text{rel},j}(t)}{1 + \frac{I_{u,j}(t)}{K_{i,\text{rev},j}}}$$
$$CL(t) = CL_{\text{baseline}} \cdot CL_{\text{scale}}(t)$$

#### 5. Analytical Toxicity Threshold Crossing
Uses root-finding event functions with `direction = 1.0`:
$$g(t, y) = \frac{A_c(t)}{V_c} - C_{\text{tox}} = 0$$

---

## Target Compounds & Parameter Provenance

All parameters in `data/temporal_catalog.json` feature explicit `provenance` records including:
- `source_reference`: Peer-reviewed clinical literature (Clin Pharmacokinet, Drug Metab Dispos, etc.)
- `evidence_level`: `clinical_trial`, `in_vivo_human`, `in_vitro`, or `estimated`
- `parameter_status`: `sourced`, `estimated`, or `illustrative`
- `uncertainty_range`: Low/high bounds with confidence ratings
- `population_assumptions` & `limitations`

### Starter Drugs
1. **Paracetamol** (Linear, CYP2E1/glucuronidation)
2. **Ibuprofen** (Linear, CYP2C9)
3. **Warfarin** (Narrow-therapeutic-index victim, CYP2C9/CYP3A4)
4. **Phenytoin** (Michaelis-Menten non-linear saturation, CYP2C9/CYP2C19)
5. **Metformin** (Renal elimination, OCT1/OCT2 transporters)
6. **Amlodipine** (2-compartment open disposition, CYP3A4)
7. **Midazolam** (Sensitive index substrate, CYP3A4)
8. **Atorvastatin** (2-compartment disposition, CYP3A4 + OATP1B1)
9. **Digoxin** (Narrow-therapeutic-index, P-glycoprotein/renal)

### Starter Herbs & Botanicals
1. **Piperine** (Black Pepper bioenhancer; potent CYP3A4/CYP2C9 MBI and P-gp inhibitor)
2. **Ashwagandha** (*Withania somnifera*; moderate CYP3A4/CYP2C19 modulator)
3. **Curcumin** (*Curcuma longa*; CYP2C9/CYP3A4 inhibitor)
4. **Ginger** (*Zingiber officinale*; mild CYP2C9/CYP3A4 modulator)
5. **Garlic** (*Allium sativum*; CYP3A4 modulator and antiplatelet synergy)

---

## Benchmark Validation Scenarios

Standardized validation scenarios are defined in `data/temporal_scenarios.json`:
- **Scenario A**: Warfarin + Piperine (MBI recovery hysteresis and post-washout bleeding hazard)
- **Scenario B**: Phenytoin + Curcumin (Michaelis-Menten non-linear saturation and runaway accumulation)
- **Scenario C**: Midazolam + Piperine (Acute reversible probe substrate exposure elevation)
- **Scenario D**: Paracetamol Q6H (Multi-dose steady-state accumulation and washout)
- **Scenario E**: Digoxin (Narrow-therapeutic-index clearance reduction and TMR decline)

---

## Installation & Running

```bash
# Clone and install with testing & UI dependencies
git clone https://github.com/priyamsolanki643-glitch/biomy.git
cd biomy
python -m pip install -e ".[app,test]"

# Run full research test suite (56 unit, regression, and scenario tests)
python -m pytest

# Launch interactive Streamlit interface
python -m streamlit run app.py
```

---

## Research Disclaimer

> **IMPORTANT CLINICAL NOTICE**
> BioMy and its Temporal Pharmacokinetic State Machine (TPSM) are computational research prototypes developed for scientific inquiry, systems biology modeling, and educational demonstration.
>
> - TPSM models biophysical mechanisms but **is NOT an approved medical device or prescribing system**.
> - It does NOT recommend doses, predict individual patient toxicities, or advise altering medical regimens.
> - Parameters represent reference populations and require formal clinical therapeutic drug monitoring (TDM) and toxicological oversight for actual patients.
