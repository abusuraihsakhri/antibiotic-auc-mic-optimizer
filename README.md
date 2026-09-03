# Antibiotic AUC/MIC Precision Pharmacokinetics Optimizer

A pure Python clinical pharmacokinetics, therapeutic drug monitoring (TDM), and precision antimicrobial dosing engine implementing:
- **Cockcroft-Gault Creatinine Clearance & Devine Anthropometrics:**
  $$\text{CrCl} = \frac{(140 - \text{Age}) \times \text{Weight}_{\text{dosing}}}{72 \times \text{SCr}} \times (0.85 \text{ if female})$$
  - Utilizes Devine Ideal Body Weight (IBW) and Adjusted Body Weight ($ABW_{0.4}$) when actual weight exceeds $120\%$ IBW.
- **Vancomycin $\text{AUC}_{24}/\text{MIC}$ Bayesian/Analytical Precision Solver:**
  - 2020 ASHP/IDSA/PIDS/SIDP Consensus Guidelines targeting steady-state $\text{AUC}_{24}/\text{MIC}_{\text{BMD}} = 400 - 600\text{ mg}\cdot\text{h/L}$.
  - Sawchuk-Zaske paired level elimination rate: $k_{el} = \frac{\ln(C_{\text{peak}} / C_{\text{trough}})}{\Delta t}$.
  - Matzke population clearance fallback: $k_{el} = 0.00083 \times \text{CrCl} + 0.0044\text{ h}^{-1}$, $V_d = 0.70\text{ L/kg}$.
  - Acute Kidney Injury (AKI) Nephrotoxicity risk stratification: Low ($<400$), Optimal ($400-600$), Elevated ($601-800$), Severe ($>800$).
  - Automatic dose titration scaling regimens to the $500\text{ mg}\cdot\text{h/L}$ consensus midpoint in $250\text{ mg}$ steps.
- **Aminoglycoside Concentration-Dependent Optimization:**
  - Hartford high-dose extended-interval nomogram targeting $C_{\text{max}}/\text{MIC} \ge 8 - 10$ with complete trough clearance.
- **Beta-Lactam Time-Dependent Target Attainment:**
  - Evaluates $\%fT_{>\text{MIC}}$ and calculates continuous infusion delivery rates: $R_{\text{inf}} = \frac{(\text{Target} \times \text{MIC}) \times \text{CL}}{f_u}$.
- **Monte Carlo Probability of Target Attainment (PTA):**
  - Log-normal stochastic population parameter sampling across patient cohorts to quantify target attainment percentages.
- **High-Throughput Batch Patient CSV Processing:** Ingests electronic health record datasets for ward and ICU surveillance.

Requires Python standard library only (zero external runtime dependencies).

---

## Precision PK/PD Consensus Targets

| Drug Class | Pharmacodynamic Driver | Clinical Target | Clinical Endpoint |
|:-----------|:-----------------------|:----------------|:------------------|
| **Vancomycin** | $\text{AUC}_{24}/\text{MIC}$ | $400 - 600\text{ mg}\cdot\text{h/L}$ | MRSA bacteremia eradication & nephrotoxicity minimization |
| **Aminoglycosides** | $C_{\text{max}}/\text{MIC}$ | $\ge 8 - 10$ | Concentration-dependent bactericidal kill & post-antibiotic effect |
| **Beta-Lactams** | $\%fT_{>\text{MIC}}$ | $40\% - 100\%$ of dosing interval | Time-dependent cell wall synthesis inhibition |
| **Fluoroquinolones** | $\text{AUC}_{24}/\text{MIC}$ | $\ge 125$ (Gram-negative), $\ge 30$ (S. pneumoniae) | Bacterial eradication & prevention of resistance |

---

## Features

- **ASHP/IDSA 2020 Consensus Compliance:** Replaces trough-only monitoring with rigorous area-under-the-curve optimization.
- **Linear & Log-Linear Trapezoidal Integration:** Exact non-compartmental integration over empirical serum time-concentration curves.
- **Nephrotoxicity Risk Mitigation:** Proactively alerts clinicians to steep increases in odds of acute kidney injury when exposure exceeds $600\text{ mg}\cdot\text{h/L}$.
- **Batch CSV Cohort Processing:** Rapidly processes cohorts of ICU/ward patients with automated regimen recommendations.

---

## Installation & Requirements

- Python 3.10+ (tested on 3.10, 3.11, 3.12)
- Zero external runtime dependencies. `pytest` is optional for running tests.

```bash
git clone https://github.com/abusuraihsakhri/antibiotic-auc-mic-optimizer.git
cd antibiotic-auc-mic-optimizer
```

---

## CLI Usage

### 1. Run Complete PK/PD Benchmark Demonstration
```bash
python cli.py --demo --json
```

### 2. Individual Vancomycin Regimen Assessment
```bash
# Format: <age> <gender M/F> <weight_kg> <serum_creatinine> <dose_mg> <interval_hours> [mic]
python cli.py --vanco 65 M 78.0 1.2 1250 12 1.0
```

### 3. Batch Process Patient CSV
```bash
python cli.py --batch sample.csv results.csv
```

---

## Python API Quickstart

```python
from auc_mic import PatientProfile, optimize_vancomycin

patient = PatientProfile(
    patient_id="ICU-402",
    age_years=62,
    gender="Male",
    weight_kg=80.0,
    height_cm=178.0,
    serum_creatinine_mg_dl=1.1,
)

res = optimize_vancomycin(patient, dose_mg=1250.0, interval_hours=12.0, mic_mg_l=1.0)

print(f"Calculated CrCl: {res.renal_crcl_ml_min} mL/min")
print(f"Steady-state AUC24: {res.auc24_mg_h_l} mg*h/L")
print(f"Target Achieved: {res.target_achieved}")
print(f"Risk Band: {res.nephrotoxicity_risk}")
print(f"Recommended Regimen: {res.recommended_regimen}")
```

---

## Running Tests

Run the test suite using standard `unittest` or `pytest`:

```bash
pytest -v
```
