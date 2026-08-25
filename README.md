# Antibiotic AUC/MIC Optimizer
*Clinical Pharmacokinetics & Precision Antibiotic Dosing Optimization Engine*

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Python: 3.10+](https://img.shields.io/badge/Python-3.10+-brightgreen.svg)](https://python.org)
[![Tests: 100% Pass](https://img.shields.io/badge/Tests-21%20Passed-success.svg)]()

Antibiotic AUC/MIC Optimizer is a clinical pharmacokinetics/pharmacodynamics (PK/PD) decision support system designed to optimize exposure, maximize bactericidal efficacy, and prevent acute kidney injury (AKI) across critical antimicrobial classes.

---

## Pharmacokinetic Formulations & Mathematical Models

### 1. Vancomycin AUC24/MIC Precision Dosing
Engineered according to the **2020 Consensus Guidelines by ASHP, IDSA, PIDS, and SIDP** (Target $\text{AUC}_{24}/\text{MIC} = 400 - 600\text{ mg}\cdot\text{h/L}$ for MRSA):

- **Sawchuk-Zaske Paired Level Model** (Patient-Specific First-Order Elimination):
  $$k_{\text{el}} = \frac{\ln(C_{\text{peak}} / C_{\text{trough}})}{\Delta t}$$
- **Matzke Population Clearance Model** (Empiric / Level-Free):
  $$k_{\text{el}} = 0.00083 \times \text{CrCl} + 0.0044 \quad (\text{hr}^{-1})$$
  $$V_d = 0.70\text{ L/kg} \times \text{Weight}, \quad \text{CL} = k_{\text{el}} \times V_d$$
- **Steady-State Area Under the Curve**:
  $$\text{AUC}_{24} = \frac{\text{Total Daily Dose}}{\text{CL}}$$
- **Dose Adjustment Solver** (Scales to 500 mg*h/L midpoint in 250 mg increments):
  $$\text{Target Daily Dose} = \text{round}\left( \frac{500 \times \text{MIC} \times \text{CL}}{250} \right) \times 250$$
- **Nephrotoxicity Risk Grading**:
  - Subtherapeutic: $< 400\text{ mg}\cdot\text{h/L}$
  - Target Window: $400 - 600\text{ mg}\cdot\text{h/L}$
  - Elevated Risk: $601 - 800\text{ mg}\cdot\text{h/L}$
  - High AKI Risk: $> 800\text{ mg}\cdot\text{h/L}$

---

### 2. Trapezoidal AUC Integration
Calculates systemic exposure from serial serum concentration measurements:
- **Linear Trapezoidal Rule** (Ascending / Steady Phase):
  $$\text{AUC}_{1-2} = \frac{C_1 + C_2}{2} \times \Delta t$$
- **Log-Linear Trapezoidal Rule** (Descending Elimination Phase):
  $$\text{AUC}_{1-2} = \frac{C_1 - C_2}{\ln(C_1 / C_2)} \times \Delta t$$

---

### 3. Extended-Interval Aminoglycosides (Hartford Nomogram)
Concentration-dependent killing optimization:
- Target peak: $C_{\text{max}} / \text{MIC} \ge 8 - 10$
- High-dose single daily infusion ($5 - 7\text{ mg/kg}$ for Gentamicin/Tobramycin, $15\text{ mg/kg}$ for Amikacin).
- Interval selection ($\tau \in \{24\text{h}, 36\text{h}, 48\text{h}\}$) based on Cockcroft-Gault CrCl.
- Safety check confirming estimated trough $< 1.0\ \mu\text{g/mL}$.

---

### 4. Beta-Lactam %fT>MIC & Continuous Infusion Solver
Time-dependent bactericidal action:
- **Time Above MIC**:
  $$t_{>\text{MIC}} = \frac{\ln(f_u \cdot C_{\text{max}} / \text{MIC})}{k_{\text{el}}}$$
  $$\%fT_{>\text{MIC}} = \min\left(100.0, \frac{t_{>\text{MIC}}}{\tau} \times 100\right)$$
- **Continuous Infusion Rate Solver** (Targeting $C_{\text{ss}} = 4 \times \text{MIC}$):
  $$R_{\text{infusion}} = \frac{4 \times \text{MIC} \times \text{CL}}{f_u} \quad (\text{mg/hr})$$

---

### 5. Monte Carlo Probability of Target Attainment (PTA)
Simulates population PK variability ($V_d$ and $\text{CL}$ sampled from log-normal distributions, $\text{CV} = 25\%$) over $N \ge 1,000$ patient cohorts across discrete MIC breakpoints ($0.25 - 4.0\text{ mg/L}$).

---

## CLI Usage

### 1. Run Complete Benchmark PK/PD Dossier
```bash
python cli.py --demo
```
Output:
```text
============================================================================
      ANTIBIOTIC AUC/MIC OPTIMIZER - CLINICAL PK/PD DOSSIER          
============================================================================
 Patient ID : PT-BENCH-01 (60yo Male, 78.0kg, SCr: 1.2 mg/dL)
 Renal CrCl : 72.2 mL/min (IBW: 73.18 kg, Dosing Wt: 75.11 kg)
----------------------------------------------------------------------------
 1. VANCOMYCIN AUC24/MIC OPTIMIZATION (Target: 400 - 600 mg*h/L):
----------------------------------------------------------------------------
  Current Regimen    : 1250 mg q12h (Daily: 2500 mg)
  Steady-State AUC24 : 708.8 mg*h/L (Ratio: 708.8)
  Target Attained    : False
  Estimated Peak/Trg : Peak: 35.8 mg/L, Trough: 16.6 mg/L
  Risk Profile       : Elevated Risk (601 - 800 mg*h/L)
  Dose Recommendation: 1000 mg IV every 12 hours (Target Daily: 1750 mg)
----------------------------------------------------------------------------
 2. EXTENDED-INTERVAL AMINOGLYCOSIDE (Target Cmax/MIC >= 8-10):
----------------------------------------------------------------------------
  Regimen            : Gentamicin 526 mg (7.0 mg/kg) q24h
  Cmax / MIC Ratio   : 28.0 (Peak: 28.0 mg/L, Trough: 0.16 mg/L)
  Target Achieved    : True (Trough Safe: True)
----------------------------------------------------------------------------
 3. BETA-LACTAM TIME-DEPENDENT TARGET (%fT>MIC):
----------------------------------------------------------------------------
  Regimen            : Meropenem 1000 mg q8h
  %fT > MIC          : 59.8% (Target: >=40.0%)
  Continuous Infusion: 23.8 mg/hr (Loading Dose: 98 mg)
============================================================================
```

### 2. Calculate Patient Vancomycin AUC24/MIC
```bash
# Syntax: python cli.py --vanco <age> <gender: M/F> <weight_kg> <scr> <dose_mg> <interval_h> [mic]
python cli.py --vanco 60 M 78 1.2 1250 12 1.0
```

### 3. Batch Process Clinical CSV Scenarios
```bash
python cli.py --batch input_cases.csv output_evaluated.csv
```

### 4. Interactive Terminal
```bash
python cli.py --interactive
```

---

## Test Suite Execution

Run the complete 21-test suite with pure Python standard library:

```bash
python -m unittest test_auc_mic.py
python -m unittest discover -s tests
```

---

## License
MIT License. Developed for clinical pharmacokinetics research and antimicrobial stewardship decision support.
