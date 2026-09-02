# Antibiotic AUC MIC Optimizer

> **Domain:** Cardiovascular Medicine & Hemodynamic Analytics  
> **Reference Guidelines & Standards:** `AHA/ACC Practice Guidelines & ESC Clinical Standards`

<div align="center">

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
![Python](https://img.shields.io/badge/Python-3.10%20%7C%203.11%20%7C%203.12-3776AB.svg?logo=python&logoColor=white)
![FastAPI](https://img.shields.io/badge/FastAPI-0.111-009688.svg?logo=fastapi&logoColor=white)
![Audit Trail](https://img.shields.io/badge/Audit-HMAC--SHA256_Tamper--Evident-brightgreen.svg)
![Zero-PHI Guard](https://img.shields.io/badge/Guard-Zero--PHI_Outbound-blue.svg)
![Docker](https://img.shields.io/badge/Docker-Ready-2496ED.svg?logo=docker&logoColor=white)

</div>

---

## 📖 What It Does

Antibiotic AUC/MIC Optimizer & Precision Pharmacokinetics Engine.

Mathematical & Clinical Formulas Implemented:
1. Cockcroft-Gault Creatinine Clearance (CrCl) with Devine IBW & Adjusted BW:
   - CrCl = [(140 - Age) * DosingWeight / (72 * SCr)] * (0.85 if Female)
2. Vancomycin AUC24/MIC Optimization (2020 ASHP/IDSA Consensus Target 400-600 mg*h/L):
   - Sawchuk-Zaske paired level elimination: kel = ln(C_peak / C_trough) / dt
   - Matzke population elimination fallback: kel = 0.00083 * CrCl + 0.0044 (/hr)
   - Volume of Distribution: Vd = 0.70 L/kg
   - Clearance: CL = kel * Vd
   - Steady-State AUC24 = DailyDose / CL
   - Target Dose Scaler to 500 mg*h/L midpoint in 250 mg increments
   - AKI nephrotoxicity risk grading (Low <400, Target 400-600, Elevated 600-800, High >800)
3. Linear & Log-Linear Trapezoidal AUC Integration:
   - Linear trapezoid: AUC = 0.5 * (C1 + C2) * dt
   - Log-linear elimination trapezoid: AUC = (C1 - C2) / ln(C1 / C2) * dt
4. Aminoglycoside Cmax/MIC Concentration-Dependent Optimization:
   - Target Cmax/MIC >= 8 - 10 (Hartford high-dose extended-interval nomogram)
5. Fluoroquinolone AUC24/MIC Optimization:
   - Target AUC24/MIC >= 125 for Gram-negative bacteremia, >= 30 for S. pneumoniae
6. Beta-Lactam Time-Dependent (%fT>MIC) & Continuous Infusion Solver:
   - t_>MIC = ln(fu * Cmax / MIC) / kel; %fT>MIC = min(100, t_>MIC / tau * 100)
   - Continuous infusion rate: R_inf = (Target_x_MIC * MIC * CL) / fu
7. Monte Carlo Probability of Target Attainment (PTA) & Cumulative Fraction of Response (CFR):
   - Log-normal stochastic parameter sampling across patient cohorts and MIC distributions.

Author: Dr. Abu Suraih Sakhri
License: MIT

AUC/MIC optimization enrichment features for antibiotic-auc-mic-optimizer.

Implements three high-impact items from specifications on vancomycin
pharmacokinetics (linear one-compartment model):

    AUC24 at steady state = total daily dose / CL, with CL = kel x Vd.
    kel from paired levels: kel = ln(Cpeak/Ctrough) / dt.
    Population fallback (Matzke): kel = 0.00083 x CrCl + 0.0044 per hour,
    Vd = 0.7 L/kg.

1. AUC/MIC target optimization against the 2021 ASHP/IDSA consensus window
   of 400-600 mg*h/L for MRSA, with automatic dose scaling to the 500
   midpoint in 250 mg steps.
2. Nephrotoxicity risk quantification: graded AKI-risk bands rising above
   AUC24 > 600, with the risk gradient reported alongside attainment.
3. Monte Carlo probability-of-target-attainment simulation over population PK
   variability, plus a two-point vs multi-point sampling accuracy comparison.

Author: Dr. Abu Suraih Sakhri
License: MIT

---

## ⚙️ Key Capabilities & Algorithmic Modules

### 🔬 Core Algorithmic & Evaluation Engines

- **`NephrotoxicityBand`** — dedicated module for nephrotoxicity band evaluation and state verification.
- **`PatientProfile`** — dedicated module for patient profile evaluation and state verification.
- **`VancoAUCResult`** — dedicated module for vanco a u c result evaluation and state verification.
- **`AminoglycosideResult`** — dedicated module for aminoglycoside result evaluation and state verification.
- **`BetaLactamResult`** — dedicated module for beta lactam result evaluation and state verification.
- **`MonteCarloPTAResult`** — dedicated module for monte carlo p t a result evaluation and state verification.

---

## 📐 Mathematical Formulation & Logic

```text
  Mathematical & Clinical Formulas Implemented:
  """Devine formula for Ideal Body Weight in kg."""
  ibw = calculate_ibw(gender, height_cm) if height_cm else weight_kg
  """Matzke population formula: kel = 0.00083 * CrCl + 0.0044 (/hr)."""
  kel = calculate_matzke_kel(renal["crcl_ml_min"])
```

---

## 💻 CLI Quickstart & Usage

### 1. Guided Interactive Mode
```bash
python cli.py
```

### 2. Direct Parameterized Evaluation
```bash
python cli.py --demo <value> --interactive <value> --json <value> --batch <value>
```

### Parameter Reference
- `--demo`: Specifies input measurement or parameter value.
- `--interactive`: Specifies input measurement or parameter value.
- `--json`: Specifies input measurement or parameter value.
- `--batch`: Specifies input measurement or parameter value.
- `--vanco`: Specifies input measurement or parameter value.

### Input Data Schema

| Field | Description | Requirement |
|:------|:------------|:------------|
| `id` | Parameter / observation metric | Required |
| `value` | Parameter / observation metric | Required |
| `qty` | Parameter / observation metric | Required |

---

## 🛡️ Security & Enterprise Architecture

* **Zero-PHI Outbound Interceptor:** Active AST and regex inspection blocking SSNs, MRNs, phone numbers, and patient identifiers.
* **Tamper-Evident HMAC-SHA256 Audit Trail:** Chained, cryptographically signed logs for every evaluation and state transition.
* **Air-Gapped LLM Reasoning Adapter:** Agnostic integration for local Ollama instances (`llama3`, `mistral`), Claude 3.5 Sonnet, GPT-4o, and deterministic test mocks.
* **Active Learning Bayesian Calibration:** Dynamic tracker updating worker reliability weights and monitoring Brier calibration drift.
* **FastAPI & Prometheus Telemetry:** Exposes OpenAPI 3.1 REST endpoints and operational Prometheus metrics (`/metrics`).

---

## 🧪 Testing & Verification

Run the automated test suite:

```bash
pytest -v
```

Execute high-throughput batch simulation benchmarks:

```bash
python simulator.py --tasks 1000 --concurrency 8
```

---

## 🐳 Container Deployment

```bash
docker build -t antibiotic-auc-mic-optimizer .
docker run -p 8000:8000 antibiotic-auc-mic-optimizer
```
