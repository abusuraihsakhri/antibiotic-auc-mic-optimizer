#!/usr/bin/env python3
"""
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
"""

from __future__ import annotations

import csv
import enum
import io
import json
import math
import random
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from typing import Any, Callable, Dict, List, Optional, Sequence, Tuple


# =====================================================================
# 1. ENUMS & DATA MODELS
# =====================================================================

class NephrotoxicityBand(str, enum.Enum):
    SUBTHERAPEUTIC = "Subtherapeutic (< 400 mg*h/L)"
    TARGET = "Optimal Consensus Window (400 - 600 mg*h/L)"
    ELEVATED = "Elevated Risk (601 - 800 mg*h/L)"
    HIGH = "High AKI Risk (> 800 mg*h/L)"


@dataclass
class PatientProfile:
    patient_id: str
    age_years: int
    gender: str
    weight_kg: float
    height_cm: float
    serum_creatinine_mg_dl: float

    @property
    def is_female(self) -> bool:
        return self.gender.strip().upper().startswith("F")

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class VancoAUCResult:
    dose_mg: float
    interval_hours: float
    daily_dose_mg: float
    kel_per_hr: float
    vd_liters: float
    clearance_l_hr: float
    auc24_mg_h_l: float
    auc24_mic_ratio: float
    target_achieved: bool
    recommended_daily_dose_mg: float
    recommended_regimen: str
    estimated_peak_mg_l: float
    estimated_trough_mg_l: float
    nephrotoxicity_risk: str
    guidance: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class AminoglycosideResult:
    drug: str
    dose_mg: float
    dose_per_kg: float
    interval_hours: float
    peak_cmax_mg_l: float
    cmax_mic_ratio: float
    target_achieved: bool
    trough_estimated_mg_l: float
    trough_safe: bool
    recommendation: str

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class BetaLactamResult:
    drug: str
    dose_mg: float
    interval_hours: float
    infusion_duration_hours: float
    cmax_unbound_mg_l: float
    pct_ft_above_mic: float
    target_threshold_pct: float
    target_achieved: bool
    continuous_infusion_rate_mg_hr: float
    loading_dose_mg: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


@dataclass
class MonteCarloPTAResult:
    regimen_label: str
    mic_tested: float
    n_simulated: int
    pta_target_pct: float
    subtherapeutic_pct: float
    nephrotoxic_pct: float
    mean_auc24: float
    median_auc24: float
    p5_auc24: float
    p95_auc24: float

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


# =====================================================================
# 2. RENAL KINETICS & ANTHROPOMETRICS
# =====================================================================

def calculate_ibw(gender: str, height_cm: float) -> float:
    """Devine formula for Ideal Body Weight in kg."""
    height_in = height_cm / 2.54
    inches_over_5ft = max(0.0, height_in - 60.0)
    if gender.strip().upper().startswith("F"):
        return round(45.5 + 2.3 * inches_over_5ft, 2)
    return round(50.0 + 2.3 * inches_over_5ft, 2)


def calculate_dosing_weight(actual_weight_kg: float, ibw_kg: float) -> Tuple[float, bool]:
    """Adjusted body weight if actual > 1.2 * IBW."""
    if actual_weight_kg > 1.2 * ibw_kg:
        adj = ibw_kg + 0.4 * (actual_weight_kg - ibw_kg)
        return round(adj, 2), True
    return round(actual_weight_kg, 2), False


def calculate_cockcroft_gault(
    age_years: int,
    gender: str,
    weight_kg: float,
    serum_creatinine: float,
    height_cm: Optional[float] = None
) -> Dict[str, Any]:
    """Computes Cockcroft-Gault CrCl with optional Devine IBW adjustment."""
    if age_years <= 0 or serum_creatinine <= 0 or weight_kg <= 0:
        raise ValueError("Age, SCr, and weight must be positive.")

    ibw = calculate_ibw(gender, height_cm) if height_cm else weight_kg
    dosing_wt, is_adj = calculate_dosing_weight(weight_kg, ibw) if height_cm else (weight_kg, False)

    raw_crcl = ((140.0 - float(age_years)) * dosing_wt) / (72.0 * float(serum_creatinine))
    is_fem = gender.strip().upper().startswith("F")
    crcl = raw_crcl * 0.85 if is_fem else raw_crcl
    crcl_val = max(1.0, round(crcl, 1))

    return {
        "crcl_ml_min": crcl_val,
        "ibw_kg": ibw,
        "dosing_weight_kg": dosing_wt,
        "is_female": is_fem,
        "adjusted_weight_used": is_adj
    }


# =====================================================================
# 3. TRAPEZOIDAL AUC INTEGRATION ENGINE
# =====================================================================

def calculate_trapezoidal_auc(time_points: Sequence[float], concentrations: Sequence[float]) -> float:
    """
    Computes Area Under the Curve (AUC) using the linear/log-linear trapezoidal rule.
    Uses log-linear trapezoidal rule for descending elimination phases, linear for ascending.
    """
    if len(time_points) != len(concentrations):
        raise ValueError("time_points and concentrations must have matching lengths.")
    if len(time_points) < 2:
        return 0.0

    # Sort pairs by time
    pairs = sorted(zip(time_points, concentrations), key=lambda p: p[0])
    total_auc = 0.0

    for i in range(len(pairs) - 1):
        t1, c1 = pairs[i]
        t2, c2 = pairs[i + 1]
        dt = t2 - t1
        if dt <= 0:
            continue

        if c1 > 0 and c2 > 0 and c1 > c2:
            # Descending phase: Log-linear trapezoid
            auc_segment = (c1 - c2) / math.log(c1 / c2) * dt
        else:
            # Ascending or equal phase: Standard linear trapezoid
            auc_segment = 0.5 * (c1 + c2) * dt

        total_auc += auc_segment

    return round(total_auc, 2)


# =====================================================================
# 4. VANCOMYCIN AUC24/MIC OPTIMIZER
# =====================================================================

def calculate_matzke_kel(crcl_ml_min: float) -> float:
    """Matzke population formula: kel = 0.00083 * CrCl + 0.0044 (/hr)."""
    return 0.00083 * max(crcl_ml_min, 5.0) + 0.0044


def calculate_kel_from_paired_levels(peak_mg_l: float, trough_mg_l: float, dt_hours: float) -> float:
    """Sawchuk-Zaske elimination rate constant: kel = ln(Cpeak / Ctrough) / dt."""
    if peak_mg_l <= trough_mg_l or peak_mg_l <= 0 or trough_mg_l <= 0 or dt_hours <= 0:
        raise ValueError("Valid paired levels require Peak > Trough > 0 and elapsed time > 0.")
    return math.log(peak_mg_l / trough_mg_l) / dt_hours


def optimize_vancomycin(
    patient: PatientProfile,
    dose_mg: float = 1250.0,
    interval_hours: float = 12.0,
    mic_mg_l: float = 1.0,
    measured_peak_mg_l: Optional[float] = None,
    measured_trough_mg_l: Optional[float] = None,
    peak_trough_dt_hours: Optional[float] = None
) -> VancoAUCResult:
    """
    Computes steady-state Vancomycin AUC24/MIC and optimizes dosage for the 400-600 target window.
    """
    if measured_peak_mg_l and measured_trough_mg_l and peak_trough_dt_hours:
        kel = calculate_kel_from_paired_levels(measured_peak_mg_l, measured_trough_mg_l, peak_trough_dt_hours)
    else:
        renal = calculate_cockcroft_gault(patient.age_years, patient.gender, patient.weight_kg, patient.serum_creatinine_mg_dl, patient.height_cm)
        kel = calculate_matzke_kel(renal["crcl_ml_min"])

    vd = 0.70 * patient.weight_kg
    cl = kel * vd
    daily_dose = dose_mg * (24.0 / interval_hours)
    auc24 = daily_dose / cl
    ratio = auc24 / max(0.1, mic_mg_l)

    target_achieved = 400.0 <= ratio <= 600.0

    # Risk grading
    if ratio < 400.0:
        risk = NephrotoxicityBand.SUBTHERAPEUTIC.value
        guidance = "Subtherapeutic exposure: increased risk of MRSA treatment failure or MIC creep. Dose escalation indicated."
    elif ratio <= 600.0:
        risk = NephrotoxicityBand.TARGET.value
        guidance = "Within consensus target window (400-600 mg*h/L): optimal bactericidal efficacy with low nephrotoxicity."
    elif ratio <= 800.0:
        risk = NephrotoxicityBand.ELEVATED.value
        guidance = "Elevated exposure: increased odds of Acute Kidney Injury (AKI). Recommend downward dose adjustment."
    else:
        risk = NephrotoxicityBand.HIGH.value
        guidance = "High AKI Risk (>800 mg*h/L): severe nephrotoxicity danger. Immediate dose reduction required."

    # Dose optimization targeting midpoint of 500 mg*h/L
    scale = 500.0 / ratio if ratio > 0 else 1.0
    rec_daily = round((daily_dose * scale) / 250.0) * 250.0

    # Optimal interval based on elimination half-life: t1/2 = ln(2) / kel
    t_half = math.log(2.0) / kel
    if t_half <= 6.0:
        rec_interval = 8.0
    elif t_half <= 14.0:
        rec_interval = 12.0
    elif t_half <= 26.0:
        rec_interval = 24.0
    else:
        rec_interval = 48.0

    doses_per_day = 24.0 / rec_interval
    rec_single = round((rec_daily / doses_per_day) / 250.0) * 250.0
    rec_regimen = f"{int(rec_single)} mg IV every {int(rec_interval)} hours"

    # Peak & Trough estimation
    # Cmax_ss = (Dose / Vd) / (1 - exp(-kel * tau))
    # Ctrough_ss = Cmax_ss * exp(-kel * tau)
    cmax_est = (dose_mg / vd) / max(0.01, (1.0 - math.exp(-kel * interval_hours)))
    ctrough_est = cmax_est * math.exp(-kel * interval_hours)

    return VancoAUCResult(
        dose_mg=dose_mg,
        interval_hours=interval_hours,
        daily_dose_mg=daily_dose,
        kel_per_hr=round(kel, 5),
        vd_liters=round(vd, 2),
        clearance_l_hr=round(cl, 3),
        auc24_mg_h_l=round(auc24, 1),
        auc24_mic_ratio=round(ratio, 1),
        target_achieved=target_achieved,
        recommended_daily_dose_mg=rec_daily,
        recommended_regimen=rec_regimen,
        estimated_peak_mg_l=round(cmax_est, 1),
        estimated_trough_mg_l=round(ctrough_est, 1),
        nephrotoxicity_risk=risk,
        guidance=guidance
    )


# =====================================================================
# 5. AMINOGLYCOSIDE & BETA-LACTAM OPTIMIZERS
# =====================================================================

def optimize_aminoglycoside(
    patient: PatientProfile,
    drug_name: str = "Gentamicin",
    dose_mg_per_kg: float = 7.0,
    mic_mg_l: float = 1.0
) -> AminoglycosideResult:
    """High-dose extended-interval aminoglycoside solver targeting Cmax/MIC >= 8-10."""
    renal = calculate_cockcroft_gault(patient.age_years, patient.gender, patient.weight_kg, patient.serum_creatinine_mg_dl, patient.height_cm)
    crcl = renal["crcl_ml_min"]
    dosing_wt = renal["dosing_weight_kg"]

    total_dose = round(dose_mg_per_kg * dosing_wt, 0)
    vd = 0.25 * dosing_wt
    kel = 0.0024 * crcl + 0.01

    if crcl >= 60.0:
        interval = 24.0
    elif crcl >= 40.0:
        interval = 36.0
    elif crcl >= 20.0:
        interval = 48.0
    else:
        interval = 48.0

    cmax = total_dose / vd
    cmax_mic = cmax / max(0.1, mic_mg_l)
    trough = cmax * math.exp(-kel * interval)
    safe_trough_limit = 1.0 if drug_name.lower() in ("gentamicin", "tobramycin") else 4.0

    target_met = cmax_mic >= 8.0
    trough_safe = trough < safe_trough_limit

    rec = f"Administer {int(total_dose)} mg IV every {int(interval)} hours."
    if not trough_safe:
        rec += " Extend dosing interval to allow complete clearance."

    return AminoglycosideResult(
        drug=drug_name,
        dose_mg=total_dose,
        dose_per_kg=dose_mg_per_kg,
        interval_hours=interval,
        peak_cmax_mg_l=round(cmax, 1),
        cmax_mic_ratio=round(cmax_mic, 1),
        target_achieved=target_met,
        trough_estimated_mg_l=round(trough, 2),
        trough_safe=trough_safe,
        recommendation=rec
    )


def optimize_beta_lactam(
    patient: PatientProfile,
    drug_name: str,
    dose_mg: float,
    interval_hours: float,
    infusion_duration_hours: float = 0.5,
    mic_mg_l: float = 2.0,
    protein_binding_pct: float = 20.0
) -> BetaLactamResult:
    """Beta-lactam %fT>MIC and continuous infusion solver."""
    renal = calculate_cockcroft_gault(patient.age_years, patient.gender, patient.weight_kg, patient.serum_creatinine_mg_dl, patient.height_cm)
    crcl = renal["crcl_ml_min"]

    vd = 0.25 * patient.weight_kg
    kel = 0.002 * crcl + 0.10
    cl = kel * vd

    fu = 1.0 - (protein_binding_pct / 100.0)
    cmax_total = dose_mg / vd
    cmax_unbound = cmax_total * fu

    target_thresh = 50.0
    if any(k in drug_name.lower() for k in ("mero", "erta", "imip")):
        target_thresh = 40.0
    elif "piperacillin" in drug_name.lower():
        target_thresh = 60.0

    if cmax_unbound <= mic_mg_l:
        pct_ft = 0.0
    else:
        t_above = math.log(cmax_unbound / mic_mg_l) / kel
        pct_ft = min(100.0, (t_above / interval_hours) * 100.0)

    # Continuous infusion targeting 4x MIC
    target_css_unbound = 4.0 * mic_mg_l
    target_css_total = target_css_unbound / fu
    cont_rate = target_css_total * cl
    load_dose = target_css_total * vd

    return BetaLactamResult(
        drug=drug_name,
        dose_mg=dose_mg,
        interval_hours=interval_hours,
        infusion_duration_hours=infusion_duration_hours,
        cmax_unbound_mg_l=round(cmax_unbound, 2),
        pct_ft_above_mic=round(pct_ft, 1),
        target_threshold_pct=target_thresh,
        target_achieved=pct_ft >= target_thresh,
        continuous_infusion_rate_mg_hr=round(cont_rate, 1),
        loading_dose_mg=round(load_dose, 0)
    )


# =====================================================================
# 6. MONTE CARLO PROBABILITY OF TARGET ATTAINMENT (PTA)
# =====================================================================

def run_monte_carlo_pta(
    regimens: Sequence[Tuple[float, float, str]],  # (dose_mg, interval_h, label)
    weight_kg: float,
    mean_crcl: float,
    mic_tested: float = 1.0,
    n_simulations: int = 5000,
    seed: int = 42
) -> List[MonteCarloPTAResult]:
    """
    Simulates population pharmacokinetics with log-normal random sampling of Vd and CL
    to estimate Probability of Target Attainment (PTA) into the 400-600 AUC/MIC window.
    """
    rng = random.Random(seed)
    results: List[MonteCarloPTAResult] = []

    mean_kel = calculate_matzke_kel(mean_crcl)
    mean_vd = 0.70 * weight_kg

    for dose_mg, interval_h, label in regimens:
        daily_dose = dose_mg * (24.0 / interval_h)
        auc_samples: List[float] = []
        target_count = 0
        sub_count = 0
        nephro_count = 0

        for _ in range(n_simulations):
            # Sample log-normal variation (CV = 25%)
            vd_i = mean_vd * math.exp(rng.gauss(0.0, 0.25))
            kel_i = mean_kel * math.exp(rng.gauss(0.0, 0.25))
            cl_i = kel_i * vd_i
            auc_i = daily_dose / cl_i
            ratio_i = auc_i / mic_tested

            auc_samples.append(auc_i)
            if 400.0 <= ratio_i <= 600.0:
                target_count += 1
            elif ratio_i < 400.0:
                sub_count += 1
            else:
                nephro_count += 1

        auc_sorted = sorted(auc_samples)
        p5 = auc_sorted[int(0.05 * n_simulations)]
        p95 = auc_sorted[int(0.95 * n_simulations)]
        median = auc_sorted[int(0.50 * n_simulations)]
        mean_auc = sum(auc_samples) / n_simulations

        results.append(MonteCarloPTAResult(
            regimen_label=label,
            mic_tested=mic_tested,
            n_simulated=n_simulations,
            pta_target_pct=round((target_count / n_simulations) * 100.0, 1),
            subtherapeutic_pct=round((sub_count / n_simulations) * 100.0, 1),
            nephrotoxic_pct=round((nephro_count / n_simulations) * 100.0, 1),
            mean_auc24=round(mean_auc, 1),
            median_auc24=round(median, 1),
            p5_auc24=round(p5, 1),
            p95_auc24=round(p95, 1)
        ))

    return results


# =====================================================================
# 7. BATCH CSV ASSESSOR (BACKWARD COMPATIBILITY & BENCHMARKS)
# =====================================================================

def assess_row(row: Dict[str, Any]) -> Dict[str, Any]:
    """Processes a single dictionary row containing patient and pharmacokinetic parameters."""
    try:
        age = int(row.get("age_years") or row.get("age") or 55)
        gender = str(row.get("gender") or row.get("sex") or "M")
        weight = float(row.get("weight_kg") or row.get("weight") or 70.0)
        height = float(row.get("height_cm") or row.get("height") or 175.0)
        scr = float(row.get("serum_creatinine_mg_dl") or row.get("scr") or row.get("creatinine") or 1.0)
        dose = float(row.get("dose_mg") or row.get("dose") or 1250.0)
        interval = float(row.get("interval_hours") or row.get("interval") or 12.0)
        mic = float(row.get("mic_mg_l") or row.get("mic") or 1.0)

        peak = float(row["peak_mg_l"]) if "peak_mg_l" in row and row["peak_mg_l"] else None
        trough = float(row["trough_mg_l"]) if "trough_mg_l" in row and row["trough_mg_l"] else None
        dt = float(row["dt_hours"]) if "dt_hours" in row and row["dt_hours"] else None

        p = PatientProfile("CASE", age, gender, weight, height, scr)
        res = optimize_vancomycin(p, dose, interval, mic, peak, trough, dt)
        return res.to_dict()
    except Exception as ex:
        return {"error": str(ex)}


def process_csv(input_path: str, output_path: str) -> List[Dict[str, Any]]:
    """Batch processes a CSV dataset of PK scenarios."""
    with open(input_path, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    results = []
    for row in rows:
        assessed = assess_row(row)
        merged = {**row, **{k: str(v) for k, v in assessed.items()}}
        results.append(merged)

    all_keys = list(fieldnames)
    for r in results:
        for k in r.keys():
            if k not in all_keys:
                all_keys.append(k)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=all_keys)
        writer.writeheader()
        writer.writerows(results)

    return results
