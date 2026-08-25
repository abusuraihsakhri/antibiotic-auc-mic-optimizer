#!/usr/bin/env python3
"""
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
"""

from __future__ import annotations

import math
import random
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple


def cockcroft_gault_crcl(age_years: int, weight_kg: float, scr_mg_dl: float,
                         female: bool = False) -> float:
    crcl = ((140 - age_years) * weight_kg) / (72.0 * scr_mg_dl)
    if female:
        crcl *= 0.85
    return max(crcl, 5.0)


def kel_from_levels(peak_mg_l: float, trough_mg_l: float, dt_hours: float) -> float:
    if peak_mg_l <= trough_mg_l or peak_mg_l <= 0 or trough_mg_l <= 0:
        raise ValueError("need peak > trough > 0")
    return math.log(peak_mg_l / trough_mg_l) / dt_hours


def matzke_kel(crcl_ml_min: float) -> float:
    return 0.00083 * crcl_ml_min + 0.0044


@dataclass(frozen=True)
class AUCEstimate:
    auc24_mg_h_per_l: float
    auc24_over_mic: float
    within_400_600: bool
    recommended_daily_dose_mg: int


def estimate_auc24(daily_dose_mg: float, clearance_l_hr: float,
                   mic_mg_l: float = 1.0) -> AUCEstimate:
    """Exact steady-state linear-PK identity: AUC_tau = Dose/CL."""
    auc24 = daily_dose_mg / clearance_l_hr
    ratio = auc24 / mic_mg_l
    within = 400.0 <= ratio <= 600.0
    scaled_daily = daily_dose_mg * (500.0 * mic_mg_l / auc24)
    recommended = int(round(scaled_daily / 250.0) * 250.0)
    return AUCEstimate(
        auc24_mg_h_per_l=round(auc24, 1),
        auc24_over_mic=round(ratio, 1),
        within_400_600=within,
        recommended_daily_dose_mg=recommended,
    )


def nephrotoxicity_risk(auc24_mg_h_per_l: float) -> Dict[str, object]:
    """Consensus gradient: AKI odds climb as AUC24 exceeds ~600 mg*h/L."""
    if auc24_mg_h_per_l < 400:
        band, note = "low", "subtherapeutic range; efficacy risk dominates"
    elif auc24_mg_h_per_l <= 600:
        band, note = "acceptable", "consensus target window; monitor SCr q48h"
    elif auc24_mg_h_per_l <= 800:
        band, note = ("elevated", "above 600: increased AKI odds; tighten monitoring "
                                    "and reassess dose necessity")
    else:
        band, note = "high", "AUC24 > 800 strongly associated with AKI; reduce exposure"
    return {"band": band, "guidance": note}


@dataclass(frozen=True)
class MonteCarloResult:
    regimen_label: str
    pta_pct: float
    nephrotoxic_exposure_pct: float


def monte_carlo_pta(regimens: Sequence[Tuple[float, str]], weight_kg: float,
                    mean_crcl: float, mic_mg_l: float,
                    n_patients: int = 10_000, seed: int = 42) -> List[MonteCarloResult]:
    """Sample Vd and CL variability; report PTA into the 400-600 ratio window."""
    rng = random.Random(seed)
    results: List[MonteCarloResult] = []
    mean_cl = matzke_kel(mean_crcl) * 0.7 * weight_kg

    for daily_dose, label in regimens:
        attained = 0
        toxic_exposure = 0
        for _ in range(n_patients):
            vd = rng.gauss(0.7 * weight_kg, 0.07 * weight_kg)
            cl = max(rng.gauss(mean_cl, 0.30 * mean_cl), 0.5)
            auc24 = daily_dose / cl
            if 400.0 <= auc24 / mic_mg_l <= 600.0:
                attained += 1
            if auc24 > 600.0:
                toxic_exposure += 1
        results.append(MonteCarloResult(
            regimen_label=label,
            pta_pct=round(100.0 * attained / n_patients, 1),
            nephrotoxic_exposure_pct=round(100.0 * toxic_exposure / n_patients, 1),
        ))
    return sorted(results, key=lambda r: r.pta_pct, reverse=True)


def multipoint_vs_twopoint_bias(true_kel: float, true_vd_liters: float,
                                daily_dose_mg: float, mic_mg_l: float,
                                seed: int = 7) -> Dict[str, float]:
    """Simulate a monoexponential profile; compare 2-point vs 3-point k recovery."""
    rng = random.Random(seed)
    cmax = (daily_dose_mg / true_vd_liters)

    def level(t: float, noise_frac: float = 0.05) -> float:
        truth = cmax * math.exp(-true_kel * t)
        return truth * (1.0 + rng.gauss(0.0, noise_frac))

    def estimate_kel(times: List[float]) -> float:
        xs = times
        ys = [math.log(level(t)) for t in times]
        mean_x = sum(xs) / len(xs)
        mean_y = sum(ys) / len(ys)
        sxx = sum((x - mean_x) ** 2 for x in xs)
        sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
        return -sxy / sxx

    cl_true = true_kel * true_vd_liters
    auc_true = daily_dose_mg / cl_true

    k2 = estimate_kel([2.0, 8.0])
    k3 = estimate_kel([1.0, 4.0, 9.0])

    def pct_error(k_est: float) -> float:
        auc_est = daily_dose_mg / (k_est * true_vd_liters)
        return round((auc_est - auc_true) / auc_true * 100.0, 2)

    return {
        "true_auc24": round(auc_true, 1),
        "twopoint_auc_bias_pct": pct_error(k2),
        "threepoint_auc_bias_pct": pct_error(k3),
    }


def _demo() -> None:
    crcl = cockcroft_gault_crcl(age_years=64, weight_kg=72, scr_mg_dl=1.0)
    kel = kel_from_levels(28.0, 11.0, 6.0)
    vd = 0.7 * 72.0
    cl = kel * vd
    est = estimate_auc24(daily_dose_mg=3000.0, clearance_l_hr=cl, mic_mg_l=1.0)
    print({"crcl": round(crcl), "kel": round(kel, 4), "auc24": est.auc24_mg_h_per_l,
           "auc_over_mic": est.auc24_over_mic, "in_target": est.within_400_600,
           "recommended_daily_mg": est.recommended_daily_dose_mg})
    print(nephrotoxicity_risk(est.auc24_mg_h_per_l))

    regimens = [(1500.0, "750 mg q12h"), (2000.0, "1 g q12h"),
                (2500.0, "1.25 g q12h"), (3000.0, "1.5 g q12h")]
    for mc in monte_carlo_pta(regimens, weight_kg=72.0, mean_crcl=crcl, mic_mg_l=1.0):
        print({"regimen": mc.regimen_label, "pta_pct": mc.pta_pct,
               "toxic_exposure_pct": mc.nephrotoxic_exposure_pct})

    print(multipoint_vs_twopoint_bias(true_kel=kel, true_vd_liters=vd,
                                      daily_dose_mg=3000.0, mic_mg_l=1.0))


if __name__ == "__main__":
    _demo()
