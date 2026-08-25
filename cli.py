#!/usr/bin/env python3
"""
Antibiotic AUC/MIC Optimizer - Command Line Interface (CLI)

Precision dosing optimization for Vancomycin, Aminoglycosides, Beta-Lactams,
Trapezoidal AUC integration, Monte Carlo PTA simulations, and batch CSV processing.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from typing import List

from auc_mic import (
    NephrotoxicityBand,
    PatientProfile,
    assess_row,
    calculate_cockcroft_gault,
    calculate_trapezoidal_auc,
    optimize_aminoglycoside,
    optimize_beta_lactam,
    optimize_vancomycin,
    process_csv,
    run_monte_carlo_pta,
)


def run_demo_evaluation(json_out: bool = False) -> None:
    """Executes benchmark clinical scenarios demonstrating the PK/PD optimizer capabilities."""
    patient = PatientProfile(
        patient_id="PT-BENCH-01",
        age_years=60,
        gender="Male",
        weight_kg=78.0,
        height_cm=178.0,
        serum_creatinine_mg_dl=1.2
    )

    vanco = optimize_vancomycin(patient, dose_mg=1250.0, interval_hours=12.0, mic_mg_l=1.0)
    genta = optimize_aminoglycoside(patient, drug_name="Gentamicin", dose_mg_per_kg=7.0, mic_mg_l=1.0)
    mero = optimize_beta_lactam(patient, drug_name="Meropenem", dose_mg=1000.0, interval_hours=8.0, mic_mg_l=1.0)
    pta = run_monte_carlo_pta(
        regimens=[(1000, 12, "1000mg q12h"), (1250, 12, "1250mg q12h"), (1500, 12, "1500mg q12h")],
        weight_kg=78.0,
        mean_crcl=65.0,
        mic_tested=1.0,
        n_simulations=2000
    )

    if json_out:
        out = {
            "patient": patient.to_dict(),
            "vancomycin": vanco.to_dict(),
            "aminoglycoside": genta.to_dict(),
            "beta_lactam": mero.to_dict(),
            "monte_carlo_pta": [p.to_dict() for p in pta]
        }
        print(json.dumps(out, indent=2))
        return

    print("=" * 76)
    print("      ANTIBIOTIC AUC/MIC OPTIMIZER - CLINICAL PK/PD DOSSIER          ")
    print("=" * 76)
    print(f" Patient ID : {patient.patient_id} ({patient.age_years}yo {patient.gender}, {patient.weight_kg}kg, SCr: {patient.serum_creatinine_mg_dl} mg/dL)")
    crcl_info = calculate_cockcroft_gault(patient.age_years, patient.gender, patient.weight_kg, patient.serum_creatinine_mg_dl, patient.height_cm)
    print(f" Renal CrCl : {crcl_info['crcl_ml_min']} mL/min (IBW: {crcl_info['ibw_kg']} kg, Dosing Wt: {crcl_info['dosing_weight_kg']} kg)")
    print("-" * 76)
    print(" 1. VANCOMYCIN AUC24/MIC OPTIMIZATION (Target: 400 - 600 mg*h/L):")
    print("-" * 76)
    print(f"  Current Regimen    : {vanco.dose_mg:.0f} mg q{vanco.interval_hours:.0f}h (Daily: {vanco.daily_dose_mg:.0f} mg)")
    print(f"  Steady-State AUC24 : {vanco.auc24_mg_h_l:.1f} mg*h/L (Ratio: {vanco.auc24_mic_ratio:.1f})")
    print(f"  Target Attained    : {vanco.target_achieved}")
    print(f"  Estimated Peak/Trg : Peak: {vanco.estimated_peak_mg_l} mg/L, Trough: {vanco.estimated_trough_mg_l} mg/L")
    print(f"  Risk Profile       : {vanco.nephrotoxicity_risk}")
    print(f"  Dose Recommendation: {vanco.recommended_regimen} (Target Daily: {vanco.recommended_daily_dose_mg:.0f} mg)")
    print("-" * 76)
    print(" 2. EXTENDED-INTERVAL AMINOGLYCOSIDE (Target Cmax/MIC >= 8-10):")
    print("-" * 76)
    print(f"  Regimen            : {genta.drug} {genta.dose_mg:.0f} mg ({genta.dose_per_kg} mg/kg) q{genta.interval_hours:.0f}h")
    print(f"  Cmax / MIC Ratio   : {genta.cmax_mic_ratio:.1f} (Peak: {genta.peak_cmax_mg_l} mg/L, Trough: {genta.trough_estimated_mg_l} mg/L)")
    print(f"  Target Achieved    : {genta.target_achieved} (Trough Safe: {genta.trough_safe})")
    print("-" * 76)
    print(" 3. BETA-LACTAM TIME-DEPENDENT TARGET (%fT>MIC):")
    print("-" * 76)
    print(f"  Regimen            : {mero.drug} {mero.dose_mg:.0f} mg q{mero.interval_hours:.0f}h")
    print(f"  %fT > MIC          : {mero.pct_ft_above_mic:.1f}% (Target: >={mero.target_threshold_pct}%)")
    print(f"  Continuous Infusion: {mero.continuous_infusion_rate_mg_hr} mg/hr (Loading Dose: {mero.loading_dose_mg:.0f} mg)")
    print("-" * 76)
    print(" 4. MONTE CARLO PROBABILITY OF TARGET ATTAINMENT (PTA, n=2,000):")
    print("-" * 76)
    for p in pta:
        print(f"  Regimen {p.regimen_label:14s} | PTA: {p.pta_target_pct:5.1f}% | Sub: {p.subtherapeutic_pct:4.1f}% | AKI: {p.nephrotoxic_pct:4.1f}% | Med AUC: {p.median_auc24:.0f}")
    print("=" * 76)


def interactive_cli() -> None:
    print("Antibiotic AUC/MIC Optimizer Interactive Shell. Type 'help' for commands, 'exit' to quit.\n")
    while True:
        try:
            line = input("auc-mic> ").strip()
        except (EOFError, KeyboardInterrupt):
            print("\nExiting.")
            break

        if not line:
            continue
        if line.lower() in ("exit", "quit"):
            break
        elif line.lower() == "help":
            print("Commands:")
            print("  demo                      - Run comprehensive clinical PK/PD dossier")
            print("  vanco <age> <M/F> <wt> <scr> <dose> <tau> [mic] - Optimize Vancomycin")
            print("  trapezoid <t1,c1> <t2,c2> ... - Compute trapezoidal AUC")
            print("  exit                      - Quit CLI")
        elif line.lower() == "demo":
            run_demo_evaluation(False)
        elif line.lower().startswith("vanco "):
            parts = line.split()[1:]
            if len(parts) < 6:
                print("Usage: vanco <age> <M/F> <wt_kg> <scr> <dose_mg> <interval_h> [mic]")
                continue
            try:
                p = PatientProfile("CLI", int(parts[0]), parts[1], float(parts[2]), 175.0, float(parts[3]))
                mic = float(parts[6]) if len(parts) > 6 else 1.0
                res = optimize_vancomycin(p, float(parts[4]), float(parts[5]), mic)
                print(f"AUC24: {res.auc24_mg_h_l} mg*h/L (Target: {res.target_achieved}) -> Rec: {res.recommended_regimen}")
            except Exception as ex:
                print(f"Error: {ex}")
        elif line.lower().startswith("trapezoid "):
            pairs_str = line.split()[1:]
            try:
                ts, cs = [], []
                for p in pairs_str:
                    t, c = p.split(",")
                    ts.append(float(t))
                    cs.append(float(c))
                auc = calculate_trapezoidal_auc(ts, cs)
                print(f"Computed Trapezoidal AUC: {auc} mg*h/L")
            except Exception as ex:
                print(f"Error parsing trapezoid points: {ex}")
        else:
            print(f"Unknown command: {line}. Type 'help'.")


def main() -> None:
    parser = argparse.ArgumentParser(description="Antibiotic AUC/MIC Optimizer & Precision Pharmacokinetics Engine")
    parser.add_argument("--demo", action="store_true", help="Run complete benchmark demonstration")
    parser.add_argument("--interactive", "-i", action="store_true", help="Launch interactive shell")
    parser.add_argument("--json", action="store_true", help="Output in JSON format")
    parser.add_argument("--batch", nargs=2, metavar=("INPUT_CSV", "OUTPUT_CSV"), help="Batch process patient CSV")
    parser.add_argument("--vanco", nargs="+", help="Vancomycin: <age> <M/F> <wt_kg> <scr> <dose_mg> <interval_h> [mic]")

    args = parser.parse_args()

    if args.batch:
        inp, out = args.batch
        res = process_csv(inp, out)
        print(f"Successfully processed {len(res)} cases from {inp} -> {out}")
    elif args.vanco:
        if len(args.vanco) < 6:
            print("Error: Required arguments: <age> <M/F> <wt_kg> <scr> <dose_mg> <interval_h> [mic]")
            sys.exit(1)
        age = int(args.vanco[0])
        gender = args.vanco[1]
        wt = float(args.vanco[2])
        scr = float(args.vanco[3])
        dose = float(args.vanco[4])
        interval = float(args.vanco[5])
        mic = float(args.vanco[6]) if len(args.vanco) > 6 else 1.0

        p = PatientProfile("CLI", age, gender, wt, 175.0, scr)
        res = optimize_vancomycin(p, dose, interval, mic)
        if args.json:
            print(json.dumps(res.to_dict(), indent=2))
        else:
            print(f"Vancomycin AUC24: {res.auc24_mg_h_l} mg*h/L (Target: 400-600)")
            print(f"Target Achieved : {res.target_achieved}")
            print(f"Risk Band       : {res.nephrotoxicity_risk}")
            print(f"Recommendation  : {res.recommended_regimen}")
    elif args.interactive:
        interactive_cli()
    else:
        run_demo_evaluation(json_out=args.json)


if __name__ == "__main__":
    main()
