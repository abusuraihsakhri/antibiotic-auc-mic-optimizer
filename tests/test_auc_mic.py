#!/usr/bin/env python3
"""
Comprehensive Unit Test Suite for Antibiotic AUC/MIC Optimizer & Pharmacokinetics Engine.
"""

import math
import os
import tempfile
import unittest
from auc_mic import (
    NephrotoxicityBand,
    PatientProfile,
    assess_row,
    calculate_cockcroft_gault,
    calculate_dosing_weight,
    calculate_ibw,
    calculate_kel_from_paired_levels,
    calculate_matzke_kel,
    calculate_trapezoidal_auc,
    optimize_aminoglycoside,
    optimize_beta_lactam,
    optimize_vancomycin,
    process_csv,
    run_monte_carlo_pta,
)


class TestRenalAndAnthropometrics(unittest.TestCase):
    def test_devine_ibw_male(self):
        # Male 6ft (72 inches) -> 50 + 2.3 * 12 = 77.6 kg
        ibw = calculate_ibw("Male", 182.88)
        self.assertAlmostEqual(ibw, 77.6, delta=0.1)

    def test_devine_ibw_female(self):
        # Female 5ft 4in (64 inches) -> 45.5 + 2.3 * 4 = 54.7 kg
        ibw = calculate_ibw("Female", 162.56)
        self.assertAlmostEqual(ibw, 54.7, delta=0.1)

    def test_dosing_weight_obesity(self):
        ibw = 60.0
        wt_obese = 100.0  # > 1.2 * 60 = 72
        # AdjBW = 60 + 0.4 * (100 - 60) = 76.0 kg
        dosing_wt, is_adj = calculate_dosing_weight(wt_obese, ibw)
        self.assertEqual(dosing_wt, 76.0)
        self.assertTrue(is_adj)

    def test_dosing_weight_normal(self):
        ibw = 60.0
        wt_norm = 65.0
        dosing_wt, is_adj = calculate_dosing_weight(wt_norm, ibw)
        self.assertEqual(dosing_wt, 65.0)
        self.assertFalse(is_adj)

    def test_cockcroft_gault_normal_male(self):
        res = calculate_cockcroft_gault(age_years=50, gender="Male", weight_kg=72.0, serum_creatinine=1.0)
        # (140-50)*72 / 72 = 90 mL/min
        self.assertAlmostEqual(res["crcl_ml_min"], 90.0, delta=0.5)

    def test_cockcroft_gault_female_multiplier(self):
        res = calculate_cockcroft_gault(age_years=50, gender="Female", weight_kg=72.0, serum_creatinine=1.0)
        # 90 * 0.85 = 76.5 mL/min
        self.assertAlmostEqual(res["crcl_ml_min"], 76.5, delta=0.5)

    def test_cockcroft_gault_invalid_inputs(self):
        with self.assertRaises(ValueError):
            calculate_cockcroft_gault(-10, "M", 70, 1.0)
        with self.assertRaises(ValueError):
            calculate_cockcroft_gault(50, "M", 70, 0.0)


class TestTrapezoidalAUC(unittest.TestCase):
    def test_linear_trapezoid_single_interval(self):
        # Time 0 to 2, C 10 to 10 -> Area = 10 * 2 = 20
        auc = calculate_trapezoidal_auc([0.0, 2.0], [10.0, 10.0])
        self.assertEqual(auc, 20.0)

    def test_log_linear_trapezoid_descending(self):
        # Time 0 to 4, C 40 to 10
        # AUC = (40 - 10) / ln(40/10) * 4 = 30 / 1.38629 * 4 = 86.56
        auc = calculate_trapezoidal_auc([0.0, 4.0], [40.0, 10.0])
        self.assertAlmostEqual(auc, 86.56, delta=0.1)

    def test_trapezoid_empty_or_single_point(self):
        self.assertEqual(calculate_trapezoidal_auc([], []), 0.0)
        self.assertEqual(calculate_trapezoidal_auc([1.0], [10.0]), 0.0)

    def test_trapezoid_mismatched_lengths_raises(self):
        with self.assertRaises(ValueError):
            calculate_trapezoidal_auc([1.0, 2.0], [10.0])


class TestVancomycinOptimizer(unittest.TestCase):
    def setUp(self):
        self.patient = PatientProfile("P_TEST", age_years=55, gender="Male", weight_kg=75.0, height_cm=175.0, serum_creatinine_mg_dl=1.0)

    def test_vanco_population_in_target_window(self):
        res = optimize_vancomycin(self.patient, dose_mg=1250.0, interval_hours=12.0, mic_mg_l=1.0)
        self.assertTrue(350.0 < res.auc24_mg_h_l < 650.0)
        self.assertEqual(res.dose_mg, 1250.0)
        self.assertTrue(res.clearance_l_hr > 0)

    def test_vanco_subtherapeutic_alert(self):
        # Low dose 500mg q24h
        res = optimize_vancomycin(self.patient, dose_mg=500.0, interval_hours=24.0, mic_mg_l=1.0)
        self.assertFalse(res.target_achieved)
        self.assertEqual(res.nephrotoxicity_risk, NephrotoxicityBand.SUBTHERAPEUTIC.value)
        self.assertIn("Subtherapeutic", res.guidance)

    def test_vanco_aki_high_risk_alert(self):
        # Massive overdose: 2500mg q8h = 7500mg/day
        res = optimize_vancomycin(self.patient, dose_mg=2500.0, interval_hours=8.0, mic_mg_l=1.0)
        self.assertFalse(res.target_achieved)
        self.assertEqual(res.nephrotoxicity_risk, NephrotoxicityBand.HIGH.value)

    def test_vanco_paired_levels_kel(self):
        # Peak 32, Trough 8, dt 12h -> kel = ln(32/8)/12 = ln(4)/12 = 1.38629/12 = 0.1155 /hr
        kel = calculate_kel_from_paired_levels(32.0, 8.0, 12.0)
        self.assertAlmostEqual(kel, 0.1155, delta=0.001)

        res = optimize_vancomycin(
            self.patient, dose_mg=1000.0, interval_hours=12.0, mic_mg_l=1.0,
            measured_peak_mg_l=32.0, measured_trough_mg_l=8.0, peak_trough_dt_hours=12.0
        )
        self.assertAlmostEqual(res.kel_per_hr, 0.1155, delta=0.001)

    def test_vanco_invalid_paired_levels_raises(self):
        with self.assertRaises(ValueError):
            calculate_kel_from_paired_levels(10.0, 20.0, 12.0)  # Peak < Trough
        with self.assertRaises(ValueError):
            calculate_kel_from_paired_levels(20.0, 10.0, 0.0)   # dt <= 0


class TestAminoglycosideOptimizer(unittest.TestCase):
    def test_gentamicin_hartford_dosing(self):
        patient = PatientProfile("P_GENT", age_years=45, gender="Male", weight_kg=70.0, height_cm=175.0, serum_creatinine_mg_dl=0.9)
        res = optimize_aminoglycoside(patient, drug_name="Gentamicin", dose_mg_per_kg=7.0, mic_mg_l=1.0)
        self.assertEqual(res.dose_mg, 490.0)  # 7 * 70
        self.assertEqual(res.interval_hours, 24.0)
        self.assertTrue(res.target_achieved)
        self.assertTrue(res.trough_safe)


class TestBetaLactamOptimizer(unittest.TestCase):
    def test_meropenem_ft_mic(self):
        patient = PatientProfile("P_MERO", age_years=50, gender="Male", weight_kg=80.0, height_cm=175.0, serum_creatinine_mg_dl=1.0)
        res = optimize_beta_lactam(patient, drug_name="Meropenem", dose_mg=1000.0, interval_hours=8.0, mic_mg_l=1.0)
        self.assertTrue(res.pct_ft_above_mic > 0)
        self.assertEqual(res.target_threshold_pct, 40.0)
        self.assertTrue(res.continuous_infusion_rate_mg_hr > 0)
        self.assertTrue(res.loading_dose_mg > 0)


class TestMonteCarloPTA(unittest.TestCase):
    def test_monte_carlo_pta_simulation(self):
        regimens = [
            (1000.0, 12.0, "1000mg q12h"),
            (1500.0, 12.0, "1500mg q12h")
        ]
        results = run_monte_carlo_pta(regimens, weight_kg=70.0, mean_crcl=80.0, mic_tested=1.0, n_simulations=1000, seed=123)
        self.assertEqual(len(results), 2)
        for r in results:
            self.assertEqual(r.n_simulated, 1000)
            self.assertTrue(0.0 <= r.pta_target_pct <= 100.0)
            self.assertTrue(r.mean_auc24 > 0)
            self.assertTrue(r.p5_auc24 < r.p95_auc24)


class TestCSVBatchAssessor(unittest.TestCase):
    def test_assess_row_valid(self):
        row = {"age": "60", "sex": "M", "weight": "75", "scr": "1.0", "dose": "1250", "interval": "12"}
        assessed = assess_row(row)
        self.assertIn("auc24_mg_h_l", assessed)
        self.assertIn("recommended_regimen", assessed)

    def test_process_csv_pipeline(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            inp = os.path.join(tmpdir, "input.csv")
            out = os.path.join(tmpdir, "output.csv")
            with open(inp, "w", newline="", encoding="utf-8") as f:
                f.write("age,sex,weight,scr,dose,interval\n")
                f.write("50,M,70,1.0,1000,12\n")
                f.write("65,F,60,1.5,1250,12\n")

            results = process_csv(inp, out)
            self.assertEqual(len(results), 2)
            self.assertTrue(os.path.exists(out))


class TestCLI(unittest.TestCase):
    def test_cli_demo_json(self):
        from cli import main
        import io
        from unittest.mock import patch
        with patch("sys.stdout", new=io.StringIO()) as fake_out:
            code = main(["--demo", "--json"])
            self.assertEqual(code, 0)
            self.assertIn("patient", fake_out.getvalue())

    def test_cli_vanco(self):
        from cli import main
        code = main(["--vanco", "60", "M", "75", "1.1", "1250", "12", "1.0"])
        self.assertEqual(code, 0)

    def test_cli_batch(self):
        from cli import main
        import tempfile
        sample_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "sample.csv")
        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = os.path.join(tmpdir, "out_batch.csv")
            code = main(["--batch", sample_path, out_file])
            self.assertEqual(code, 0)
            self.assertTrue(os.path.exists(out_file))


if __name__ == "__main__":
    unittest.main()


