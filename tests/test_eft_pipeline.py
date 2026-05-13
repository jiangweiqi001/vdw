import csv
import math
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from eft_alpha import alpha0_from_osc, self_c6_from_osc


class EftPipelineTests(unittest.TestCase):
    def test_compare_alpha_c6_reports_percent_errors(self):
        from compare_alpha_c6 import compare_tables

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            eft_path = tmp_path / "alpha_c6_table.csv"
            ref_path = tmp_path / "reference_alpha_c6.csv"

            self._write_csv(
                eft_path,
                ["atom", "alpha0_au", "C6_self_au", "n_channels"],
                [{"atom": "Na", "alpha0_au": "110.0", "C6_self_au": "210.0", "n_channels": "1"}],
            )
            self._write_csv(
                ref_path,
                ["atom", "alpha0_ref", "C6_ref", "source"],
                [{"atom": "Na", "alpha0_ref": "100.0", "C6_ref": "200.0", "source": "unit"}],
            )

            rows = compare_tables(eft_path, ref_path)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Na")
        self.assertAlmostEqual(rows[0]["err_alpha_pct"], 10.0)
        self.assertAlmostEqual(rows[0]["err_C6_pct"], 5.0)

    def test_calibrated_channels_reproduce_reference_alpha0_and_self_c6(self):
        from build_eft_channels import build_calibrated_channels

        with TemporaryDirectory() as tmp:
            ref_path = Path(tmp) / "reference_alpha_c6.csv"
            self._write_csv(
                ref_path,
                ["atom", "alpha0_ref", "C6_ref", "source"],
                [{"atom": "Ar", "alpha0_ref": "11.1", "C6_ref": "64.3", "source": "unit"}],
            )

            rows = build_calibrated_channels(ref_path)

        self.assertEqual(len(rows), 1)
        row = rows[0]
        delta = float(row["delta_Ha"])
        osc = float(row["osc"])
        self.assertEqual(row["atom"], "Ar")
        self.assertEqual(row["channel"], "calibrated_single_oscillator")
        self.assertEqual(row["d2"], "")
        self.assertAlmostEqual(alpha0_from_osc([delta], [osc]), 11.1, places=6)
        self.assertAlmostEqual(self_c6_from_osc([delta], [osc]), 64.3, places=5)

    def test_spectral_builder_converts_orbital_specs_to_channel_rows(self):
        from build_eft_channels_spectral import build_spectral_channels

        rows = build_spectral_channels(
            atom="Na",
            occupied_core_orbitals=[
                {
                    "label": "2p",
                    "energy_Ha": -1.25,
                    "occupation": 6.0,
                    "is_frozen_core": True,
                }
            ],
            virtual_orbitals=[
                {
                    "label": "3s",
                    "energy_Ha": -0.25,
                    "dipole": [1.0, 2.0, 0.0],
                }
            ],
        )

        self.assertEqual(len(rows), 1)
        row = rows[0]
        self.assertEqual(row["atom"], "Na")
        self.assertEqual(row["channel"], "2p_to_3s")
        self.assertEqual(row["is_core"], "true")
        self.assertAlmostEqual(float(row["delta_Ha"]), 1.0)
        self.assertAlmostEqual(float(row["d2"]), 5.0)
        self.assertAlmostEqual(float(row["osc"]), (2.0 / 3.0) * 1.0 * 6.0 * 5.0)

    def test_loader_accepts_consistent_d2_and_osc_columns(self):
        from eft_alpha import load_channels_csv

        with TemporaryDirectory() as tmp:
            channels_path = Path(tmp) / "atomic_channels.csv"
            self._write_csv(
                channels_path,
                ["atom", "channel", "delta_Ha", "d2", "osc", "is_core"],
                [
                    {
                        "atom": "Na",
                        "channel": "2p_to_3s",
                        "delta_Ha": "1.0",
                        "d2": "5.0",
                        "osc": str((2.0 / 3.0) * 1.0 * 6.0 * 5.0),
                        "is_core": "true",
                    }
                ],
            )

            data = load_channels_csv(channels_path)

        self.assertEqual(data["Na"]["delta"].tolist(), [1.0])
        self.assertEqual(data["Na"]["osc"].tolist(), [20.0])

    def test_compare_pair_c6_reports_heteronuclear_percent_errors(self):
        from compare_pair_c6 import compare_pair_tables

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            model_path = tmp_path / "c6_table.csv"
            ref_path = tmp_path / "reference_pair_c6.csv"

            self._write_csv(
                model_path,
                ["A", "B", "C6_au"],
                [{"A": "Na", "B": "K", "C6_au": "2200.0"}],
            )
            self._write_csv(
                ref_path,
                ["A", "B", "C6_ref", "source"],
                [{"A": "K", "B": "Na", "C6_ref": "2000.0", "source": "unit"}],
            )

            rows = compare_pair_tables(model_path, ref_path)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["pair"], "K-Na")
        self.assertAlmostEqual(rows[0]["err_pct"], 10.0)

    def test_spectral_builder_reads_orbital_and_radial_csv_inputs(self):
        from build_eft_channels_spectral import build_spectral_channels_from_csv

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            orbitals_path = tmp_path / "atomic_spectral_input.csv"
            dipoles_path = tmp_path / "radial_dipoles.csv"
            self._write_csv(
                orbitals_path,
                ["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
                [
                    {"atom": "Na", "orbital": "2s", "type": "core", "n": "2", "l": "0", "occupation": "2", "energy_Ha": "-2.0"},
                    {"atom": "Na", "orbital": "2p", "type": "core", "n": "2", "l": "1", "occupation": "6", "energy_Ha": "-1.0"},
                    {"atom": "Na", "orbital": "3s", "type": "virtual", "n": "3", "l": "0", "occupation": "0", "energy_Ha": "0.2"},
                    {"atom": "Na", "orbital": "3p", "type": "virtual", "n": "3", "l": "1", "occupation": "0", "energy_Ha": "0.5"},
                ],
            )
            self._write_csv(
                dipoles_path,
                ["atom", "from_orbital", "to_orbital", "d2"],
                [
                    {"atom": "Na", "from_orbital": "2p", "to_orbital": "3s", "d2": "5.0"},
                    {"atom": "Na", "from_orbital": "2s", "to_orbital": "3p", "d2": "3.0"},
                ],
            )

            rows = build_spectral_channels_from_csv(orbitals_path, dipoles_path)

        by_channel = {row["channel"]: row for row in rows}
        self.assertEqual(sorted(by_channel), ["2p_to_3s", "2s_to_3p"])
        self.assertAlmostEqual(float(by_channel["2p_to_3s"]["delta_Ha"]), 1.2)
        self.assertAlmostEqual(float(by_channel["2p_to_3s"]["osc"]), (2.0 / 3.0) * 1.2 * 6.0 * 5.0)
        self.assertAlmostEqual(float(by_channel["2s_to_3p"]["delta_Ha"]), 2.5)
        self.assertAlmostEqual(float(by_channel["2s_to_3p"]["osc"]), (2.0 / 3.0) * 2.5 * 2.0 * 3.0)

    def test_sum_rule_summary_reports_core_strength_and_response(self):
        from check_sum_rules import summarize_sum_rules

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            orbitals_path = tmp_path / "atomic_spectral_input.csv"
            channels_path = tmp_path / "atomic_channels.csv"
            self._write_csv(
                orbitals_path,
                ["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
                [
                    {"atom": "Na", "orbital": "2s", "type": "core", "n": "2", "l": "0", "occupation": "2", "energy_Ha": "-2.0"},
                    {"atom": "Na", "orbital": "2p", "type": "core", "n": "2", "l": "1", "occupation": "6", "energy_Ha": "-1.0"},
                ],
            )
            self._write_csv(
                channels_path,
                ["atom", "channel", "delta_Ha", "d2", "osc", "is_core"],
                [{"atom": "Na", "channel": "2p_to_3s", "delta_Ha": "1.0", "d2": "5.0", "osc": "4.0", "is_core": "true"}],
            )

            rows = summarize_sum_rules(orbitals_path, channels_path)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Na")
        self.assertAlmostEqual(rows[0]["N_core"], 8.0)
        self.assertAlmostEqual(rows[0]["sum_osc"], 4.0)
        self.assertAlmostEqual(rows[0]["sum_osc_over_N_core"], 0.5)
        self.assertAlmostEqual(rows[0]["alpha0"], 4.0)

    def test_spectral_builder_adds_residual_oscillator_to_saturate_core_sum(self):
        from build_eft_channels_spectral import build_spectral_channels_from_csv

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            orbitals_path = tmp_path / "atomic_spectral_input.csv"
            dipoles_path = tmp_path / "radial_dipoles.csv"
            residual_path = tmp_path / "residual_oscillators.csv"
            self._write_csv(
                orbitals_path,
                ["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
                [
                    {"atom": "Na", "orbital": "2s", "type": "core", "n": "2", "l": "0", "occupation": "2", "energy_Ha": "-2.0"},
                    {"atom": "Na", "orbital": "3p", "type": "virtual", "n": "3", "l": "1", "occupation": "0", "energy_Ha": "-1.0"},
                ],
            )
            self._write_csv(
                dipoles_path,
                ["atom", "from_orbital", "to_orbital", "d2"],
                [{"atom": "Na", "from_orbital": "2s", "to_orbital": "3p", "d2": "0.3"}],
            )
            self._write_csv(
                residual_path,
                ["atom", "delta_missing_Ha"],
                [{"atom": "Na", "delta_missing_Ha": "1.5"}],
            )

            rows = build_spectral_channels_from_csv(
                orbitals_path,
                dipoles_path,
                add_residual_oscillator=True,
                residual_path=residual_path,
            )

        missing = [row for row in rows if row["channel"] == "missing_core_continuum"]
        self.assertEqual(len(missing), 1)
        self.assertEqual(missing[0]["atom"], "Na")
        self.assertEqual(missing[0]["delta_Ha"], "1.500000000000")
        self.assertAlmostEqual(sum(float(row["osc"]) for row in rows), 2.0)

    def test_compare_modes_formats_mode_error_rows(self):
        from compare_modes import format_mode_comparison_rows

        rows = format_mode_comparison_rows(
            "spectral",
            [
                {
                    "atom": "Na",
                    "alpha0_eft": 80.0,
                    "alpha0_ref": 100.0,
                    "err_alpha_pct": -20.0,
                    "C6_eft": 150.0,
                    "C6_ref": 200.0,
                    "err_C6_pct": -25.0,
                }
            ],
        )

        self.assertEqual(rows[0]["mode"], "spectral")
        self.assertEqual(rows[0]["atom"], "Na")
        self.assertAlmostEqual(rows[0]["err_alpha_pct"], -20.0)
        self.assertAlmostEqual(rows[0]["err_C6_pct"], -25.0)

    def test_fit_residual_deltas_matches_reference_alpha0(self):
        from fit_residual_oscillators import fit_residual_deltas

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            orbitals_path = tmp_path / "atomic_spectral_input.csv"
            dipoles_path = tmp_path / "radial_dipoles.csv"
            ref_path = tmp_path / "reference_alpha_c6.csv"
            self._write_csv(
                orbitals_path,
                ["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
                [
                    {"atom": "Na", "orbital": "2s", "type": "core", "n": "2", "l": "0", "occupation": "2", "energy_Ha": "-2.0"},
                    {"atom": "Na", "orbital": "3p", "type": "virtual", "n": "3", "l": "1", "occupation": "0", "energy_Ha": "-1.0"},
                ],
            )
            self._write_csv(
                dipoles_path,
                ["atom", "from_orbital", "to_orbital", "d2"],
                [{"atom": "Na", "from_orbital": "2s", "to_orbital": "3p", "d2": "0.75"}],
            )
            self._write_csv(
                ref_path,
                ["atom", "alpha0_ref", "C6_ref", "source"],
                [{"atom": "Na", "alpha0_ref": "2.0", "C6_ref": "1.0", "source": "unit"}],
            )

            rows = fit_residual_deltas(orbitals_path, dipoles_path, ref_path)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Na")
        self.assertAlmostEqual(rows[0]["delta_missing_Ha"], 1.0)

    def test_radial_orbitals_generate_shell_averaged_dipoles(self):
        from build_radial_dipoles import build_radial_dipoles_from_orbitals

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            orbitals_path = tmp_path / "atomic_spectral_input.csv"
            radial_path = tmp_path / "radial_orbitals.csv"
            self._write_csv(
                orbitals_path,
                ["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
                [
                    {"atom": "X", "orbital": "1s", "type": "core", "n": "1", "l": "0", "occupation": "2", "energy_Ha": "-1.0"},
                    {"atom": "X", "orbital": "2p", "type": "virtual", "n": "2", "l": "1", "occupation": "0", "energy_Ha": "0.5"},
                ],
            )
            self._write_csv(
                radial_path,
                ["atom", "orbital", "r_bohr", "u"],
                [
                    {"atom": "X", "orbital": "1s", "r_bohr": "0.0", "u": "1.0"},
                    {"atom": "X", "orbital": "1s", "r_bohr": "1.0", "u": "1.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "0.0", "u": "1.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "1.0", "u": "1.0"},
                ],
            )

            rows = build_radial_dipoles_from_orbitals(orbitals_path, radial_path, normalize=False)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "X")
        self.assertEqual(rows[0]["from_orbital"], "1s")
        self.assertEqual(rows[0]["to_orbital"], "2p")
        self.assertAlmostEqual(rows[0]["radial_integral"], 0.5)
        self.assertAlmostEqual(rows[0]["angular_factor"], 1.0)
        self.assertAlmostEqual(rows[0]["d2"], 0.25)

    def test_radial_orbital_dipoles_skip_forbidden_delta_l(self):
        from build_radial_dipoles import build_radial_dipoles_from_orbitals

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            orbitals_path = tmp_path / "atomic_spectral_input.csv"
            radial_path = tmp_path / "radial_orbitals.csv"
            self._write_csv(
                orbitals_path,
                ["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
                [
                    {"atom": "X", "orbital": "1s", "type": "core", "n": "1", "l": "0", "occupation": "2", "energy_Ha": "-1.0"},
                    {"atom": "X", "orbital": "2s", "type": "virtual", "n": "2", "l": "0", "occupation": "0", "energy_Ha": "0.5"},
                ],
            )
            self._write_csv(
                radial_path,
                ["atom", "orbital", "r_bohr", "u"],
                [
                    {"atom": "X", "orbital": "1s", "r_bohr": "0.0", "u": "1.0"},
                    {"atom": "X", "orbital": "1s", "r_bohr": "1.0", "u": "1.0"},
                    {"atom": "X", "orbital": "2s", "r_bohr": "0.0", "u": "1.0"},
                    {"atom": "X", "orbital": "2s", "r_bohr": "1.0", "u": "1.0"},
                ],
            )

            rows = build_radial_dipoles_from_orbitals(orbitals_path, radial_path)

        self.assertEqual(rows, [])

    def test_check_radial_orbitals_reports_norm_and_monotonic_grid(self):
        from check_radial_orbitals import check_radial_orbitals

        with TemporaryDirectory() as tmp:
            radial_path = Path(tmp) / "radial_orbitals.csv"
            self._write_csv(
                radial_path,
                ["atom", "orbital", "r_bohr", "u"],
                [
                    {"atom": "X", "orbital": "1s", "r_bohr": "0.0", "u": "1.0"},
                    {"atom": "X", "orbital": "1s", "r_bohr": "1.0", "u": "1.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "0.0", "u": "1.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "0.5", "u": "1.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "1.0", "u": "1.0"},
                ],
            )

            report = check_radial_orbitals(radial_path, norm_tol=1e-8)

        self.assertEqual(len(report["orbitals"]), 2)
        by_orbital = {row["orbital"]: row for row in report["orbitals"]}
        self.assertTrue(by_orbital["1s"]["monotonic_grid"])
        self.assertAlmostEqual(by_orbital["1s"]["norm"], 1.0)
        self.assertTrue(by_orbital["1s"]["norm_ok"])
        self.assertEqual(report["coverage_warnings"], [])

    def test_check_radial_orbitals_flags_bad_norm_grid_and_coverage(self):
        from check_radial_orbitals import check_radial_orbitals

        with TemporaryDirectory() as tmp:
            radial_path = Path(tmp) / "radial_orbitals.csv"
            self._write_csv(
                radial_path,
                ["atom", "orbital", "r_bohr", "u"],
                [
                    {"atom": "X", "orbital": "1s", "r_bohr": "0.0", "u": "2.0"},
                    {"atom": "X", "orbital": "1s", "r_bohr": "1.0", "u": "2.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "0.0", "u": "1.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "0.75", "u": "1.0"},
                    {"atom": "X", "orbital": "2p", "r_bohr": "0.5", "u": "1.0"},
                ],
            )

            report = check_radial_orbitals(radial_path, norm_tol=1e-8)

        by_orbital = {row["orbital"]: row for row in report["orbitals"]}
        self.assertFalse(by_orbital["1s"]["norm_ok"])
        self.assertFalse(by_orbital["2p"]["monotonic_grid"])
        self.assertEqual(len(report["coverage_warnings"]), 1)

    def test_convert_atomic_solver_output_multiplies_R_by_r(self):
        from convert_atomic_solver_output import convert_rows

        rows = [
            {"atom": "Ar", "orbital": "3p", "r_bohr": "0.5", "R": "4.0"},
            {"atom": "Ar", "orbital": "3p", "r_bohr": "1.0", "R": "3.0"},
        ]

        converted = convert_rows(rows, input_kind="R", normalize=False)

        self.assertEqual(converted[0]["atom"], "Ar")
        self.assertEqual(converted[0]["orbital"], "3p")
        self.assertAlmostEqual(converted[0]["u"], 2.0)
        self.assertAlmostEqual(converted[1]["u"], 3.0)

    def test_convert_atomic_solver_output_preserves_u_and_can_normalize(self):
        from convert_atomic_solver_output import convert_rows

        rows = [
            {"atom": "Ar", "orbital": "3p", "r_bohr": "0.0", "u": "2.0"},
            {"atom": "Ar", "orbital": "3p", "r_bohr": "1.0", "u": "2.0"},
        ]

        converted = convert_rows(rows, input_kind="u", normalize=True)

        self.assertAlmostEqual(converted[0]["u"], 1.0)
        self.assertAlmostEqual(converted[1]["u"], 1.0)

    def test_analyze_channels_reports_alpha_and_c6_contributions(self):
        from analyze_channels import analyze_channels

        with TemporaryDirectory() as tmp:
            channels_path = Path(tmp) / "atomic_channels.csv"
            self._write_csv(
                channels_path,
                ["atom", "channel", "delta_Ha", "d2", "osc", "is_core"],
                [
                    {"atom": "Ar", "channel": "low", "delta_Ha": "1.0", "d2": "", "osc": "2.0", "is_core": "true"},
                    {"atom": "Ar", "channel": "high", "delta_Ha": "2.0", "d2": "", "osc": "4.0", "is_core": "true"},
                ],
            )

            rows = analyze_channels(channels_path)

        by_channel = {row["channel"]: row for row in rows}
        self.assertAlmostEqual(by_channel["low"]["alpha0_contribution"], 2.0)
        self.assertAlmostEqual(by_channel["high"]["alpha0_contribution"], 1.0)
        self.assertAlmostEqual(by_channel["low"]["alpha0_fraction"], 2.0 / 3.0)
        self.assertAlmostEqual(by_channel["high"]["single_channel_c6"], 1.5)
        self.assertGreater(by_channel["low"]["single_channel_c6_fraction"], by_channel["high"]["single_channel_c6_fraction"])

    def test_ar_basis_convergence_summary_extracts_key_channels(self):
        from run_basis_convergence_ar import summarize_basis_result

        alpha_rows = [{"atom": "Ar", "alpha0_au": "8.0", "C6_self_au": "70.0", "n_channels": "3"}]
        channel_rows = [
            {
                "atom": "Ar",
                "channel": "3p_to_3d",
                "delta_Ha": 1.3,
                "osc": 10.0,
                "alpha0_contribution": 5.0,
                "alpha0_fraction": 0.625,
                "cross_inclusive_c6_fraction": 0.7,
                "is_residual": False,
            },
            {
                "atom": "Ar",
                "channel": "missing_core_continuum",
                "delta_Ha": 1.5,
                "osc": 2.0,
                "alpha0_contribution": 1.0,
                "alpha0_fraction": 0.125,
                "cross_inclusive_c6_fraction": 0.2,
                "is_residual": True,
            },
        ]

        summary = summarize_basis_result("cc-pVTZ", alpha_rows, channel_rows)

        self.assertEqual(summary["basis"], "cc-pVTZ")
        self.assertAlmostEqual(summary["alpha0"], 8.0)
        self.assertAlmostEqual(summary["C6"], 70.0)
        self.assertAlmostEqual(summary["sum_osc_discrete"], 10.0)
        self.assertAlmostEqual(summary["sum_osc_residual"], 2.0)
        self.assertAlmostEqual(summary["3p_to_3d_delta"], 1.3)
        self.assertAlmostEqual(summary["3p_to_3d_osc"], 10.0)
        self.assertAlmostEqual(summary["3p_to_3d_alpha_fraction"], 0.625)
        self.assertAlmostEqual(summary["residual_fraction"], 0.125)

    def test_sum_rules_by_shell_flags_large_shell_ratio(self):
        from check_sum_rules_by_shell import summarize_by_shell

        with TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            orbitals_path = tmp_path / "atomic_spectral_input.csv"
            channels_path = tmp_path / "atomic_channels.csv"
            self._write_csv(
                orbitals_path,
                ["atom", "orbital", "type", "n", "l", "occupation", "energy_Ha"],
                [
                    {"atom": "Ar", "orbital": "3p", "type": "core", "n": "3", "l": "1", "occupation": "6", "energy_Ha": "-0.5"},
                ],
            )
            self._write_csv(
                channels_path,
                ["atom", "channel", "delta_Ha", "d2", "osc", "is_core"],
                [
                    {"atom": "Ar", "channel": "3p_to_3d", "delta_Ha": "1.0", "d2": "", "osc": "7.2", "is_core": "true"},
                    {"atom": "Ar", "channel": "missing_core_continuum", "delta_Ha": "1.5", "d2": "", "osc": "0.5", "is_core": "true"},
                ],
            )

            rows = summarize_by_shell(channels_path, orbitals_path, warn_ratio=1.1)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["from_orbital"], "3p")
        self.assertAlmostEqual(rows[0]["occupation"], 6.0)
        self.assertAlmostEqual(rows[0]["sum_osc"], 7.2)
        self.assertAlmostEqual(rows[0]["ratio_to_occupation"], 1.2)
        self.assertTrue(rows[0]["warn"])

    def test_basis_summary_marks_invalid_prediction_for_trk_excess(self):
        from run_basis_convergence_ar import summarize_basis_result

        alpha_rows = [{"atom": "Ar", "alpha0_au": "8.0", "C6_self_au": "70.0", "n_channels": "1"}]
        channel_rows = [
            {
                "atom": "Ar",
                "channel": "3p_to_3d",
                "delta_Ha": 1.0,
                "osc": 20.0,
                "alpha0_fraction": 1.0,
                "cross_inclusive_c6_fraction": 1.0,
                "is_residual": False,
            }
        ]

        summary = summarize_basis_result("aug-cc-pVQZ", alpha_rows, channel_rows, n_core=18.0)

        self.assertTrue(summary["invalid_for_prediction"])
        self.assertIn("sum_osc_discrete/N_core", summary["invalid_reason"])

    def test_mo_oscillator_channels_from_arrays(self):
        from pyscf_export_ar_mo_oscillators import mo_oscillator_channels_from_arrays

        mo_energy = [0.0, 2.0]
        mo_occ = [2.0, 0.0]
        dipole_mo = np.zeros((3, 2, 2))
        dipole_mo[0, 0, 1] = 3.0

        rows = mo_oscillator_channels_from_arrays(mo_energy, mo_occ, dipole_mo, atom="Ar")

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Ar")
        self.assertEqual(rows[0]["channel"], "mo_occ_000_to_virt_001")
        self.assertAlmostEqual(rows[0]["delta_Ha"], 2.0)
        self.assertAlmostEqual(rows[0]["osc"], (2.0 / 3.0) * 2.0 * 2.0 * 9.0)
        self.assertEqual(rows[0]["source"], "PySCF_MO")

    def test_alpha_table_can_read_custom_input_path(self):
        from run_alpha_table import build_alpha_rows

        with TemporaryDirectory() as tmp:
            channels_path = Path(tmp) / "ar_mo_channels.csv"
            self._write_csv(
                channels_path,
                ["atom", "channel", "delta_Ha", "osc", "is_core", "source"],
                [{"atom": "Ar", "channel": "mo", "delta_Ha": "2.0", "osc": "8.0", "is_core": "true", "source": "unit"}],
            )

            rows = build_alpha_rows(channels_path)

        self.assertEqual(rows[0]["atom"], "Ar")
        self.assertEqual(rows[0]["n_channels"], 1)
        self.assertAlmostEqual(float(rows[0]["alpha0_au"]), 2.0)

    def test_tdhf_channel_rows_from_arrays(self):
        from pyscf_export_ar_tdhf_oscillators import tdhf_channel_rows_from_arrays

        rows = tdhf_channel_rows_from_arrays([0.5, 1.0], [0.0, 2.0], min_osc=1e-12)

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Ar")
        self.assertEqual(rows[0]["channel"], "tdhf_002")
        self.assertAlmostEqual(rows[0]["delta_Ha"], 1.0)
        self.assertAlmostEqual(rows[0]["osc"], 2.0)
        self.assertEqual(rows[0]["source"], "PySCF_TDHF")

    @staticmethod
    def _write_csv(path, fieldnames, rows):
        with open(path, "w", newline="", encoding="utf-8") as fp:
            writer = csv.DictWriter(fp, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
