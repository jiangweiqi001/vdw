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

    def test_tdhf_nstates_convergence_row_includes_errors_and_sum_rule(self):
        from run_tdhf_nstates_convergence_ar import summarize_nstates_result

        summary = summarize_nstates_result(
            basis="aug-cc-pVTZ",
            nstates=100,
            tdhf_summary={"sum_osc": 9.0, "n_channels": 3},
            comparison={
                "alpha0_eft": 10.0,
                "C6_eft": 60.0,
                "err_alpha_pct": -10.0,
                "err_C6_pct": -5.0,
            },
            n_electrons=18.0,
        )

        self.assertEqual(summary["basis"], "aug-cc-pVTZ")
        self.assertEqual(summary["nstates"], 100)
        self.assertAlmostEqual(summary["sum_osc"], 9.0)
        self.assertAlmostEqual(summary["sum_osc_over_N"], 0.5)
        self.assertAlmostEqual(summary["alpha0_err"], -10.0)
        self.assertAlmostEqual(summary["C6_err"], -5.0)

    def test_ar2_tail_row_uses_minus_c6_over_r6(self):
        from run_ar2_tail import tail_row

        row = tail_row(
            r_bohr=2.0,
            c6_ref=64.0,
            c6_tdhf=60.0,
            c6_mo=80.0,
            c6_calibrated=64.0,
        )

        self.assertAlmostEqual(row["E_ref_Ha"], -1.0)
        self.assertAlmostEqual(row["E_tdhf_Ha"], -60.0 / 64.0)
        self.assertAlmostEqual(row["err_tdhf_pct"], -6.25)
        self.assertAlmostEqual(row["err_mo_pct"], 25.0)

    def test_generic_tdhf_atom_channel_rows_from_arrays(self):
        from run_tdhf_atom import tdhf_rows_from_arrays

        rows = tdhf_rows_from_arrays(atom="Ne", energies=[0.4, 0.8], oscillator_strengths=[0.0, 1.2])

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Ne")
        self.assertEqual(rows[0]["channel"], "tdhf_002")
        self.assertAlmostEqual(rows[0]["delta_Ha"], 0.8)
        self.assertAlmostEqual(rows[0]["osc"], 1.2)

    def test_noble_gas_summary_row_uses_reference_errors(self):
        from run_noble_gas_tdhf import summarize_atom_result

        row = summarize_atom_result(
            atom="Ne",
            alpha0_ref=2.0,
            c6_ref=6.0,
            alpha_row={"alpha0_au": "1.8", "C6_self_au": "5.4", "n_channels": "3"},
        )

        self.assertEqual(row["atom"], "Ne")
        self.assertAlmostEqual(row["alpha0_err"], -10.0)
        self.assertAlmostEqual(row["C6_err"], -10.0)

    def test_ar_c6_method_rows_mark_unavailable_optional_methods(self):
        from run_ar_c6_method_comparison import method_row

        row = method_row("D4", None, 64.3, "not_available")

        self.assertEqual(row["method"], "D4")
        self.assertEqual(row["C6_ArAr"], "")
        self.assertEqual(row["error_pct"], "")
        self.assertEqual(row["note"], "not_available")

    def test_ne_convergence_summary_handles_unavailable_basis(self):
        from run_ne_tdhf_convergence import unavailable_row

        row = unavailable_row("aug-cc-pV5Z", 200, "missing basis")

        self.assertEqual(row["atom"], "Ne")
        self.assertEqual(row["basis"], "aug-cc-pV5Z")
        self.assertEqual(row["nstates"], 200)
        self.assertEqual(row["status"], "unavailable")
        self.assertIn("missing basis", row["note"])

    def test_ne_convergence_summary_row_has_errors(self):
        from run_ne_tdhf_convergence import summarize_case

        row = summarize_case(
            basis="aug-cc-pVQZ",
            nstates=300,
            alpha0=2.4,
            c6=5.8,
            n_channels=10,
            sum_osc=1.5,
            alpha0_ref=3.0,
            c6_ref=6.0,
        )

        self.assertEqual(row["status"], "ok")
        self.assertAlmostEqual(row["alpha0_err"], -20.0)
        self.assertAlmostEqual(row["C6_err"], -100.0 * 0.2 / 6.0)

    def test_tdhf_component_rows_can_label_non_ar_atom(self):
        from pyscf_export_ar_tdhf_decomposed import tdhf_component_rows

        rows = tdhf_component_rows(
            energies=[1.0],
            all_dipoles=[[2.0, 0.0, 0.0]],
            valence_dipoles=[[1.0, 0.0, 0.0]],
            core_dipoles=[[0.5, 0.0, 0.0]],
            atom="Kr",
            min_osc=0.0,
        )

        self.assertEqual(rows["all"][0]["atom"], "Kr")
        self.assertEqual(rows["all"][0]["component"], "all")
        self.assertAlmostEqual(rows["all"][0]["osc"], (2.0 / 3.0) * 4.0)
        self.assertAlmostEqual(rows["valence"][0]["osc"], (2.0 / 3.0) * 1.0)
        self.assertAlmostEqual(rows["core"][0]["osc"], (2.0 / 3.0) * 0.25)

    def test_partition_definition_parser_reads_semicolon_shells(self):
        from run_partition_decomposition import load_partition_definitions

        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "partition_definitions.csv"
            self._write_csv(
                path,
                ["atom", "partition", "shells"],
                [
                    {"atom": "Ca", "partition": "valence", "shells": "4s"},
                    {"atom": "Ca", "partition": "semicore", "shells": "3s;3p"},
                ],
            )

            partitions = load_partition_definitions(path, "Ca")

        self.assertEqual(partitions["valence"], {"4s"})
        self.assertEqual(partitions["semicore"], {"3s", "3p"})

    def test_partition_combinations_include_valence_plus_semicore(self):
        from run_partition_decomposition import partition_combinations

        combos = partition_combinations({"valence": {"4s"}, "semicore": {"3s", "3p"}, "deep_core": {"1s"}})

        self.assertEqual(combos["valence_plus_semicore"], {"4s", "3s", "3p"})
        self.assertEqual(combos["valence_plus_semicore_plus_deep_core"], {"4s", "3s", "3p", "1s"})

    def test_all_e_rpa_channel_rows_from_arrays(self):
        from run_all_e_rpa_atom import response_rows_from_arrays

        rows = response_rows_from_arrays("Ar", [0.2, 0.4], [0.0, 1.5], source="PySCF_TDDFT_LDA")

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Ar")
        self.assertEqual(rows[0]["channel"], "rpa_002")
        self.assertAlmostEqual(rows[0]["delta_Ha"], 0.4)
        self.assertAlmostEqual(rows[0]["osc"], 1.5)
        self.assertEqual(rows[0]["source"], "PySCF_TDDFT_LDA")

    def test_all_e_rpa_summary_row_uses_reference_errors(self):
        from run_all_e_rpa_atom import summarize_response

        row = summarize_response(
            atom="Ar",
            xc="lda",
            basis="aug-cc-pVTZ",
            nstates=20,
            method="TDDFT",
            alpha_row={"alpha0_au": "10.0", "C6_self_au": "60.0", "n_channels": "3"},
            alpha0_ref=11.0,
            c6_ref=66.0,
        )

        self.assertEqual(row["atom"], "Ar")
        self.assertEqual(row["method"], "TDDFT")
        self.assertAlmostEqual(row["alpha0_err"], -100.0 / 11.0)
        self.assertAlmostEqual(row["C6_err"], -100.0 / 11.0)

    def test_all_e_rpa_convergence_row_records_status(self):
        from run_all_e_rpa_convergence_ar import convergence_row

        row = convergence_row(
            basis="aug-cc-pVTZ",
            nstates=100,
            summary={
                "alpha0": 10.0,
                "C6": 60.0,
                "alpha0_err": -10.0,
                "C6_err": -5.0,
                "n_channels": 12,
            },
        )

        self.assertEqual(row["basis"], "aug-cc-pVTZ")
        self.assertEqual(row["nstates"], 100)
        self.assertEqual(row["status"], "ok")
        self.assertAlmostEqual(row["C6_err"], -5.0)

    def test_all_e_method_comparison_row_computes_deltas(self):
        from run_all_e_method_comparison import comparison_row

        row = comparison_row(
            atom="Ar",
            xc="pbe",
            ks_alpha0=11.7,
            ks_c6=68.0,
            hf_alpha0=10.6,
            hf_c6=60.7,
            ref_alpha0=11.1,
            ref_c6=64.3,
            basis="aug-cc-pVQZ",
        )

        self.assertEqual(row["atom"], "Ar")
        self.assertAlmostEqual(row["delta_C6_KS_minus_HF"], 7.3)
        self.assertAlmostEqual(row["KS_C6_err"], 100.0 * (68.0 - 64.3) / 64.3)

    def test_psp_rpa_summary_row_records_valence_metadata(self):
        from run_psp_rpa_atom import summarize_psp_response

        row = summarize_psp_response(
            atom="Ar",
            psp="gth-lda",
            basis="gth-dzvp",
            xc="lda",
            nstates=20,
            method="TDDFT",
            nelec=8,
            alpha_row={"alpha0_au": "7.0", "C6_self_au": "40.0", "n_channels": "12"},
        )

        self.assertEqual(row["atom"], "Ar")
        self.assertEqual(row["active_electrons"], 8)
        self.assertEqual(row["status"], "ok")
        self.assertAlmostEqual(row["C6_psp"], 40.0)

    def test_psp_rpa_summary_infers_mg_q2_active_shell(self):
        from run_psp_rpa_atom import summarize_psp_response

        row = summarize_psp_response(
            atom="Mg",
            psp="GTH-PBE-q2",
            basis="TZV2P-MOLOPT-SR-GTH-q2",
            xc="pbe",
            nstates=100,
            method="TDDFT",
            nelec=2,
            alpha_row={"alpha0_au": "72.0", "C6_self_au": "638.0", "n_channels": "6"},
        )

        self.assertEqual(row["active_shells"], "3s")

    def test_dipole_validation_uses_clean_mg_q2_path(self):
        from run_eft_core_dipole_validation import BEST_PSP

        self.assertIn("GTH-PBE-q2_TZV2P-MOLOPT-SR-GTH-q2_pbe_tddft", BEST_PSP["Mg"])
        self.assertNotIn("placeholder", BEST_PSP["Mg"])

    def test_psp_rpa_unavailable_row(self):
        from run_psp_rpa_atom import unavailable_row

        row = unavailable_row("Ca", "gth-lda", "gth-dzvp", "lda", 20, "missing basis")

        self.assertEqual(row["atom"], "Ca")
        self.assertEqual(row["status"], "unavailable")
        self.assertIn("missing basis", row["note"])

    def test_cp2k_named_block_extraction(self):
        from run_psp_rpa_atom import extract_cp2k_named_block

        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "BASIS"
            path.write_text(
                "# comment\n"
                "Ca SMALL-BASIS\n"
                "1\n"
                "2 0 0 1 1\n"
                "1.0 1.0\n"
                "Ca TARGET-BASIS TARGET-ALIAS\n"
                "1\n"
                "2 0 0 1 1\n"
                "2.0 1.0\n",
                encoding="utf-8",
            )

            block = extract_cp2k_named_block(path, "Ca", "TARGET-ALIAS")

        self.assertIn("Ca TARGET-BASIS TARGET-ALIAS", block)
        self.assertIn("2.0 1.0", block)

    def test_all_e_vs_psp_row_computes_missing_c6(self):
        from run_all_e_vs_psp_rpa_summary import comparison_row

        row = comparison_row(
            atom="Ar",
            all_e_c6=68.0,
            psp_c6=16.0,
            psp_status="ok",
            partition="gth-lda/gth-dzvp",
        )

        self.assertEqual(row["atom"], "Ar")
        self.assertAlmostEqual(row["delta_C6_missing"], 52.0)
        self.assertAlmostEqual(row["relative_missing_pct"], 100.0 * 52.0 / 68.0)

    def test_orbital_hartree_potential_and_self_coulomb(self):
        from compute_core_wilson import orbital_hartree_potential, self_coulomb

        r = np.array([1.0, 2.0, 3.0])
        u = np.array([1.0, 0.0, 0.0])

        vh = orbital_hartree_potential(r, u)
        j = self_coulomb(r, u)

        self.assertEqual(vh.shape, r.shape)
        self.assertGreaterEqual(vh[0], vh[-1])
        self.assertGreater(j, 0.0)

    def test_scalar_proxy_channel_uses_f0_over_delta_squared(self):
        from build_eft_core_correction import scalar_proxy_channel

        row = scalar_proxy_channel(
            {
                "atom": "Ca",
                "shell": "3s",
                "Delta_E_Ha": "2.0",
                "occupation": "2",
                "f0": "1.0",
            }
        )

        self.assertEqual(row["atom"], "Ca")
        self.assertEqual(row["channel"], "eft_core_scalar_proxy_3s")
        self.assertAlmostEqual(float(row["osc"]), 0.5)
        self.assertEqual(row["source"], "EFT_CORE_SCALAR_PROXY")

    def test_dipole_wilson_channels_from_arrays_selects_core_shells(self):
        from compute_dipole_wilson import dipole_wilson_channels_from_arrays

        mo_energy = [0.0, 1.0, 2.0]
        mo_occ = [2.0, 2.0, 0.0]
        dipole_mo = np.zeros((3, 3, 3))
        dipole_mo[0, 1, 2] = 2.0
        mo_to_shell = {0: "1s", 1: "2s", 2: "3p"}

        rows = dipole_wilson_channels_from_arrays(
            atom="Mg",
            mo_energy=mo_energy,
            mo_occ=mo_occ,
            dipole_mo=dipole_mo,
            mo_to_shell=mo_to_shell,
            selected_shells={"2s"},
        )

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["from_shell"], "2s")
        self.assertEqual(rows[0]["channel"], "eft_dipole_2s_occ_001_to_virt_002")
        self.assertAlmostEqual(rows[0]["delta_Ha"], 1.0)
        self.assertAlmostEqual(rows[0]["osc"], (2.0 / 3.0) * 1.0 * 2.0 * 4.0)
        self.assertAlmostEqual(rows[0]["occ_energy_Ha"], 1.0)
        self.assertAlmostEqual(rows[0]["virt_energy_Ha"], 2.0)
        self.assertEqual(rows[0]["source"], "EFT_CORE_DIPOLE_WILSON_MO_APPROX")

    def test_core_tdhf_wilson_rows_from_arrays(self):
        from compute_core_tdhf_wilson import core_tdhf_rows_from_arrays

        rows = core_tdhf_rows_from_arrays("Mg", [0.5, 1.0], [0.0, 2.0])

        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["atom"], "Mg")
        self.assertEqual(rows[0]["channel"], "core_tdhf_002")
        self.assertAlmostEqual(rows[0]["delta_Ha"], 1.0)
        self.assertAlmostEqual(rows[0]["osc"], 2.0)
        self.assertEqual(rows[0]["source"], "EFT_CORE_DIPOLE_WILSON_CORE_TDHF")

    def test_transition_density_dipole_reconstructs_oscillator_strength(self):
        from compute_multipole_core_wilson import oscillator_from_transition_dipole

        row = oscillator_from_transition_dipole(
            atom="Mg",
            state_index=1,
            energy=2.0,
            dipole_vector=np.array([3.0, 0.0, 0.0]),
            source="unit",
        )

        self.assertEqual(row["channel"], "multipole_tdhf_001")
        self.assertAlmostEqual(row["d2"], 9.0)
        self.assertAlmostEqual(row["osc"], (2.0 / 3.0) * 2.0 * 9.0)

    def test_transition_density_matrix_to_dipole_contracts_ov_integrals(self):
        from compute_multipole_core_wilson import transition_density_dipole

        dipole_ov = np.zeros((3, 1, 1))
        dipole_ov[2, 0, 0] = 4.0
        amplitudes = np.array([[0.5]])

        dipole = transition_density_dipole(dipole_ov, amplitudes)

        self.assertEqual(dipole.tolist(), [0.0, 0.0, 4.0])

    def test_double_counting_guard_blocks_overlapping_shells(self):
        from run_eft_core_dipole_validation import double_counting_status

        status = double_counting_status({"2s", "2p"}, {"2s", "2p", "3s"})

        self.assertEqual(status, "diagnostic_double_counting")

    def test_double_counting_guard_allows_clean_shells(self):
        from run_eft_core_dipole_validation import double_counting_status

        status = double_counting_status({"2s", "2p"}, {"3s"})

        self.assertEqual(status, "clean")

    def test_screened_pair_energy_recovers_bare_c6_tail(self):
        from screened_pairwise_vdw import screened_pair_energy

        self.assertAlmostEqual(screened_pair_energy(64.0, 2.0, model="bare"), -1.0)

    def test_dielectric_screening_scales_tail_by_epsilon_squared(self):
        from screened_pairwise_vdw import screened_pair_energy

        bare = screened_pair_energy(64.0, 2.0, model="bare")
        screened = screened_pair_energy(64.0, 2.0, model="dielectric", epsilon=2.0)

        self.assertAlmostEqual(screened, bare / 4.0)

    def test_logdet_vdw_second_order_matches_pair_tail_at_long_range(self):
        from screened_eft_vdw import logdet_vdw_energy

        channels = [
            {"atom": "A", "delta_Ha": "2.0", "osc": "8.0"},
            {"atom": "B", "delta_Ha": "2.0", "osc": "8.0"},
        ]
        positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 20.0]]
        energy = logdet_vdw_energy(channels, positions, atom_order=["A", "B"], n_quad=80, expansion_order=2)

        # Single oscillator has alpha0=2 and self C6=6 for identical atoms.
        self.assertAlmostEqual(energy, -6.0 / 20.0**6, delta=1e-8)

    def test_logdet_vdw_dielectric_reduces_pair_tail(self):
        from screened_eft_vdw import logdet_vdw_energy

        channels = [
            {"atom": "A", "delta_Ha": "2.0", "osc": "8.0"},
            {"atom": "B", "delta_Ha": "2.0", "osc": "8.0"},
        ]
        positions = [[0.0, 0.0, 0.0], [0.0, 0.0, 20.0]]
        bare = logdet_vdw_energy(channels, positions, atom_order=["A", "B"], n_quad=80, expansion_order=2)
        screened = logdet_vdw_energy(
            channels,
            positions,
            atom_order=["A", "B"],
            n_quad=80,
            expansion_order=2,
            screening={"model": "dielectric", "epsilon": 2.0},
        )

        self.assertAlmostEqual(screened, bare / 4.0, delta=1e-10)

    def test_mg_q2_benchmark_row_computes_closure(self):
        from run_mg_q2_benchmark import benchmark_row

        row = benchmark_row(
            psp_row={
                "atom": "Mg",
                "psp": "GTH-PBE-q2",
                "basis": "TZV2P-MOLOPT-SR-GTH-q2",
                "xc": "pbe",
                "method": "TDDFT",
                "active_electrons": "2",
                "active_shells": "3s",
                "C6_psp": "638.62015545",
                "alpha0_psp": "72.11681004",
            },
            all_e_row={
                "atom": "Mg",
                "xc": "pbe",
                "basis": "aug-cc-pVQZ",
                "method": "TDDFT",
                "C6": "647.58810039",
                "alpha0": "73.64242187",
            },
            eft_row={
                "atom": "Mg",
                "C6_psp": "638.62015545",
                "C6_psp_plus_dipole": "647.60794451",
                "Delta_C6_dipole": "8.98778906",
                "correction_shells": "2p;2s",
                "psp_explicit_valence_shells": "3s",
                "double_counting_status": "clean",
            },
        )

        self.assertEqual(row["benchmark_status"], "clean_candidate")
        self.assertAlmostEqual(row["closure_fraction"], 1.0022127834)
        self.assertAlmostEqual(row["residual_C6"], -0.01984412)

    def test_mg_q2_audit_flags_placeholder_paths(self):
        from run_mg_q2_benchmark import audit_row

        row = audit_row(
            psp_path="results/psp_rpa/mg/GTH-PBE-q2_TZV2P-MOLOPT-SR-GTH-q2_pbe_tddft/mg_psp_channels.csv",
            psp_row={"active_electrons": "2", "active_shells": "3s", "xc": "pbe", "method": "TDDFT"},
            all_e_row={"xc": "pbe", "method": "TDDFT"},
            eft_row={
                "correction_shells": "2p;2s",
                "psp_explicit_valence_shells": "3s",
                "double_counting_status": "clean",
            },
        )

        self.assertEqual(row["audit_status"], "pass")
        self.assertEqual(row["placeholder_path_used"], "false")
        self.assertEqual(row["shell_overlap"], "false")

    def test_mg_q2_sensitivity_rows_extract_shell_and_cutoff_cases(self):
        from run_mg_q2_benchmark import sensitivity_rows

        rows = sensitivity_rows(
            shell_rows=[
                {"atom": "Mg", "shell": "2s", "Delta_C6_shell": "0.33", "C6_psp": "638.62"},
                {"atom": "Ca", "shell": "3s", "Delta_C6_shell": "1.48", "C6_psp": "1972.56"},
            ],
            virtual_rows=[
                {"atom": "Mg", "max_delta_Ha": "5.0", "Delta_C6_cutoff": "8.63", "C6_psp": "638.62"},
            ],
        )

        self.assertEqual(len(rows), 2)
        self.assertEqual(rows[0]["case_type"], "shell")
        self.assertEqual(rows[0]["case"], "2s")
        self.assertEqual(rows[1]["case_type"], "virtual_cutoff")
        self.assertEqual(rows[1]["case"], "5.0")

    def test_mg_q2_stability_row_flags_tolerance(self):
        from run_mg_q2_benchmark import stability_row

        passing = stability_row(
            category="all_e_nstates",
            case="nstates_150",
            value=647.0,
            reference=648.0,
            tolerance_pct=0.2,
            quantity="C6_all_e",
        )
        review = stability_row(
            category="psp_nstates",
            case="nstates_20",
            value=630.0,
            reference=638.0,
            tolerance_pct=0.2,
            quantity="C6_psp",
        )

        self.assertEqual(passing["status"], "pass")
        self.assertEqual(review["status"], "review")
        self.assertAlmostEqual(passing["delta_pct"], -100.0 / 648.0)

    def test_mg_q2_selects_nstates_100_psp_baseline(self):
        from run_mg_q2_benchmark import select_mg_q2_psp_row

        row = select_mg_q2_psp_row([
            {"atom": "Mg", "psp": "GTH-PBE-q2", "basis": "TZV2P-MOLOPT-SR-GTH-q2", "xc": "pbe", "method": "TDDFT", "nstates": "20"},
            {"atom": "Mg", "psp": "GTH-PBE-q2", "basis": "TZV2P-MOLOPT-SR-GTH-q2", "xc": "pbe", "method": "TDDFT", "nstates": "100"},
        ])

        self.assertEqual(row["nstates"], "100")

    def test_mg_q2_selects_nstates_200_all_e_baseline(self):
        from run_mg_q2_benchmark import select_mg_pbe_all_e_row

        row = select_mg_pbe_all_e_row([
            {"atom": "Mg", "xc": "pbe", "basis": "aug-cc-pVQZ", "method": "TDDFT", "nstates": "100"},
            {"atom": "Mg", "xc": "pbe", "basis": "aug-cc-pVQZ", "method": "TDDFT", "nstates": "200"},
        ])

        self.assertEqual(row["nstates"], "200")

    def test_clean_candidate_specs_define_kr_q8_and_ca_q10_shells(self):
        from run_non_q2_clean_benchmarks import BENCHMARK_SPECS

        kr = BENCHMARK_SPECS["Kr_q8"]
        ca = BENCHMARK_SPECS["Ca_q10"]

        self.assertEqual(kr.atom, "Kr")
        self.assertEqual(kr.psp, "GTH-PBE-q8")
        self.assertEqual(kr.explicit_shells, {"4s", "4p"})
        self.assertEqual(kr.correction_shells, {"1s", "2s", "2p", "3s", "3p", "3d"})
        self.assertTrue(kr.correction_shells.isdisjoint(kr.explicit_shells))

        self.assertEqual(ca.atom, "Ca")
        self.assertEqual(ca.psp, "GTH-PBE-q10")
        self.assertEqual(ca.explicit_shells, {"3s", "3p", "4s"})
        self.assertEqual(ca.correction_shells, {"1s", "2s", "2p"})
        self.assertTrue(ca.correction_shells.isdisjoint(ca.explicit_shells))

    def test_ca_q2_pbe_adapted_spec_records_basis_provenance(self):
        from run_non_q2_clean_benchmarks import BENCHMARK_SPECS

        ca = BENCHMARK_SPECS["Ca_q2_PBE_adapted"]

        self.assertEqual(ca.atom, "Ca")
        self.assertEqual(ca.psp, "GTH-PBE-q2")
        self.assertEqual(ca.pseudo_file, "external_data/cp2k/POTENTIAL_UZH_CASR_Q2")
        self.assertEqual(ca.psp_basis, "TZV2P-MOLOPT-PBE-GTH-q10")
        self.assertEqual(ca.explicit_shells, {"4s"})
        self.assertEqual(ca.correction_shells, {"3s", "3p"})
        self.assertIn("adapted", ca.role)

    def test_be_q2_spec_uses_two_s_valence_and_one_s_correction(self):
        from run_non_q2_clean_benchmarks import BENCHMARK_SPECS
        from run_psp_rpa_atom import infer_active_shells

        be = BENCHMARK_SPECS["Be_q2"]

        self.assertEqual(be.atom, "Be")
        self.assertEqual(be.psp, "GTH-LDA-q2")
        self.assertEqual(be.psp_basis, "TZV2P-MOLOPT-PBE-GTH-q2")
        self.assertEqual(be.explicit_shells, {"2s"})
        self.assertEqual(be.correction_shells, {"1s"})
        self.assertEqual(infer_active_shells("Be", 2), "2s")

    def test_be_q2_lda_spec_uses_lda_backend_consistently(self):
        from run_non_q2_clean_benchmarks import BENCHMARK_SPECS

        be = BENCHMARK_SPECS["Be_q2_LDA"]

        self.assertEqual(be.atom, "Be")
        self.assertEqual(be.psp, "GTH-LDA-q2")
        self.assertEqual(be.psp_basis, "TZV2P-MOLOPT-PBE-GTH-q2")
        self.assertEqual(be.xc, "lda")
        self.assertEqual(be.method, "TDDFT")
        self.assertEqual(be.explicit_shells, {"2s"})
        self.assertEqual(be.correction_shells, {"1s"})

    def test_summary_row_without_external_reference_keeps_blank_ref_fields(self):
        from run_non_q2_clean_benchmarks import summarize_all_e_without_reference

        row = summarize_all_e_without_reference(
            atom="Be",
            xc="pbe",
            basis="aug-cc-pVQZ",
            nstates=100,
            method="TDDFT",
            alpha_row={"alpha0_au": "37.5", "C6_self_au": "210.0", "n_channels": "12"},
        )

        self.assertEqual(row["atom"], "Be")
        self.assertEqual(row["alpha0_ref"], "")
        self.assertEqual(row["C6_ref"], "")
        self.assertAlmostEqual(row["C6"], 210.0)

    def test_non_q2_benchmark_row_classifies_clean_and_closure(self):
        from run_non_q2_clean_benchmarks import BENCHMARK_SPECS, benchmark_row

        row = benchmark_row(
            spec=BENCHMARK_SPECS["Kr_q8"],
            psp_row={
                "atom": "Kr",
                "psp": "GTH-PBE-q8",
                "basis": "TZV2P-MOLOPT-PBE-GTH-q8",
                "xc": "pbe",
                "method": "TDDFT",
                "active_electrons": "8",
                "active_shells": "4s;4p",
                "C6_psp": "70.0",
            },
            all_e_row={
                "atom": "Kr",
                "xc": "pbe",
                "basis": "aug-cc-pVQZ",
                "method": "TDDFT",
                "C6": "100.0",
            },
            eft_row={
                "atom": "Kr",
                "C6_psp_plus_dipole": "85.0",
                "correction_shells": "1s;2s;2p;3s;3p;3d",
                "psp_explicit_valence_shells": "4s;4p",
                "double_counting_status": "clean",
            },
        )

        self.assertEqual(row["benchmark_status"], "clean_candidate")
        self.assertAlmostEqual(row["closure_fraction"], 0.5)
        self.assertAlmostEqual(row["residual_C6"], 15.0)

    def test_non_q2_audit_rejects_correction_overlap(self):
        from run_non_q2_clean_benchmarks import BENCHMARK_SPECS, audit_row

        row = audit_row(
            spec=BENCHMARK_SPECS["Ca_q10"],
            psp_path="results/psp_rpa/ca/GTH-PBE-q10_TZV2P-MOLOPT-PBE-GTH-q10_pbe_tddft/ca_psp_channels.csv",
            psp_row={"active_electrons": "10", "active_shells": "3s;3p;4s", "xc": "pbe", "method": "TDDFT"},
            all_e_row={"xc": "pbe", "method": "TDDFT"},
            eft_row={
                "correction_shells": "2p;3s",
                "psp_explicit_valence_shells": "3s;3p;4s",
                "double_counting_status": "diagnostic_double_counting",
            },
        )

        self.assertEqual(row["audit_status"], "fail")
        self.assertEqual(row["shell_overlap"], "true")

    def test_q2_candidate_status_requires_matched_basis(self):
        from probe_large_core_q2_candidates import candidate_status

        row = candidate_status(
            atom="Ca",
            pseudo_name="GTH-LDA-q2",
            basis_name="",
            can_build=False,
            can_run_rks=False,
            can_run_tddft=False,
            note="no basis",
        )

        self.assertEqual(row["candidate_status"], "no_matched_q2_basis")
        self.assertEqual(row["active_shells"], "4s")

    def test_q2_candidate_status_promotes_tddft_smoke(self):
        from probe_large_core_q2_candidates import candidate_status

        row = candidate_status(
            atom="Mg",
            pseudo_name="GTH-PBE-q2",
            basis_name="TZV2P-MOLOPT-SR-GTH-q2",
            can_build=True,
            can_run_rks=True,
            can_run_tddft=True,
            note="",
        )

        self.assertEqual(row["candidate_status"], "tddft_smoke_ok")
        self.assertEqual(row["active_electrons"], 2)

    def test_cp2k_header_discovery_finds_q2_aliases(self):
        from probe_large_core_q2_candidates import discover_headers

        headers = discover_headers(
            lines=[
                "# comment",
                "Mg GTH-PBE-q2",
                "2 4",
                "Ca GTH-PBE-q10 GTH-PBE",
                "2 4",
                "Mg TZV2P-MOLOPT-SR-GTH-q2",
                "1",
            ],
            atom="Mg",
            required_token="q2",
        )

        self.assertEqual(headers, ["GTH-PBE-q2", "TZV2P-MOLOPT-SR-GTH-q2"])

    def test_cp2k_header_discovery_does_not_treat_q20_as_q2(self):
        from probe_large_core_q2_candidates import discover_headers

        headers = discover_headers(
            lines=[
                "Cd QZVPP-MOLOPT-PBE-GTH-q20",
                "Cd GTH-LDA-q2",
            ],
            atom="Cd",
            required_token="q2",
        )

        self.assertEqual(headers, ["GTH-LDA-q2"])

    def test_q2_basis_rank_prefers_larger_basis(self):
        from probe_large_core_q2_candidates import preferred_bases

        ranked = preferred_bases([
            "DZVP-MOLOPT-SR-GTH-q2",
            "TZV2P-MOLOPT-SR-GTH-q2",
            "SZV-MOLOPT-SR-GTH-q2",
        ])

        self.assertEqual(ranked[0], "TZV2P-MOLOPT-SR-GTH-q2")

    def test_imported_q2_candidate_basis_records_provenance(self):
        from probe_large_core_q2_candidates import load_imported_basis_candidates

        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "candidates.csv"
            self._write_csv(
                path,
                ["atom", "basis_file", "basis_name", "basis_label", "basis_provenance", "basis_note"],
                [
                    {
                        "atom": "Ca",
                        "basis_file": "external_data/cp2k/BASIS_MOLOPT_UZH",
                        "basis_name": "TZV2P-MOLOPT-PBE-GTH-q10",
                        "basis_label": "TZV2P-MOLOPT-PBE-GTH-q10-as-q2-adapted",
                        "basis_provenance": "adapted_from_q10",
                        "basis_note": "q10 basis used as q2 candidate",
                    }
                ],
            )

            candidates = load_imported_basis_candidates(path, atom="Ca")

        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]["basis_name"], "TZV2P-MOLOPT-PBE-GTH-q10")
        self.assertEqual(candidates[0]["basis_label"], "TZV2P-MOLOPT-PBE-GTH-q10-as-q2-adapted")
        self.assertEqual(candidates[0]["basis_provenance"], "adapted_from_q10")

    def test_q2_candidate_status_includes_basis_provenance(self):
        from probe_large_core_q2_candidates import candidate_status

        row = candidate_status(
            atom="Ca",
            pseudo_name="GTH-LDA-q2",
            basis_name="TZV2P-MOLOPT-PBE-GTH-q10-as-q2-adapted",
            can_build=True,
            can_run_rks=True,
            can_run_tddft=True,
            note="smoke ok",
            basis_provenance="adapted_from_q10",
        )

        self.assertEqual(row["candidate_status"], "tddft_smoke_ok")
        self.assertEqual(row["basis_provenance"], "adapted_from_q10")

    @staticmethod
    def _write_csv(path, fieldnames, rows):
        with open(path, "w", newline="", encoding="utf-8") as fp:
            writer = csv.DictWriter(fp, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
