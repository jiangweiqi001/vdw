import csv
import math
import sys
import unittest
from pathlib import Path
from tempfile import TemporaryDirectory

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

    @staticmethod
    def _write_csv(path, fieldnames, rows):
        with open(path, "w", newline="", encoding="utf-8") as fp:
            writer = csv.DictWriter(fp, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(rows)


if __name__ == "__main__":
    unittest.main()
