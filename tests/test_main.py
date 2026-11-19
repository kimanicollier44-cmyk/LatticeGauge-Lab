"""Tests for LatticeGauge-Lab. Every test grades a computed number
against a closed form or an independent implementation in main.py."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import main as lg

TOL = 1e-9


def test_jacobi_known_spectra():
    assert all(abs(a - b) < 1e-12 for a, b in
               zip(lg.hermitian_eigenvalues(lg.mat([[2, 1], [1, 2]])),
                   (1.0, 3.0)))
    assert all(abs(a - b) < 1e-12 for a, b in
               zip(lg.hermitian_eigenvalues(lg.mat([[0, 1], [1, 0]])),
                   (-1.0, 1.0)))


def test_plaquette_spectrum_J0_binomial():
    vals, _ = lg.plaquette_spectrum(1.0, 0.0)
    expect = sorted(-(4 - 2 * k) for k in range(5)
                    for _ in range(math.comb(4, k)))
    assert all(abs(a - b) < 1e-10 for a, b in zip(vals, expect))


def test_plaquette_spectrum_h0_two_levels():
    vals, _ = lg.plaquette_spectrum(0.0, 1.0)
    assert all(abs(a - b) < 1e-10 for a, b in zip(vals, [-1.0] * 8 + [1.0] * 8))


def test_plaquette_ground_closed_form():
    for J in (0.0, 0.5, 1.0, 2.0):
        vals, _ = lg.plaquette_spectrum(1.0, J)
        assert abs(vals[0] + math.sqrt(16.0 + J * J)) < 1e-9


def test_plaquette_trace_identities():
    vals, _ = lg.plaquette_spectrum(1.0, 1.0)
    assert abs(sum(vals)) < 1e-10
    assert abs(sum(v * v for v in vals) - 16 * 5) < 1e-8


def test_gauss_commutators_zero():
    for L, plaqs, verts in lg.LATTICES.values():
        for v in verts:
            assert lg.commutator_max_dev(L, plaqs, v, 1.0, 0.7) < 1e-12


def test_sector_dims_two_ways():
    assert len(lg.sector_states(lg.PLAQ_LINKS, lg.PLAQ_VERTICES)) == 2
    assert lg.sector_dimension_formula(lg.PLAQ_LINKS, lg.PLAQ_VERTICES) == 2
    assert len(lg.sector_states(lg.LADDER_LINKS, lg.LADDER_VERTICES)) == 4
    assert lg.sector_dimension_formula(
        lg.LADDER_LINKS, lg.LADDER_VERTICES) == 4


def test_ground_states_physical():
    _, g = lg.plaquette_spectrum(1.0, 1.0)
    assert all(abs(v - 1.0) < 1e-9
               for v in lg.gauss_expectations(g, lg.PLAQ_VERTICES))
    _, g = lg.ladder_ground(1.0, 1.0)
    assert all(abs(v - 1.0) < 1e-9
               for v in lg.gauss_expectations(g, lg.LADDER_VERTICES))


def test_lanczos_refereed_on_plaquette():
    e_lanczos, _ = lg.lanczos_ground(lg.PLAQ_LINKS, lg.PLAQ_PLAQUETTES,
                                     1.0, 1.0)
    vals, _ = lg.plaquette_spectrum(1.0, 1.0)
    assert abs(e_lanczos - vals[0]) < 1e-9


def test_ladder_lanczos_matches_sector_jacobi():
    e_lanczos, _ = lg.ladder_ground(1.0, 1.0)
    e_sector = lg.sector_ground_energy(lg.LADDER_LINKS, lg.LADDER_PLAQUETTES,
                                       lg.LADDER_VERTICES, 1.0, 1.0)
    assert abs(e_lanczos - e_sector) < 1e-8


def test_ladder_decoupled_limits():
    assert abs(lg.ladder_ground(1.0, 0.0)[0] + 7.0) < 1e-9
    assert abs(lg.ladder_ground(0.0, 1.0)[0] + 2.0) < 1e-9


def test_ladder_weak_series():
    e, _ = lg.ladder_ground(1.0, 0.1)
    assert abs(e - lg.ladder_weak_series(1.0, 0.1)) < 1e-5


def test_wilson_limits_and_monotonic():
    sweep = lg.wilson_sweep()
    w = {row["J_over_h"]: row["wilson_W"] for row in sweep.values()}
    assert abs(w[0.0]) < 1e-12
    assert w[20.0] > 0.98
    rs = lg.WILSON_SWEEP
    assert all(w[a] < w[b] for a, b in zip(rs, rs[1:]))


def test_wilson_perturbation_slope():
    assert abs(lg.wilson_slope(0.05) - 0.25) < 0.005


def test_wilson_hellmann_feynman():
    eps = 1e-4
    e_plus, _ = lg.plaquette_spectrum(1.0, 1.0 + eps)
    e_minus, _ = lg.plaquette_spectrum(1.0, 1.0 - eps)
    _, g = lg.plaquette_spectrum(1.0, 1.0)
    W, _, _ = lg.plaquette_expectations(g)
    assert abs((e_plus[0] - e_minus[0]) / (2 * eps) + W) < 1e-4


def test_thooft_dual_behaviour():
    sweep = lg.wilson_sweep()
    t = {row["J_over_h"]: row["thooft_X0"] for row in sweep.values()}
    assert abs(t[0.0] - 1.0) < 1e-12
    assert t[20.0] < 0.2
    rs = lg.WILSON_SWEEP
    assert all(t[a] > t[b] for a, b in zip(rs, rs[1:]))


def test_charge_costs_rise_with_separation():
    conf = lg.confinement_report()
    costs = [row["energy_cost"] for row in conf["pairs"].values()]
    assert all(v > 0 for v in costs)
    assert costs[0] < costs[1] < costs[2]


def test_charge_sector_referee():
    e_sector = lg.sector_ground_energy(
        lg.LADDER_LINKS, lg.LADDER_PLAQUETTES, lg.LADDER_VERTICES,
        1.0, 1.0, charges=frozenset((0, 3)))
    e_lanczos = lg.sector_ground_lanczos(
        lg.LADDER_LINKS, lg.LADDER_PLAQUETTES, lg.LADDER_VERTICES,
        1.0, 1.0, charges=frozenset((0, 3)))
    assert abs(e_sector - e_lanczos) < 1e-8


def test_classifier_scenarios():
    rep = lg.classifier_report()
    assert rep["deep_electric_r0.1"]["knn"] == "electric-dominated"
    assert rep["crossover_r1.0"]["knn"] == "crossover"
    assert rep["magnetic_r4.0"]["knn"] == "magnetic-dominated"
    assert rep["deep_magnetic_r10.0"]["knn"] == "magnetic-dominated"
    assert all(row["knn"] == row["policy"] for row in rep.values())


def test_classifier_grid_consistency():
    for feats, label in lg.TRAINING_GRID:
        assert lg.classify_features(feats) == label


def test_hold_gates():
    assert len(lg.HOLD_GATES) == 14
    assert lg.gate("qcd_or_continuum_claim") == "HOLD"
    assert lg.gate("confinement_proof_claim") == "HOLD"
    assert lg.gate("model_as_measurement_claim") == "HOLD"
    assert lg.gate("live_hardware_control") == "HOLD"
    assert lg.gate("investment_claim") == "HOLD"
    assert lg.gate("certification_claim") == "HOLD"
    assert lg.gate("run_simulation") == "PROCEED"


def test_seats_never_override():
    brains = lg.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["model_id"] == "gemini-argon-4"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = lg.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]


def test_cli_hold_demo_exit_2():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "hold-demo"],
        capture_output=True, text=True)
    assert proc.returncode == 2
    assert "HOLD" in proc.stdout


def test_cli_checks_command():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "checks"],
        capture_output=True, text=True)
    assert proc.returncode == 0
    assert "CHECKS " in proc.stdout
    assert "FAIL" not in proc.stdout


def test_demo_deterministic():
    main_py = str(Path(__file__).resolve().parent.parent / "main.py")
    a = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    b = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    assert a.stdout == b.stdout and a.returncode == 0
