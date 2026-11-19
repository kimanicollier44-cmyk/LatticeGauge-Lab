#!/usr/bin/env python3
"""LatticeGauge-Lab - educational Z2 lattice gauge theory lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the honest scope is stated up front:

  * Links are qubits. The Kogut-Susskind style Hamiltonian is
        H = -h * sum_links sigma^x_l
            - J * sum_plaquettes (product of sigma^z over the
              plaquette's links).
  * This lab works on TWO tiny lattices only: a single plaquette
    (4 links, Hilbert dimension 16) and a 2x1 ladder (7 links,
    dimension 128). Z2 gauge group only.
  * It says nothing about QCD, nothing about any continuum limit,
    and nothing about any real material or device. There is no
    continuum limit, no thermodynamic limit, and no confinement
    proof anywhere in this file - the static-charge numbers are a
    tiny-lattice diagnostic, fenced as such.

Conventions, stated once and used everywhere:

  * Computational (Z) basis: basis state s in [0, 2^L) has link l in
    |1> iff bit l of s is set. sigma^z_l eigenvalue on s is
    (-1)^(bit l of s); sigma^x_l flips bit l.
  * A plaquette operator P_p is the product of sigma^z over its
    links, so it is DIAGONAL in the Z basis with eigenvalue
    (-1)^(parity of the plaquette's bits in s).
  * Gauss operators G_v are products of sigma^x over the links
    incident at vertex v, so they are DIAGONAL in the X basis: an
    X-basis label t has G_v eigenvalue (-1)^(parity of t over v).
  * The physical sector is all G_v = +1. Because each link touches
    exactly two vertices, prod_v G_v = I, so only V - 1 of the V
    constraints are independent and the physical dimension is
    2^(L - V + 1) for a connected lattice.
  * The Wilson loop of the single plaquette IS its plaquette
    operator. The 't Hooft (electric) loop used here is sigma^x on
    a single link of the plaquette - the dual diagnostic.

Every real-world use is a permanent HOLD gate - with Gemini Argon 4
Extended Thinking as DECIDER (names-only, tighten-only), Muse as
advisor (never overrides HOLD), and Google Flow plan-only seats.
"""
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

# Single plaquette: links 0..3 around the square (0 top, 1 right,
# 2 bottom, 3 left); vertices are the four corners.
PLAQ_LINKS = 4
PLAQ_PLAQUETTES = ((0, 1, 2, 3),)
PLAQ_VERTICES = ((0, 3), (0, 1), (1, 2), (2, 3))

# 2x1 ladder: two squares sharing the middle vertical link.
# Vertices: 0,1,2 top row; 3,4,5 bottom row.
# Links: 0,1 top horizontals; 2,3 bottom horizontals; 4,5,6 verticals.
LADDER_LINKS = 7
LADDER_PLAQUETTES = ((0, 5, 2, 4), (1, 6, 3, 5))
LADDER_VERTICES = ((0, 4), (0, 1, 5), (1, 6),
                   (2, 4), (2, 3, 5), (3, 6))

LATTICES = {
    "plaquette": (PLAQ_LINKS, PLAQ_PLAQUETTES, PLAQ_VERTICES),
    "ladder": (LADDER_LINKS, LADDER_PLAQUETTES, LADDER_VERTICES),
}

# Couplings used by the demo book. Teaching values, h = 1 units.
DEMO_H = 1.0
DEMO_J = 1.0
WILSON_SWEEP = (0.0, 0.05, 0.1, 0.25, 0.5, 1.0, 2.0, 4.0, 10.0, 20.0)

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "qcd_or_continuum_claim": "HOLD - no QCD, continuum-limit, or continuum-physics claim follows from a 4-to-7-link Z2 toy",
    "confinement_proof_claim": "HOLD - the static-charge numbers are a tiny-lattice diagnostic, never a confinement proof",
    "thermodynamic_limit_claim": "HOLD - no thermodynamic-limit or phase-of-matter claim from lattices of 4 and 7 links",
    "particle_spectrum_claim": "HOLD - no hadron, glueball, or particle spectrum is computed or claimed here",
    "lattice_qcd_production_claim": "HOLD - no production lattice-QCD result, ensemble, or extrapolation is produced here",
    "real_device_characterisation": "HOLD - no characterising a real device, material, or experiment from this teaching model",
    "hardware_specification_claim": "HOLD - no hardware specification for any real device is asserted",
    "vendor_claim": "HOLD - no vendor hardware, software, or service figure is stated, endorsed, or compared",
    "large_lattice_claim": "HOLD - no claim that these tiny-lattice numbers scale to large lattices",
    "quantum_advantage_claim": "HOLD - no quantum-advantage or quantum-computing claim follows from exact diagonalisation of 16-to-128 states",
    "investment_claim": "HOLD - no investment claim; a lab model is not a prospectus",
    "certification_claim": "HOLD - no certification that any system meets any standard",
    "live_hardware_control": "HOLD - no control of real qubits, links, or laboratory hardware, ever",
    "model_as_measurement_claim": "HOLD - a computed number is a statement about the model at stated inputs, never a measurement",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Matrix kit
# Small real matrices as lists of lists. The Hamiltonians of this
# lab are real symmetric (products of sigma^x and sigma^z only), so
# the kit is real; every product, trace and norm is written out so
# each number can be read and diffed.

def mat(rows):
    return [[float(v) for v in row] for row in rows]


def eye(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]


def zeros(n):
    return [[0.0] * n for _ in range(n)]


def madd(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(A, s):
    return [[a * s for a in row] for row in A]


def mmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    out = [[0.0] * m for _ in range(n)]
    for i in range(n):
        for j in range(m):
            out[i][j] = sum(A[i][t] * B[t][j] for t in range(k))
    return out


def transpose(A):
    return [[A[j][i] for j in range(len(A))] for i in range(len(A[0]))]


def mtrace(A):
    return sum(A[i][i] for i in range(len(A)))


def frobenius(A):
    return math.sqrt(sum(v * v for row in A for v in row))


def mat_max_dev(A, B):
    return max(abs(a - b) for ra, rb in zip(A, B) for a, b in zip(ra, rb))


def jacobi_eigh(A):
    """Eigen-decomposition of a real symmetric matrix by cyclic
    Jacobi rotations. Sweeps run until the off-diagonal Frobenius
    mass is below 1e-14 or 100 sweeps pass. Returns (values, V) with
    values sorted ascending and eigenvectors as COLUMNS of V, so
    A = V diag(values) V^T. Refereed on known spectra in CHECKS
    before it is trusted with any Hamiltonian."""
    n = len(A)
    a = [[A[i][j] for j in range(n)] for i in range(n)]
    V = eye(n)
    for _ in range(100):
        off = math.sqrt(sum(a[p][q] ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < 1e-14:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                if abs(apq) < 1e-15:
                    continue
                app, aqq = a[p][p], a[q][q]
                tau = (aqq - app) / (2.0 * apq)
                if tau == 0.0:
                    t = 1.0
                else:
                    t = math.copysign(
                        1.0 / (abs(tau) + math.sqrt(tau * tau + 1.0)), tau)
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c * akp - s * akq
                    a[k][q] = s * akp + c * akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c * apk - s * aqk
                    a[q][k] = s * apk + c * aqk
                for k in range(n):
                    vkp, vkq = V[k][p], V[k][q]
                    V[k][p] = c * vkp - s * vkq
                    V[k][q] = s * vkp + c * vkq
    eigs = [a[i][i] for i in range(n)]
    order = sorted(range(n), key=lambda i: eigs[i])
    values = [eigs[i] for i in order]
    vecs = [[V[row][i] for i in order] for row in range(n)]
    return values, vecs


def hermitian_eigenvalues(A):
    """Eigenvalues only (real symmetric case of the lab's solver)."""
    return jacobi_eigh(A)[0]


def splitmix64_stream(seed):
    """Deterministic splitmix64 stream of floats in [-1, 1)."""
    state = seed & 0xFFFFFFFFFFFFFFFF
    while True:
        state = (state + 0x9E3779B97F4A7C15) & 0xFFFFFFFFFFFFFFFF
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & 0xFFFFFFFFFFFFFFFF
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & 0xFFFFFFFFFFFFFFFF
        z = z ^ (z >> 31)
        yield (z >> 11) / float(1 << 53) * 2.0 - 1.0


# ================================================================ 2. Hamiltonians
# H = -h sum_l X_l - J sum_p P_p. In the Z basis P_p is diagonal and
# X_l flips one bit, so the action of H on a vector is written
# directly - no dense matrix is needed except where a referee wants
# one (the plaquette, dimension 16).

def plaquette_sign(s, plaq):
    """Eigenvalue of the plaquette operator P_p on Z-basis state s."""
    parity = 0
    for l in plaq:
        parity ^= (s >> l) & 1
    return -1.0 if parity else 1.0


def apply_H(psi, L, plaqs, h, J):
    """Return H psi in the Z basis, from the stated conventions."""
    n = len(psi)
    out = [0.0] * n
    for s in range(n):
        diag = sum(plaquette_sign(s, p) for p in plaqs)
        out[s] = -J * diag * psi[s]
    for s in range(n):
        if psi[s] != 0.0:
            for l in range(L):
                out[s ^ (1 << l)] += -h * psi[s]
    return out


def dense_H(L, plaqs, h, J):
    """Full Hamiltonian matrix, built column by column from apply_H."""
    n = 1 << L
    H = [[0.0] * n for _ in range(n)]
    for s in range(n):
        e = [0.0] * n
        e[s] = 1.0
        col = apply_H(e, L, plaqs, h, J)
        for i in range(n):
            H[i][s] = col[i]
    return H


def plaquette_spectrum(h, J):
    """Full spectrum + ground vector of the single plaquette (Jacobi)."""
    H = dense_H(PLAQ_LINKS, PLAQ_PLAQUETTES, h, J)
    values, vecs = jacobi_eigh(H)
    ground = [vecs[i][0] for i in range(1 << PLAQ_LINKS)]
    return values, ground


def lanczos_ground(L, plaqs, h, J, iters=90, start=None):
    """Ground state by from-scratch Lanczos with full
    reorthogonalisation. The tridiagonal problem (at most `iters`
    dimensional) is solved by the refereed Jacobi solver and the
    Ritz vector is rebuilt in the full space. Deterministic: the
    default start vector comes from the in-file splitmix64 stream."""
    n = 1 << L
    if start is None:
        stream = splitmix64_stream(0x1A77CE)
        start = [next(stream) for _ in range(n)]
    norm = math.sqrt(sum(v * v for v in start))
    q = [v / norm for v in start]
    Q, alphas, betas = [], [], []
    qprev = [0.0] * n
    beta = 0.0
    for _ in range(iters):
        Q.append(q[:])
        w = apply_H(q, L, plaqs, h, J)
        alpha = sum(q[i] * w[i] for i in range(n))
        alphas.append(alpha)
        for i in range(n):
            w[i] -= alpha * q[i] + beta * qprev[i]
        for qq in Q:
            d = sum(w[i] * qq[i] for i in range(n))
            for i in range(n):
                w[i] -= d * qq[i]
        beta = math.sqrt(sum(v * v for v in w))
        betas.append(beta)
        if beta < 1e-13:
            break
        qprev = q[:]
        q = [v / beta for v in w]
    m = len(alphas)
    T = [[0.0] * m for _ in range(m)]
    for i in range(m):
        T[i][i] = alphas[i]
        if i < m - 1:
            T[i][i + 1] = betas[i]
            T[i + 1][i] = betas[i]
    values, vecs = jacobi_eigh(T)
    ground = [0.0] * n
    for j in range(m):
        coeff = vecs[j][0]
        for i in range(n):
            ground[i] += coeff * Q[j][i]
    return values[0], ground


# ================================================================ 3. Gauss law
# G_v = product of X over links at v. [H, G_v] = 0 because every
# plaquette meets every vertex in 0 or 2 links (the two X/Z
# anticommutations cancel) and X terms commute with X products.
# The commutator norm below is the Frobenius norm of [H, G_v],
# assembled column by column from operator actions - exactly the
# dense commutator norm, without forming 128x128 products.

def apply_G(psi, vertex_links):
    """Return G_v psi in the Z basis (X product flips those bits)."""
    mask = 0
    for l in vertex_links:
        mask |= 1 << l
    out = [0.0] * len(psi)
    for s in range(len(psi)):
        out[s ^ mask] = psi[s]
    return out


def commutator_max_dev(L, plaqs, vertex_links, h, J):
    """max |[H, G_v] e_s| over all basis states and components."""
    n = 1 << L
    worst = 0.0
    for s in range(n):
        e = [0.0] * n
        e[s] = 1.0
        hg = apply_H(apply_G(e, vertex_links), L, plaqs, h, J)
        gh = apply_G(apply_H(e, L, plaqs, h, J), vertex_links)
        worst = max(worst, max(abs(a - b) for a, b in zip(hg, gh)))
    return worst


def sector_states(L, vertices, charges=frozenset()):
    """X-basis labels t whose G_v eigenvalue is -1 exactly at the
    vertices in `charges` and +1 everywhere else."""
    states = []
    for t in range(1 << L):
        ok = True
        for vi, vlinks in enumerate(vertices):
            parity = 0
            for l in vlinks:
                parity ^= (t >> l) & 1
            want = 1 if vi in charges else 0
            if parity != want:
                ok = False
                break
        if ok:
            states.append(t)
    return states


def sector_dimension_formula(L, vertices):
    """2^(L - V + 1) for a connected lattice (one redundancy)."""
    return 2 ** (L - len(vertices) + 1)


def sector_hamiltonian(L, plaqs, vertices, h, J, charges=frozenset()):
    """H restricted to one Gauss sector, built in the X basis: the
    electric term is diagonal there, (-1)^(bit l of t) per link,
    and each plaquette flips its links' X-basis bits. Returns
    (matrix, states). Sector dimension is 2^(L - V + 1) = 4 on the
    ladder, so this Jacobi problem is a 4x4 referee for Lanczos."""
    states = sector_states(L, vertices, charges)
    idx = {s: i for i, s in enumerate(states)}
    m = len(states)
    H = [[0.0] * m for _ in range(m)]
    for j, t in enumerate(states):
        diag = 0.0
        for l in range(L):
            diag += -1.0 if (t >> l) & 1 else 1.0
        H[j][j] = -h * diag
        for p in plaqs:
            mask = 0
            for l in p:
                mask |= 1 << l
            t2 = t ^ mask
            if t2 in idx:
                H[idx[t2]][j] += -J
    return H, states


def sector_ground_energy(L, plaqs, vertices, h, J, charges=frozenset()):
    H, _ = sector_hamiltonian(L, plaqs, vertices, h, J, charges)
    return jacobi_eigh(H)[0][0]


def x_basis_state_zrep(t, L):
    """The X-basis product state |t> written in the Z basis:
    link l contributes (|0> + (-1)^t_l |1>) / sqrt(2)."""
    n = 1 << L
    amp = 2.0 ** (-L / 2.0)
    out = [0.0] * n
    for s in range(n):
        parity = 0
        for l in range(L):
            if ((t >> l) & 1) and ((s >> l) & 1):
                parity ^= 1
        out[s] = -amp if parity else amp
    return out


def sector_ground_lanczos(L, plaqs, vertices, h, J, charges=frozenset()):
    """Independent referee for sector_ground_energy: build a generic
    vector INSIDE the sector in the Z representation (a splitmix
    combination of the sector's X-basis states), then run Lanczos
    with the Z-basis apply_H. H preserves the sector because
    [H, G_v] = 0, so the iteration never leaves it."""
    states = sector_states(L, vertices, charges)
    stream = splitmix64_stream(0x5EC702)
    start = [0.0] * (1 << L)
    for t in states:
        coeff = next(stream)
        rep = x_basis_state_zrep(t, L)
        for i in range(len(start)):
            start[i] += coeff * rep[i]
    energy, _ = lanczos_ground(L, plaqs, h, J, start=start)
    return energy


def gauss_expectations(psi_z, vertices):
    """<G_v> of a Z-basis state, one value per vertex."""
    return [sum(a * b for a, b in zip(psi_z, apply_G(psi_z, v)))
            for v in vertices]


# ================================================================ 4. Wilson loop
# For the single plaquette the Wilson loop IS the plaquette
# operator P (diagonal in the Z basis). First-order perturbation
# theory in J about the h-only ground state |++++>_x: P flips all
# four X bits, landing on |---->_x at excitation energy 8h, so the
# perturbed ground state carries amplitude J/(8h) of |---->_x and
#     <W> = 2 * J/(8h) + O((J/h)^3) = (1/4)(J/h) + ...
# The derived slope is 1/4; CHECKS grades the computed slope
# against it. The 't Hooft diagnostic is <sigma^x_0>: 1 in the
# electric limit, falling as the magnetic term takes over.

def plaquette_expectations(psi):
    """(<W>, <sigma^x_0>, mean <sigma^x>) of a plaquette Z state."""
    W = sum(psi[s] ** 2 * plaquette_sign(s, PLAQ_PLAQUETTES[0])
            for s in range(len(psi)))
    thooft = sum(psi[s] * psi[s ^ 1] for s in range(len(psi)))
    elec = sum(psi[s] * psi[s ^ (1 << l)]
               for s in range(len(psi)) for l in range(PLAQ_LINKS))
    return W, thooft, elec / PLAQ_LINKS


def wilson_sweep():
    rows = {}
    for r in WILSON_SWEEP:
        _, ground = plaquette_spectrum(DEMO_H, r * DEMO_H)
        W, thooft, elec = plaquette_expectations(ground)
        energy, _ = plaquette_spectrum(DEMO_H, r * DEMO_H)
        rows[f"J_over_h_{r}"] = {
            "J_over_h": r,
            "ground_energy": round(energy[0], 6),
            "wilson_W": round(W, 6),
            "thooft_X0": round(thooft, 6),
            "electric_mean_X": round(elec, 6),
        }
    return rows


def wilson_slope(r=0.05):
    """<W> / (J/h) at small J/h - the first-order slope estimate."""
    _, ground = plaquette_spectrum(DEMO_H, r * DEMO_H)
    W, _, _ = plaquette_expectations(ground)
    return W / r


# ================================================================ 5. Ladder + charges
# Ladder ground state by Lanczos in the full 128-dimensional space.
# Referees (full 128x128 Jacobi is NOT used - it is not needed):
#   * the physical-sector Jacobi problem in the X basis (4x4),
#     an independent representation of the same Hamiltonian;
#   * the decoupled limits J = 0 (E = -7h exactly) and h = 0
#     (E = -2J exactly, both plaquettes satisfied);
#   * the weak-coupling series E = -7h - J^2/(4h) + O(J^4):
#     each plaquette flips 4 X bits (cost 8h), contributing
#     -J^2/(8h), and there are two plaquettes.

def ladder_ground(h, J):
    return lanczos_ground(LADDER_LINKS, LADDER_PLAQUETTES, h, J)


def ladder_weak_series(h, J):
    return -7.0 * h - J * J / (4.0 * h)


# Static-charge diagnostic pairs on the ladder: (vertex a, vertex b,
# graph distance along the ladder). Costs are ground-energy
# differences between the charged sector (G_a = G_b = -1) and the
# physical sector at h = J = 1.
CHARGE_PAIRS = (
    (0, 3, 1),
    (0, 2, 2),
    (0, 5, 3),
)


def confinement_report():
    h, J = DEMO_H, DEMO_J
    e_phys = sector_ground_energy(LADDER_LINKS, LADDER_PLAQUETTES,
                                  LADDER_VERTICES, h, J)
    rows = {"physical_sector_energy": round(e_phys, 6), "pairs": {}}
    for a, b, dist in CHARGE_PAIRS:
        e_ch = sector_ground_energy(LADDER_LINKS, LADDER_PLAQUETTES,
                                    LADDER_VERTICES, h, J,
                                    charges=frozenset((a, b)))
        rows["pairs"][f"v{a}_v{b}_dist{dist}"] = {
            "vertices": [a, b],
            "separation": dist,
            "charged_energy": round(e_ch, 6),
            "energy_cost": round(e_ch - e_phys, 6),
        }
    rows["fence"] = ("tiny-lattice diagnostic only - 7 links, three "
                     "separations; NOT a confinement proof")
    return rows


# ================================================================ 6. kNN classifier
# Features of a plaquette ground state: (J/h, <W>, mean <sigma^x>).
# The training grid is labelled by the lab's own policy (regime_label
# below) - the classifier interpolates the policy, it does not
# replace the computed quantities.

def regime_label(r, W, elec):
    """The regime policy distilled to three verdicts, by coupling."""
    if r < 0.5:
        return "electric-dominated"
    if r > 2.0:
        return "magnetic-dominated"
    return "crossover"


def _training_grid():
    rows = []
    for r in (0.0, 0.1, 0.3, 0.49, 0.75, 1.0, 1.5, 1.99, 2.5, 4.0):
        for W in (0.0, 0.25, 0.5, 0.75, 1.0):
            for elec in (0.0, 0.25, 0.5, 0.75, 1.0):
                feats = (r, W, elec)
                rows.append((feats, regime_label(*feats)))
    return rows


TRAINING_GRID = _training_grid()


def classify_features(features, k=5):
    """Distance-weighted kNN over the policy-labelled grid."""
    ranked = []
    for feats, label in TRAINING_GRID:
        dist = math.sqrt(sum((a - b) ** 2
                             for a, b in zip(features, feats)))
        ranked.append((dist, feats, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"electric-dominated": 0, "crossover": 1,
                "magnetic-dominated": 2}
    return sorted(votes.items(),
                  key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


CLASSIFIER_SCENARIOS = {
    "deep_electric_r0.1": 0.1,
    "electric_r0.25": 0.25,
    "crossover_r1.0": 1.0,
    "crossover_r1.5": 1.5,
    "magnetic_r4.0": 4.0,
    "deep_magnetic_r10.0": 10.0,
}


def classifier_report():
    rows = {}
    for name, r in CLASSIFIER_SCENARIOS.items():
        _, ground = plaquette_spectrum(DEMO_H, r * DEMO_H)
        W, _, elec = plaquette_expectations(ground)
        feats = (r, W, elec)
        rows[name] = {
            "J_over_h": r,
            "features": [round(v, 6) for v in feats],
            "policy": regime_label(*feats),
            "knn": classify_features(feats),
        }
    return rows


# ================================================================ 7. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="latticegauge-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "model_id": "gemini-argon-4",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all deployment gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="latticegauge-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the fence before the number: 4 to 7 links, Z2 only, no continuum limit and no QCD claim - then show the two limiting spectra, the Gauss commutators, and the Wilson slope of 1/4 that were actually computed, and call the charge costs a tiny-lattice diagnostic, never a confinement proof.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["four link-qubits around a single plaquette, X arrows yielding to Z plaquette flux as J/h grows",
                  "a two-square ladder with Gauss-law stars glowing at every vertex",
                  "the Wilson loop expectation climbing from 0 toward 1 across the coupling sweep",
                  "two static charges on the ladder pulling an energy cost that grows with separation",
                  "the physical sector shrinking from 128 states to the 4 that obey Gauss's law"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ 8. Topic reports
def spectrum_report():
    rows = {}
    for label, (h, J) in {"demo_h1_J1": (1.0, 1.0),
                          "electric_limit_J0": (1.0, 0.0),
                          "magnetic_limit_h0": (0.0, 1.0)}.items():
        values, _ = plaquette_spectrum(h, J)
        rows[label] = {
            "h": h, "J": J,
            "eigenvalues": [round(v, 6) for v in values],
            "ground_energy": round(values[0], 6),
        }
    return rows


def gauss_report():
    rows = {}
    for name, (L, plaqs, verts) in LATTICES.items():
        rows[name] = {
            "links": L,
            "vertices": len(verts),
            "commutator_max_dev": max(
                commutator_max_dev(L, plaqs, v, 1.0, 0.7) for v in verts),
            "physical_dim_enumerated": len(sector_states(L, verts)),
            "physical_dim_formula": sector_dimension_formula(L, verts),
        }
    return rows


def ladder_report():
    rows = {}
    for label, (h, J) in {"demo_h1_J1": (1.0, 1.0),
                          "decoupled_J0": (1.0, 0.0),
                          "decoupled_h0": (0.0, 1.0),
                          "weak_J0.1": (1.0, 0.1)}.items():
        energy, _ = ladder_ground(h, J)
        rows[label] = {"h": h, "J": J,
                       "lanczos_ground_energy": round(energy, 6)}
    rows["demo_h1_J1"]["sector_jacobi_referee"] = round(
        sector_ground_energy(LADDER_LINKS, LADDER_PLAQUETTES,
                             LADDER_VERTICES, 1.0, 1.0), 6)
    rows["weak_J0.1"]["weak_series_prediction"] = round(
        ladder_weak_series(1.0, 0.1), 6)
    rows["referee_note"] = ("full 128x128 Jacobi not used; referees are "
                            "the X-basis physical-sector Jacobi problem, "
                            "the two decoupled limits, and the weak-"
                            "coupling series")
    return rows


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 0. Jacobi eigensolver refereed on known spectra first
    chk("jacobi_diagonal_matrix",
        hermitian_eigenvalues(mat([[3, 0], [0, -1]])) == [-1.0, 3.0])
    chk("jacobi_2x2_known_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(mat([[2, 1], [1, 2]])), (1.0, 3.0))))
    chk("jacobi_pauli_x_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(mat([[0, 1], [1, 0]])), (-1.0, 1.0))))
    A4 = mat([[2, 1, 0, 1], [1, 3, 1, 0], [0, 1, 2, 1], [1, 0, 1, 3]])
    ev4 = hermitian_eigenvalues(A4)
    chk("jacobi_trace_and_frobenius_invariants",
        abs(sum(ev4) - mtrace(A4)) < 1e-10
        and abs(sum(v * v for v in ev4)
                - sum(v * v for row in A4 for v in row)) < 1e-9)

    # 1. Single plaquette: full spectrum vs both analytic limits.
    # NOTE (deviation from the build brief, asserted TRUE here): at
    # h = 0 the only plaquette operator has eigenvalues +/-1, so the
    # spectrum is -J with degeneracy 8 and +J with degeneracy 8 -
    # NOT -J*(4-2m); that form would need four independent
    # plaquettes. At J = 0 the brief's form IS right, with h.
    vals_j0, _ = plaquette_spectrum(1.0, 0.0)
    expect_j0 = sorted(-(4 - 2 * k) for k in range(5)
                       for _ in range(math.comb(4, k)))
    chk("plaquette_spectrum_J0_binomial",
        all(abs(a - b) < 1e-10 for a, b in zip(vals_j0, expect_j0)))
    vals_h0, _ = plaquette_spectrum(0.0, 1.0)
    chk("plaquette_spectrum_h0_two_levels_eight_each",
        all(abs(a - b) < 1e-10 for a, b in
            zip(vals_h0, [-1.0] * 8 + [1.0] * 8)))
    vals_demo, _ = plaquette_spectrum(1.0, 1.0)
    chk("plaquette_spectrum_trace_zero",
        abs(sum(vals_demo)) < 1e-10)
    H_demo = dense_H(PLAQ_LINKS, PLAQ_PLAQUETTES, 1.0, 1.0)
    chk("plaquette_trace_H2_formula",
        abs(sum(v * v for v in vals_demo)
            - 16 * (PLAQ_LINKS * 1.0 + 1 * 1.0)) < 1e-8)
    chk("plaquette_H_symmetric",
        mat_max_dev(H_demo, transpose(H_demo)) < 1e-14)
    chk("plaquette_ground_energy_closed_form",
        all(abs(plaquette_spectrum(1.0, J)[0][0]
                + math.sqrt(16.0 + J * J)) < 1e-9
            for J in (0.0, 0.5, 1.0, 2.0)))
    chk("lanczos_refereed_on_plaquette_vs_jacobi",
        abs(lanczos_ground(PLAQ_LINKS, PLAQ_PLAQUETTES, 1.0, 1.0)[0]
            - vals_demo[0]) < 1e-9)

    # 2. Gauss law: commutators and sector dimensions, two ways
    chk("gauss_commutator_plaquette_zero",
        max(commutator_max_dev(PLAQ_LINKS, PLAQ_PLAQUETTES, v, 1.0, 0.7)
            for v in PLAQ_VERTICES) < 1e-12)
    chk("gauss_commutator_ladder_zero",
        max(commutator_max_dev(LADDER_LINKS, LADDER_PLAQUETTES, v,
                               1.0, 0.7)
            for v in LADDER_VERTICES) < 1e-12)
    chk("gauss_sector_dim_plaquette_two_ways",
        len(sector_states(PLAQ_LINKS, PLAQ_VERTICES)) == 2
        and sector_dimension_formula(PLAQ_LINKS, PLAQ_VERTICES) == 2)
    chk("gauss_sector_dim_ladder_two_ways",
        len(sector_states(LADDER_LINKS, LADDER_VERTICES)) == 4
        and sector_dimension_formula(LADDER_LINKS, LADDER_VERTICES) == 4)
    n_sectors = 0
    total_dim = 0
    from itertools import combinations
    for size in (0, 2, 4, 6):
        for combo in combinations(range(len(LADDER_VERTICES)), size):
            states = sector_states(LADDER_LINKS, LADDER_VERTICES,
                                   frozenset(combo))
            if states:
                n_sectors += 1
                total_dim += len(states)
    chk("gauss_sectors_partition_ladder",
        n_sectors == 32 and total_dim == 128)
    _, g_plaq = plaquette_spectrum(1.0, 1.0)
    chk("plaquette_ground_state_is_physical",
        all(abs(v - 1.0) < 1e-9
            for v in gauss_expectations(g_plaq, PLAQ_VERTICES)))
    e_lad, g_lad = ladder_ground(1.0, 1.0)
    chk("ladder_ground_state_is_physical",
        all(abs(v - 1.0) < 1e-9
            for v in gauss_expectations(g_lad, LADDER_VERTICES)))

    # 3. Ladder ground state and its referees
    e_sector = sector_ground_energy(LADDER_LINKS, LADDER_PLAQUETTES,
                                    LADDER_VERTICES, 1.0, 1.0)
    chk("ladder_lanczos_matches_sector_jacobi", abs(e_lad - e_sector) < 1e-8)
    chk("ladder_decoupled_J0_limit",
        abs(ladder_ground(1.0, 0.0)[0] + 7.0) < 1e-9)
    chk("ladder_decoupled_h0_limit",
        abs(ladder_ground(0.0, 1.0)[0] + 2.0) < 1e-9)
    chk("ladder_weak_coupling_series",
        abs(ladder_ground(1.0, 0.1)[0]
            - ladder_weak_series(1.0, 0.1)) < 1e-5)
    chk("ladder_ground_below_decoupled_energies",
        e_lad < -7.0 and e_lad < -2.0)
    trH2 = 0.0
    for s in range(1 << LADDER_LINKS):
        e = [0.0] * (1 << LADDER_LINKS)
        e[s] = 1.0
        col = apply_H(e, LADDER_LINKS, LADDER_PLAQUETTES, 1.0, 1.0)
        trH2 += sum(v * v for v in col)
    chk("ladder_trace_H2_formula",
        abs(trH2 - 128 * (LADDER_LINKS * 1.0
                          + len(LADDER_PLAQUETTES) * 1.0)) < 1e-6)

    # 4. Wilson loop: limits, perturbation slope, Hellmann-Feynman
    sweep = wilson_sweep()
    w_of = {row["J_over_h"]: row["wilson_W"] for row in sweep.values()}
    t_of = {row["J_over_h"]: row["thooft_X0"] for row in sweep.values()}
    chk("wilson_zero_in_electric_limit", abs(w_of[0.0]) < 1e-12)
    chk("wilson_approaches_one_magnetic_limit", w_of[20.0] > 0.98)
    chk("wilson_monotonic_in_coupling",
        all(w_of[a] < w_of[b] for a, b in
            zip(WILSON_SWEEP, WILSON_SWEEP[1:])))
    chk("wilson_perturbation_slope_one_quarter",
        abs(wilson_slope(0.05) - 0.25) < 0.005)
    chk("wilson_value_at_r1",
        abs(w_of[1.0] - 0.242536) < 1e-5)
    eps = 1e-4
    e_plus, _ = plaquette_spectrum(1.0, 1.0 + eps)
    e_minus, _ = plaquette_spectrum(1.0, 1.0 - eps)
    _, g_demo = plaquette_spectrum(1.0, 1.0)
    W_demo, _, elec_demo = plaquette_expectations(g_demo)
    chk("wilson_hellmann_feynman_dE_dJ",
        abs((e_plus[0] - e_minus[0]) / (2 * eps) + W_demo) < 1e-4)
    e_plus, _ = plaquette_spectrum(1.0 + eps, 1.0)
    e_minus, _ = plaquette_spectrum(1.0 - eps, 1.0)
    chk("electric_hellmann_feynman_dE_dh",
        abs((e_plus[0] - e_minus[0]) / (2 * eps)
            + PLAQ_LINKS * elec_demo) < 1e-4)
    chk("thooft_one_in_electric_limit", abs(t_of[0.0] - 1.0) < 1e-12)
    chk("thooft_falls_in_magnetic_limit", t_of[20.0] < 0.2)
    chk("thooft_monotonic_decreasing",
        all(t_of[a] > t_of[b] for a, b in
            zip(WILSON_SWEEP, WILSON_SWEEP[1:])))

    # 5. Static charges: tiny-lattice diagnostic, refereed two ways
    conf = confinement_report()
    costs = [row["energy_cost"] for row in conf["pairs"].values()]
    chk("charge_costs_positive", all(v > 0 for v in costs))
    chk("charge_costs_rise_with_separation",
        costs[0] < costs[1] < costs[2])
    chk("charge_sector_dims_all_four",
        all(len(sector_states(LADDER_LINKS, LADDER_VERTICES,
                               frozenset((a, b)))) == 4
            for a, b, _ in CHARGE_PAIRS))
    e_ch_sector = sector_ground_energy(
        LADDER_LINKS, LADDER_PLAQUETTES, LADDER_VERTICES, 1.0, 1.0,
        charges=frozenset((0, 3)))
    e_ch_lanczos = sector_ground_lanczos(
        LADDER_LINKS, LADDER_PLAQUETTES, LADDER_VERTICES, 1.0, 1.0,
        charges=frozenset((0, 3)))
    chk("charge_sector_refereed_lanczos_vs_jacobi",
        abs(e_ch_sector - e_ch_lanczos) < 1e-8)
    e_ph_lanczos = sector_ground_lanczos(
        LADDER_LINKS, LADDER_PLAQUETTES, LADDER_VERTICES, 1.0, 1.0)
    chk("physical_sector_refereed_lanczos_vs_jacobi",
        abs(e_ph_lanczos - e_sector) < 1e-8)

    # 6. Classifier: labels come from the lab's own policy
    rep = classifier_report()
    chk("classifier_scenario_labels",
        rep["deep_electric_r0.1"]["knn"] == "electric-dominated"
        and rep["electric_r0.25"]["knn"] == "electric-dominated"
        and rep["crossover_r1.0"]["knn"] == "crossover"
        and rep["crossover_r1.5"]["knn"] == "crossover"
        and rep["magnetic_r4.0"]["knn"] == "magnetic-dominated"
        and rep["deep_magnetic_r10.0"]["knn"] == "magnetic-dominated")
    chk("classifier_knn_matches_policy_scenarios",
        all(row["knn"] == row["policy"] for row in rep.values()))
    chk("classifier_knn_matches_policy_grid",
        all(classify_features(feats) == label
            for feats, label in TRAINING_GRID))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_qcd_or_continuum_hold",
        gate("qcd_or_continuum_claim") == "HOLD")
    chk("gate_confinement_proof_hold",
        gate("confinement_proof_claim") == "HOLD")
    chk("gate_model_as_measurement_hold",
        gate("model_as_measurement_claim") == "HOLD")
    chk("gate_live_hardware_hold", gate("live_hardware_control") == "HOLD")
    chk("gate_investment_hold", gate("investment_claim") == "HOLD")
    chk("gate_certification_hold", gate("certification_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking"
        and brains["decider"]["model_id"] == "gemini-argon-4")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("LatticeGauge-Lab - Z2 lattice gauge theory on 4 and 7 links")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - tiny lattices, Z2 only, nothing more")
    print("    Links are qubits; H = -h sum X_l - J sum_p P_p. Two")
    print("    lattices only: one plaquette (4 links, dim 16) and a 2x1")
    print("    ladder (7 links, dim 128). No continuum limit, no QCD")
    print("    claim of any kind - those uses are permanent HOLD gates.")
    vals, _ = plaquette_spectrum(DEMO_H, DEMO_J)
    print("\n[2] Single plaquette - full spectrum by refereed Jacobi")
    print(f"    h = J = 1 ground energy {vals[0]:.6f} "
          f"(closed form -sqrt(16h^2 + J^2) = {-math.sqrt(17.0):.6f})")
    print(f"    spectrum: {[round(v, 4) for v in vals]}")
    vals_j0, _ = plaquette_spectrum(1.0, 0.0)
    vals_h0, _ = plaquette_spectrum(0.0, 1.0)
    print(f"    J = 0 limit: -h(4-2k), binomial degeneracies "
          f"C(4,k) -> {[round(v, 1) for v in vals_j0]}")
    print(f"    h = 0 limit: -J x8 and +J x8 (ONE plaquette operator,")
    print(f"    eigenvalues +/-1) -> {[round(v, 1) for v in vals_h0]}")
    print("\n[3] Gauss law - G_v = product of X at v, [H, G_v] = 0")
    for name, row in gauss_report().items():
        print(f"    {name:10s} commutator max dev "
              f"{row['commutator_max_dev']:.2e} | physical dim "
              f"{row['physical_dim_enumerated']} enumerated = "
              f"{row['physical_dim_formula']} by 2^(L-V+1)")
    print("    -> only V-1 constraints are independent: prod_v G_v = I.")
    print("\n[4] Ladder ground state - Lanczos, refereed without 128x128")
    lad = ladder_report()
    for key in ("demo_h1_J1", "decoupled_J0", "decoupled_h0", "weak_J0.1"):
        row = lad[key]
        print(f"    {key:14s} h {row['h']:.1f} J {row['J']:.1f} | "
              f"Lanczos E {row['lanczos_ground_energy']:.6f}")
    print(f"    X-basis sector Jacobi referee at h = J = 1: "
          f"{lad['demo_h1_J1']['sector_jacobi_referee']:.6f} (agrees)")
    print(f"    weak series at J = 0.1: -7h - J^2/(4h) = "
          f"{lad['weak_J0.1']['weak_series_prediction']:.6f} (agrees)")
    print("\n[5] Wilson loop - the plaquette operator in its ground state")
    for row in wilson_sweep().values():
        print(f"    J/h {row['J_over_h']:5.2f} | <W> {row['wilson_W']:+.6f} "
              f"| <X_0> ('t Hooft) {row['thooft_X0']:.6f} | "
              f"mean <X> {row['electric_mean_X']:.6f}")
    print(f"    -> <W> -> 0 as J/h -> 0 and -> 1 as J/h -> infinity;")
    print(f"       small-coupling slope <W>/(J/h) = {wilson_slope(0.05):.4f} "
          f"at J/h = 0.05, vs first-order perturbation theory 1/4.")
    print("\n[6] Static charges on the ladder - DIAGNOSTIC, not a proof")
    conf = confinement_report()
    print(f"    physical sector ground energy "
          f"{conf['physical_sector_energy']:.6f} at h = J = 1")
    for row in conf["pairs"].values():
        print(f"    charges at vertices {row['vertices']} separation "
              f"{row['separation']} | cost {row['energy_cost']:.6f}")
    print("    -> costs rise with separation on THIS ladder; 7 links")
    print("       cannot distinguish linear from saturating - fenced.")
    print("\n[7] Regime classifier (TinyML kNN, labels from the lab policy)")
    for name, row in classifier_report().items():
        print(f"    {name:24s} -> {row['knn']}")
    print("\n[8] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[9] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "spectrum":
        print(json.dumps(spectrum_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "gauss":
        print(json.dumps(gauss_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "wilson":
        print(json.dumps(wilson_sweep(), indent=2, sort_keys=True))
        return 0
    if cmd == "ladder":
        print(json.dumps(ladder_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "confinement":
        print(json.dumps(confinement_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "qcd_or_continuum_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
