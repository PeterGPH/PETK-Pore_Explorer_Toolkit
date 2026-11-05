"""GPU backend (sem/gpu_solver.py): checks that run on the CPU with torch."""
import itertools
from types import SimpleNamespace

import numpy as np
import pytest

torch = pytest.importorskip("torch")
pytest.importorskip("dolfinx")

import scipy.sparse as sp  # noqa: E402
import scipy.sparse.linalg as spla  # noqa: E402

from sem import gpu_solver as gs  # noqa: E402
from sem.conductivity_models import SimpleConductivityModel  # noqa: E402
from sem.vertical_movement_sem import AnalyteOverlapError, VerticalMovementSEM  # noqa: E402

torch.set_default_dtype(torch.float64)


def _kuhn_path(origin, perm):
    path = [np.array(origin)]
    for a in perm:
        nxt = path[-1].copy()
        nxt[a] += 1
        path.append(nxt)
    return path


def _generic_p1(sig6, n, h):
    """P1 stiffness matrix of the Kuhn mesh, assembled element by element from coordinates."""
    nv = tuple(k + 1 for k in n)
    vid = np.arange(np.prod(nv)).reshape(nv)
    rows, cols, vals = [], [], []
    for i, j, k in itertools.product(*(range(m) for m in n)):
        for t, perm in enumerate(gs.KUHN):
            path = _kuhn_path((i, j, k), perm)
            T = np.c_[np.ones(4), np.array(path) * np.array(h)]
            G = np.linalg.inv(T)[1:].T
            Ke = sig6[0, t, i, j, k] * abs(np.linalg.det(T)) / 6 * G @ G.T
            ids = [vid[tuple(p)] for p in path]
            for a, b in itertools.product(range(4), repeat=2):
                rows.append(ids[a]); cols.append(ids[b]); vals.append(Ke[a, b])
    return sp.csr_matrix((vals, (rows, cols)), shape=(vid.size,) * 2).toarray(), vid


def test_kuhn_stencil_equals_generic_p1_assembly():
    n, h = (4, 3, 5), (1.0, 1.5, 0.8)
    sig6 = np.random.default_rng(0).uniform(0.1, 2.0, (1, 6, *n))
    K, vid = _generic_p1(sig6, n, h)
    L = gs.fem_level(torch.tensor(sig6), h)
    inner = vid[:, :, 1:-1].ravel()
    nu = (n[0] + 1, n[1] + 1, n[2] - 1)
    E = np.eye(int(np.prod(nu)))
    A = np.stack([gs.apply(L, torch.tensor(E[c]).reshape(1, *nu))[0].numpy().ravel()
                  for c in range(E.shape[0])], 1)
    np.testing.assert_allclose(A, K[np.ix_(inner, inner)], atol=1e-12)
    top, last = vid[:, :, -1].ravel(), vid[:, :, -2].ravel()
    np.testing.assert_allclose(-np.diag(K[np.ix_(last, top)]), L["gt"][0].numpy().ravel(), atol=1e-12)


@pytest.mark.parametrize("n", [(12, 12, 12), (10, 8, 7)])
def test_uniform_box_current_is_exact(n):
    h = tuple(60.0 / k for k in n)
    s6 = torch.full((1, 6, *n), 1.12)
    levels, factors = gs.setup_hierarchy(gs.fem_level(s6, h))
    b = gs.rhs(levels[0], 0.1)
    x, _ = gs.pcg(levels, factors, b, torch.zeros_like(b), 25)
    exact = 1.12 * 60.0 * 60.0 * 0.1 / 60.0
    assert float(gs.flux_current(levels[0], x, 0.1)[0]) == pytest.approx(exact, rel=1e-8)
    assert float(gs.power_current(levels[0], x, 0.1)[0]) == pytest.approx(exact, rel=1e-10)


def test_multigrid_cg_matches_direct_solve_with_membrane_and_insulator():
    n, h = (20, 20, 20), (3.0, 3.0, 3.0)
    grid = gs.BoxGrid(n, h, (-30.0, -30.0, -30.0))
    s6 = np.empty((1, 6, *n))
    for t, perm in enumerate(gs.KUHN):
        P = grid.points(gs.kuhn_centroid(perm))
        R, r = np.hypot(P[:, 0], P[:, 1]), np.linalg.norm(P, axis=1)
        s = np.where((np.abs(P[:, 2]) < 12) & (R > 12), 1e-8, 1.12)
        s6[0, t] = np.minimum(s, np.where(r < 9, 1.12e-7, 1.12)).reshape(n)
    levels, factors = gs.setup_hierarchy(gs.fem_level(torch.tensor(s6), h))
    L, b = levels[0], gs.rhs(levels[0], 0.1)
    x, _ = gs.pcg(levels, factors, b, torch.zeros_like(b), 30)
    nu = L["D"].shape[1:]
    idx = np.arange(int(np.prod(nu))).reshape(nu)
    rows, cols, vals = [], [], []
    for g, a0, a1 in ((L["gx"], idx[:-1], idx[1:]), (L["gy"], idx[:, :-1], idx[:, 1:]),
                      (L["gz"], idx[..., :-1], idx[..., 1:])):
        gi = g[0].numpy().ravel()
        rows += [a0.ravel(), a1.ravel()]; cols += [a1.ravel(), a0.ravel()]; vals += [-gi, -gi]
    A = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(idx.size,) * 2) + sp.diags(L["D"][0].numpy().ravel())
    xd = torch.tensor(spla.spsolve(A.tocsc(), b[0].numpy().ravel())).reshape(1, *nu)
    assert float(gs.flux_current(L, x, 0.1)[0]) == pytest.approx(float(gs.flux_current(L, xd, 0.1)[0]), rel=1e-8)


def test_grid_potential_is_p1_interpolation():
    rng = np.random.default_rng(1)
    grid = gs.BoxGrid((3, 4, 2), (1.0, 2.0, 1.5), (-1.0, -2.0, 0.0))
    phi = rng.normal(size=(4, 5, 3))
    pot = gs.GridPotential(grid, phi)
    x, y, z = (np.sort(rng.uniform(lo, lo + k * hh, 5)) for lo, k, hh in zip(grid.lo, grid.n, grid.h))
    got = pot.sample_on_grid(x, y, z)
    for (i, xi), (j, yj), (k, zk) in itertools.product(enumerate(x), enumerate(y), enumerate(z)):
        p = (np.array([xi, yj, zk]) - grid.lo) / grid.h
        c = np.floor(p).astype(int)
        for perm in gs.KUHN:   # find the tetrahedron that contains p and use barycentric weights
            path = _kuhn_path(c, perm)
            T = np.c_[np.ones(4), np.array(path, float)]
            lam = np.linalg.solve(T.T, np.r_[1.0, p])
            if (lam >= -1e-12).all():
                want = sum(l * phi[tuple(v)] for l, v in zip(lam, path))
                break
        assert got[i, j, k] == pytest.approx(want, abs=1e-12)
    np.testing.assert_allclose(pot.sample_on_grid([0.0], [0.0], [0.0])[0, 0, 0], phi[1, 1, 0])
    assert pot.sample_on_grid([99.0], [0.0], [0.0], fill_value=-7.0)[0, 0, 0] == -7.0


def _fake_sem(base, membrane, prevent):
    model = SimpleConductivityModel(bulk_conductivity=1.12)
    return SimpleNamespace(cutoff=5.0, conductivity_model=model, prevent_analyte_overlap=prevent,
                           membrane_conductivity=membrane)


@pytest.mark.parametrize("offset", [(0.5, 0.5, 0.5), (0.75, 0.25, 0.5)])
def test_analyte_sigma_and_overlap_match_petk(offset):
    rng = np.random.default_rng(2)
    grid = gs.BoxGrid((30, 30, 30), (2.0, 2.0, 2.0), (-30.0, -30.0, -30.0))
    pts = grid.points(offset)
    membrane = 1e-8
    base = np.where(np.abs(pts[:, 2]) > 20, membrane, 1.12)
    radii = rng.choice([1.52, 1.55, 1.70, 1.80], 150)
    model = SimpleConductivityModel(bulk_conductivity=1.12)
    sig = gs.AnalyteSigma(grid, base, radii, 5.0, model, "cpu", offset=offset, membrane_conductivity=membrane)
    clear = rng.normal(0, 4, (150, 3))                 # well inside the conducting slab
    touching = clear + np.array([0.0, 0.0, 17.5])      # hard cores reach the membrane
    got, overlap = sig.batch([clear, touching], with_overlap=True)
    for k, pos in enumerate((clear, touching)):
        fake = _fake_sem(base, membrane, prevent=True)
        try:
            want = VerticalMovementSEM.calculate_analyte_conductivity_modification(fake, pts, pos, radii, base)
            petk_overlap = False
            np.testing.assert_allclose(got[k].numpy(), want, rtol=0, atol=1e-14)
        except AnalyteOverlapError:
            petk_overlap = True
        assert bool(overlap[k]) == petk_overlap
    assert not bool(overlap[0]) and bool(overlap[1])
