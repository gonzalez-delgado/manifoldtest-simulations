"""
HHG independence test -- Heller, Heller & Gorfine (2013),
"A consistent multivariate test of association based on ranks of distances",
Biometrika 100(2):503-510.  arXiv:1201.3522.

Called by the R function run_hhg() of the package; it can also be run
directly:

    python3 hhg.py hyperbolic INPUT OUTPUT [--n-perm 199]
    python3 hhg.py grassmann  INPUT OUTPUT [--n-perm 199]
    python3 hhg.py grassmann  INPUT OUTPUT --distances N_ROWS
    python3 hhg.py --selftest

  hyperbolic : INPUT written by hyperbolic_shared_data() (tangent vectors
               u1..u3, v1..v3 at p = (1,0,0,0) of H^3);
               OUTPUT columns exp, scan, xval, rep, n, p_value
  grassmann  : INPUT written by figure4_grassmann.R (Stiefel representatives
               x1..x8, y1..y8 of points of Gr(4,2));
               OUTPUT columns scan, xval, rep, n, HHG
  --distances: write instead the geodesic distance matrix of the first
               N_ROWS rows of a Grassmann INPUT (cross-check with R)

Definition (paper, Section 2)
-----------------------------
For every ordered pair (i, j) with i != j, take x_i and y_i as centres and
the distances d(x_i, x_j) and d(y_i, y_j) as radii.  Cross-classify the
remaining N-2 observations k by

    I{ d(x_i, x_k) <= d(x_i, x_j) }   and   I{ d(y_i, y_k) <= d(y_i, y_j) },

giving the 2x2 table A = (A11, A12; A21, A22) with margins A1., A2., A.1, A.2.
The contribution of the pair is Pearson's chi-square statistic

    S(i, j) = (N-2) (A12 A21 - A11 A22)^2 / (A1. A2. A.1 A.2),

set to 0 whenever a margin vanishes.  The test statistic is

    T = sum_i sum_{j != i} S(i, j),

and the p-value is obtained by permuting the indices of the Y sample.

This module uses the standard Monte-Carlo permutation p-value
(1 + #{T_perm >= T_obs}) / (1 + n_perm), which is the convention used by
energy::dcov.test and Ball::bcov.test, so the HHG curve is directly
comparable with the other curves in the figures.
"""
import argparse
import hashlib
import multiprocessing as mp
import sys
import time

import numpy as np
import pandas as pd

N_DIM = 4
K_DIM = 2
N_TAN = N_DIM * K_DIM


# ----------------------------------------------------------------
# distance helpers
# ----------------------------------------------------------------
def dist_hyperbolic_from_tangent(U):
    """Geodesic distance matrix on H^3 from tangent vectors u_i in T_p(H^3).

    With x_i = (cosh r_i, sinh r_i * d_i) and cos(theta) = <d_i, d_j>,
        cosh d_H(x_i, x_j) = cosh r_i cosh r_j - sinh r_i sinh r_j cos(theta).
    Matches hyperbolic_dist() of the R package, including the
    max(., 1 + 1e-12) clamp on the Lorentzian inner product.
    """
    r = np.sqrt((U ** 2).sum(1))
    ch, sh = np.cosh(r), np.sinh(r)
    # unit directions (guard r = 0)
    D = U / np.where(r > 1e-12, r, 1.0)[:, None]
    cos_theta = np.clip(D @ D.T, -1.0, 1.0)
    C = ch[:, None] * ch[None, :] - sh[:, None] * sh[None, :] * cos_theta
    return np.arccosh(np.maximum(C, 1.0 + 1e-12))


def stiefel_from_rows(flat):
    """Stiefel matrices U of a sample.

    flat : (N, n*k) array.  Each row holds as.vector(U), as written by
    figure4_grassmann.R, namely the column-major flattening of an n x k
    matrix, so each row is reshaped to (k, n) and transposed.
    """
    return np.asarray(flat).reshape(-1, K_DIM, N_DIM).transpose(0, 2, 1)


def gr_dist_from_rows(flat):
    """Geodesic distance matrix on Gr(4,2): the square root of the sum of the
    squared principal angles between the subspaces."""
    U = stiefel_from_rows(flat)
    M = np.einsum("nar,mas->nmrs", U, U)
    s = np.linalg.svd(M)[1]
    th = np.arccos(np.clip(s, 0.0, 1.0))
    D = np.sqrt((th ** 2).sum(-1))
    np.fill_diagonal(D, 0.0)
    return D


# ----------------------------------------------------------------
# statistic
# ----------------------------------------------------------------
def hhg_stat(DX, DYp):
    """T = sum_i sum_{j != i} S(i, j).

    DX  : (N, N) distance matrix of the X sample
    DYp : (N, N) distance matrix of the Y sample *as paired* (already permuted)

    Point i is removed from the counts by pushing its distance to +inf, and
    the diagonal k = j is removed by subtracting one from nU, nV and A11.
    """
    N = DX.shape[0]
    m = N - 2
    T = 0.0
    for i in range(N):
        a = DX[i].copy(); a[i] = np.inf
        b = DYp[i].copy(); b[i] = np.inf
        U = a[:, None] <= a[None, :]        # U[k, j] = 1{a_k <= a_j}
        V = b[:, None] <= b[None, :]
        nU = U.sum(0) - 1.0                 # #{k != i, j : a_k <= a_j}
        nV = V.sum(0) - 1.0
        A11 = (U & V).sum(0) - 1.0
        A12 = nU - A11
        A21 = nV - A11
        A22 = m - nU - nV + A11
        den = nU * (m - nU) * nV * (m - nV)
        num = m * (A12 * A21 - A11 * A22) ** 2
        S = np.where(den > 0, num / np.where(den > 0, den, 1.0), 0.0)
        S[i] = 0.0                          # j = i is not part of the sum
        T += S.sum()
    return T


def hhg_stat_bruteforce(DX, DYp):
    """Straight O(N^3) transcription of the definition; used to validate."""
    N = DX.shape[0]
    m = N - 2
    T = 0.0
    for i in range(N):
        for j in range(N):
            if i == j:
                continue
            A = np.zeros((2, 2))
            for k in range(N):
                if k == i or k == j:
                    continue
                p = 0 if DX[i, k] <= DX[i, j] else 1
                q = 0 if DYp[i, k] <= DYp[i, j] else 1
                A[p, q] += 1.0
            A11, A12, A21, A22 = A[0, 0], A[0, 1], A[1, 0], A[1, 1]
            r1, r2 = A11 + A12, A21 + A22
            c1, c2 = A11 + A21, A12 + A22
            if r1 * r2 * c1 * c2 > 0:
                T += m * (A12 * A21 - A11 * A22) ** 2 / (r1 * r2 * c1 * c2)
    return T


def hhg_pvalue(DX, DY, n_perm, rng):
    """Monte-Carlo permutation p-value, permuting the Y indices."""
    N = DX.shape[0]
    T_obs = hhg_stat(DX, DY)
    cnt = 0
    for _ in range(n_perm):
        perm = rng.permutation(N)
        if hhg_stat(DX, DY[np.ix_(perm, perm)]) >= T_obs:
            cnt += 1
    return (1.0 + cnt) / (1.0 + n_perm), T_obs


# ----------------------------------------------------------------
# the two designs
# ----------------------------------------------------------------
# key columns of a replication, distance builder, input columns, output
# columns, and pandas options; the seed of a replication is derived from its
# key, so it does not depend on the other replications of the file
SPACES = {
    "hyperbolic": dict(
        key=["exp", "scan", "xval", "rep"],
        dist=dist_hyperbolic_from_tangent,
        xcols=["u1", "u2", "u3"], ycols=["v1", "v2", "v3"],
        out=["exp", "scan", "xval", "rep", "n", "p_value"],
        read=dict(), report=200),
    "grassmann": dict(
        key=["scan", "xval", "rep"],
        dist=gr_dist_from_rows,
        xcols=[f"x{i}" for i in range(1, N_TAN + 1)],
        ycols=[f"y{i}" for i in range(1, N_TAN + 1)],
        out=["scan", "xval", "rep", "n", "HHG"],
        # keep_default_na = False: the null scenario is spelled "null" in the
        # scan column, which read_csv would otherwise turn into a missing value
        read=dict(keep_default_na=False), report=100),
}

_G = {}


def seed_for(key):
    h = hashlib.md5("|".join(str(k) for k in key).encode()).digest()
    return int.from_bytes(h[:4], "little")


def init_worker(space, path, n_perm):
    sp = SPACES[space]
    df = pd.read_csv(path, **sp["read"])
    idx = df.groupby(sp["key"], sort=False).indices
    _G.update(sp=sp, df=df, keys=list(idx.keys()), idx=idx, n_perm=n_perm)


def run_one(k):
    sp = _G["sp"]
    key = _G["keys"][k]
    sub = _G["df"].iloc[_G["idx"][key]]
    x = sub[sp["xcols"]].to_numpy()
    y = sub[sp["ycols"]].to_numpy()
    rng = np.random.default_rng(seed_for(key))
    p, _ = hhg_pvalue(sp["dist"](x), sp["dist"](y), _G["n_perm"], rng)
    *lab, xval, rep = key
    return (*lab, float(xval), int(rep), int(len(x)), p)


def run(space, path_in, path_out, n_perm):
    sp = SPACES[space]
    n_tasks = (pd.read_csv(path_in, usecols=sp["key"], **sp["read"])
               .groupby(sp["key"], sort=False).ngroups)
    n_proc = max(1, mp.cpu_count() - 2)
    print(f"HHG on {space}: n_perm={n_perm}, workers={n_proc}, tasks={n_tasks}",
          flush=True)
    t0 = time.time()
    rows = []
    with mp.Pool(n_proc, initializer=init_worker,
                 initargs=(space, path_in, n_perm)) as pool:
        for i, r in enumerate(
                pool.imap_unordered(run_one, range(n_tasks), chunksize=1), 1):
            rows.append(r)
            if i % sp["report"] == 0:
                el = time.time() - t0
                print(f"  {i}/{n_tasks}  ({el/60:.1f} min elapsed, "
                      f"eta {el/i*(n_tasks-i)/60:.1f} min)", flush=True)
    out = pd.DataFrame(rows, columns=sp["out"])
    out.to_csv(path_out, index=False)
    print(f"wrote {path_out}  ({len(out)} rows, {(time.time()-t0)/60:.1f} min)",
          flush=True)


def selftest():
    rng = np.random.default_rng(0)
    N = 9
    A = rng.normal(size=(N, N)); A = A + A.T
    B = rng.normal(size=(N, N)); B = B + B.T
    np.fill_diagonal(A, 0); np.fill_diagonal(B, 0)
    print(f"hhg_stat vs brute force: {abs(hhg_stat(A, B) - hhg_stat_bruteforce(A, B)):.3e}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("space", nargs="?", choices=sorted(SPACES))
    ap.add_argument("input", nargs="?")
    ap.add_argument("output", nargs="?")
    ap.add_argument("--n-perm", type=int, default=199)
    ap.add_argument("--distances", type=int, metavar="N_ROWS")
    ap.add_argument("--selftest", action="store_true")
    a = ap.parse_args()
    if a.selftest:
        selftest()
        sys.exit()
    if a.space is None or a.input is None or a.output is None:
        ap.error("space, input and output are required")
    if a.distances is not None:
        if a.space != "grassmann":
            ap.error("--distances is only available for grassmann")
        df = pd.read_csv(a.input, nrows=a.distances, keep_default_na=False)
        D = gr_dist_from_rows(df[SPACES["grassmann"]["xcols"]].to_numpy())
        np.savetxt(a.output, D, delimiter=",", fmt="%.12f")
        print(f"wrote {a.output} ({a.distances} rows)")
    else:
        run(a.space, a.input, a.output, a.n_perm)
