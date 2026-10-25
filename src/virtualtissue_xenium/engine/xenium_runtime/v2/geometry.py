"""Bounded power (Laguerre) diagram used as the v2 cell geometry.

Every body i is a site p_i with radius rho_i and power weight w_i = rho_i**2.
Its simulated cell is its power cell intersected with its own disk. On the
power bisector of i and j, membership of disk i and disk j coincides, so the
clipped cells tile the union of the disks: two simulated cells can never
overlap, cells at free edges stay round, and gaps larger than the disks stay
empty. The diagram is computed once per mechanics substep from the lower
convex hull of the lifted sites (scipy/Qhull); areas, shared interface
lengths and contacts are evaluated in closed form, vectorised over edges.
"""
from dataclasses import dataclass

import numpy as np
from scipy.spatial import ConvexHull


class GeometryError(RuntimeError):
    pass


def _cross(a, b):
    return a[..., 0] * b[..., 1] - a[..., 1] * b[..., 0]


def _angle(a, b):
    return np.arctan2(_cross(a, b), (a * b).sum(-1))


def _clip_params(A, B, r):
    """Segment A->B relative to a disk centre: parameters s1 <= s2 of its inside part."""
    D = B - A
    a = (D * D).sum(-1)
    b = 2 * (A * D).sum(-1)
    c = (A * A).sum(-1) - r * r
    disc = b * b - 4 * a * c
    ok = (disc > 0) & (a > 0)
    sq = np.sqrt(np.where(ok, disc, 0.))
    safe = np.where(a > 0, a, 1.)
    s1 = np.where(ok, np.clip((-b - sq) / (2 * safe), 0, 1), 0.)
    s2 = np.where(ok, np.clip((-b + sq) / (2 * safe), 0, 1), 0.)
    return s1, s2, D


def sector_triangle_area(A, B, r):
    """Signed area of the triangle (0, A, B) intersected with the disk |x| <= r."""
    s1, s2, D = _clip_params(A, B, r)
    P1 = A + s1[..., None] * D
    P2 = A + s2[..., None] * D
    return .5 * r * r * _angle(A, P1) + .5 * _cross(P1, P2) + .5 * r * r * _angle(P2, B)


@dataclass
class PowerGeometry:
    """Geometry of n bodies. Pair arrays hold only real-real power edges."""
    xy: np.ndarray
    rho: np.ndarray
    area: np.ndarray
    hidden: np.ndarray
    edge_i: np.ndarray       # every power edge with a real first body (second may be a dummy)
    edge_j: np.ndarray       # index >= n means a bounding dummy site
    edge_a: np.ndarray       # oriented counter-clockwise for body edge_i
    edge_b: np.ndarray
    pair_i: np.ndarray       # every real-real power-adjacent pair
    pair_j: np.ndarray
    pair_length: np.ndarray  # shared interface inside both disks (um); 0 if the disks do not meet
    pair_gap: np.ndarray     # 0 if they share an interface, else d_ij - rho_i - rho_j (um)
    pair_normal: np.ndarray  # unit vector from site i to site j
    pair_distance: np.ndarray

    @property
    def n(self):
        return len(self.xy)

    def contacts(self, max_gap_um=1.):
        """Simulated contacts with the measured definition (H04): cells within ``max_gap_um``.

        Pairs sharing an interface have gap 0. For power-adjacent pairs whose
        disks do not meet, d - rho_i - rho_j is the disk gap (a lower bound
        on the gap between the clipped cells).
        """
        keep = self.pair_gap <= max_gap_um
        return self.pair_i[keep], self.pair_j[keep], self.pair_length[keep]

    def neighbour_lists(self, max_gap_um=1.):
        i, j, _ = self.contacts(max_gap_um)
        out = [[] for _ in range(self.n)]
        for a, b in zip(i.tolist(), j.tolist()):
            out[a].append(b)
            out[b].append(a)
        return out

    def polygon(self, k, quad_segs=24):
        """Shapely polygon of body k (for fidelity checks and export only)."""
        from shapely.geometry import Point, Polygon
        disk = Point(*self.xy[k]).buffer(float(self.rho[k]), quad_segs=quad_segs)
        if self.hidden[k]:
            return Polygon()
        first = self.edge_i == k
        second = self.edge_j == k
        a = np.vstack([self.edge_a[first], self.edge_b[second]])
        b = np.vstack([self.edge_b[first], self.edge_a[second]])
        mid = (a + b) / 2
        centre = np.vstack([a, b]).mean(0)
        order = np.argsort(np.arctan2(mid[:, 1] - centre[1], mid[:, 0] - centre[0]))
        ring = a[order]
        if len(ring) < 3:
            return Polygon()
        cell = Polygon(ring)
        if not cell.is_valid:
            cell = cell.buffer(0)
        return cell.intersection(disk)


def power_geometry(xy, rho, n_dummy=64):
    """Clipped power diagram of sites ``xy`` (n x 2, um) with radii ``rho`` (um)."""
    xy = np.asarray(xy, dtype=float)
    rho = np.asarray(rho, dtype=float)
    n = len(xy)
    if n < 1 or xy.shape != (n, 2) or rho.shape != (n,):
        raise GeometryError("Need n x 2 sites and n radii")
    if not (np.isfinite(xy).all() and np.isfinite(rho).all() and (rho > 0).all()):
        raise GeometryError("Nonfinite site or nonpositive radius")
    centre = xy.mean(0)
    P = xy - centre
    extent = float(np.sqrt((P * P).sum(1).max()))
    R = extent + 4 * float(rho.max()) + 50.
    theta = (np.arange(n_dummy) + .5) * 2 * np.pi / n_dummy
    X = np.vstack([P, R * np.c_[np.cos(theta), np.sin(theta)]])
    W = np.r_[rho * rho, np.zeros(n_dummy)]
    H = (X * X).sum(1) - W
    hull = ConvexHull(np.c_[X, H])
    tri = hull.simplices[hull.equations[:, 2] < 0]
    pa, pb, pc = X[tri[:, 0]], X[tri[:, 1]], X[tri[:, 2]]
    m0, m1 = 2 * (pb - pa), 2 * (pc - pa)
    r0, r1 = H[tri[:, 1]] - H[tri[:, 0]], H[tri[:, 2]] - H[tri[:, 0]]
    det = m0[:, 0] * m1[:, 1] - m0[:, 1] * m1[:, 0]
    with np.errstate(divide="ignore", invalid="ignore"):
        C = np.c_[(r0 * m1[:, 1] - m0[:, 1] * r1) / det, (m0[:, 0] * r1 - r0 * m1[:, 0]) / det]
    e = np.concatenate([tri[:, [0, 1]], tri[:, [1, 2]], tri[:, [2, 0]]])
    t = np.tile(np.arange(len(tri)), 3)
    e.sort(axis=1)
    order = np.lexsort((e[:, 1], e[:, 0]))
    e, t = e[order], t[order]
    same = (e[1:] == e[:-1]).all(1)
    dup = np.flatnonzero(same)
    single = np.ones(len(e), bool)
    single[dup] = False
    single[dup + 1] = False
    if (e[single, 0] < n).any():
        raise GeometryError("A real site lies on the unbounded hull; enlarge the dummy ring")
    I, J, T1, T2 = e[dup, 0], e[dup, 1], t[dup], t[dup + 1]
    real = I < n
    I, J, T1, T2 = I[real], J[real], T1[real], T2[real]
    a, b = C[T1], C[T2]
    if not (np.isfinite(a).all() and np.isfinite(b).all()):
        raise GeometryError("Degenerate power vertex")
    normal = X[J] - X[I]
    d = np.c_[-normal[:, 1], normal[:, 0]]
    flip = ((b - a) * d).sum(1) < 0
    a, b = np.where(flip[:, None], b, a), np.where(flip[:, None], a, b)
    area = np.bincount(I, weights=sector_triangle_area(a - X[I], b - X[I], rho[I]), minlength=n)
    realj = J < n
    Jr = J[realj]
    area += np.bincount(Jr, weights=sector_triangle_area(b[realj] - X[Jr], a[realj] - X[Jr], rho[Jr]), minlength=n)
    hidden = np.ones(n, bool)
    used = tri[tri < n]
    hidden[used] = False
    area = np.where(hidden, 0., np.maximum(area, 0.))
    s1, s2, D = _clip_params(a[realj] - X[I[realj]], b[realj] - X[I[realj]], rho[I[realj]])
    length = np.sqrt((D * D).sum(1)) * (s2 - s1)
    pi, pj = I[realj], Jr
    vec = X[pj] - X[pi]
    dist = np.sqrt((vec * vec).sum(1))
    length = np.where(length > 1e-9, length, 0.)
    gap = np.where(length > 0, 0., np.maximum(dist - rho[pi] - rho[pj], 0.))
    return PowerGeometry(xy=xy, rho=rho, area=area, hidden=hidden,
                         edge_i=I, edge_j=J, edge_a=a + centre, edge_b=b + centre,
                         pair_i=pi, pair_j=pj, pair_length=length, pair_gap=gap,
                         pair_normal=vec / dist[:, None], pair_distance=dist)
