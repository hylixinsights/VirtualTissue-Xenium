"""Conservative finite-volume reaction diffusion on a reflecting rectangular grid.

All concentrations and source strengths are SIM dimensionless quantities. The
rectangle encloses the mechanical halo; its boundary has zero normal flux.
Each quasi-steady solve satisfies D L c + k c = s. Deposition is conservative
bilinear cloud-in-cell; sinks are first order, so solutions stay nonnegative.
"""
import numpy as np
from scipy.sparse import diags, eye, kron
from scipy.sparse.linalg import factorized, spsolve


class FieldGrid:
    def __init__(self, bounds, spacing, specs):
        self.origin = np.asarray(bounds[:2], float)
        self.spacing = float(spacing)
        self.nx, self.ny = np.maximum(3, np.ceil((np.asarray(bounds[2:]) - self.origin) / spacing).astype(int) + 1)
        self.shape = (int(self.ny), int(self.nx))
        self.size = int(self.nx * self.ny)
        self.specs = specs
        def lap(n):
            return diags((-np.ones(n-1), np.r_[1., np.full(n-2, 2.), 1.], -np.ones(n-1)), (-1, 0, 1), format="csc")
        self.L = (kron(eye(self.ny), lap(self.nx)) + kron(lap(self.ny), eye(self.nx))).tocsc() / spacing**2
        self.values = {name: np.zeros(self.shape) for name in specs}
        self.solvers = {}
        self.audit = {}

    def weights(self, xy):
        u = np.clip((np.asarray(xy) - self.origin) / self.spacing, [0, 0], [self.nx-1.00000001, self.ny-1.00000001])
        ij = np.floor(u).astype(int)
        f = u - ij
        base = ij[:, 1] * self.nx + ij[:, 0]
        return np.c_[base, base+1, base+self.nx, base+self.nx+1], np.c_[(1-f[:,0])*(1-f[:,1]), f[:,0]*(1-f[:,1]), (1-f[:,0])*f[:,1], f[:,0]*f[:,1]]

    def deposit(self, xy, amounts):
        if not len(xy):
            return np.zeros(self.shape)
        idx, weights = self.weights(xy)
        return np.bincount(idx.ravel(), (weights * np.asarray(amounts)[:,None]).ravel(), self.size).reshape(self.shape) / self.spacing**2

    def sample(self, name, xy):
        idx, weights = self.weights(xy)
        return (self.values[name].ravel()[idx] * weights).sum(axis=1)

    def gradient(self, name, xy):
        gy, gx = np.gradient(self.values[name], self.spacing)
        idx, weights = self.weights(xy)
        return np.c_[(gx.ravel()[idx]*weights).sum(1), (gy.ravel()[idx]*weights).sum(1)]

    def solve(self, name, source, sink=None):
        p = self.specs[name]
        source = np.asarray(source, float).reshape(self.shape)
        sink = np.zeros(self.shape) if sink is None else np.asarray(sink, float).reshape(self.shape)
        if (source < 0).any() or (sink < 0).any() or not np.isfinite(source + sink).all():
            raise ValueError("Sources and sinks must be finite and nonnegative")
        loss = p["decay"] + sink.ravel()
        A = p["diffusion"] * self.L + diags(loss, format="csc")
        if not sink.any():
            if name not in self.solvers:
                self.solvers[name] = factorized(A)
            c = self.solvers[name](source.ravel())
        else:
            c = spsolve(A, source.ravel())
        if not np.isfinite(c).all() or c.min() < -1e-9:
            raise ArithmeticError("Invalid field solution")
        residual = A @ c - source.ravel()
        balance = float((loss*c).sum() - source.sum()) * self.spacing**2
        self.audit[name] = dict(source_integral=float(source.sum()*self.spacing**2),
                                loss_integral=float((loss*c).sum()*self.spacing**2),
                                balance_error=balance, residual_max=float(abs(residual).max()))
        self.values[name] = np.maximum(c, 0).reshape(self.shape)
        return self.values[name]

    def frame(self):
        return {k: np.round(v[::2,::2], 5).ravel().tolist() for k,v in self.values.items()}

    def metadata(self):
        return dict(origin_um=self.origin.tolist(), spacing_um=self.spacing*2, shape=list(self.values[next(iter(self.values))][::2,::2].shape),
                    solver_spacing_um=self.spacing, display_sampling="Every second solver node, rounded to five decimals; full kernel state is separate",
                    boundary="zero flux on the rectangle enclosing the halo", units="SIM dimensionless")
