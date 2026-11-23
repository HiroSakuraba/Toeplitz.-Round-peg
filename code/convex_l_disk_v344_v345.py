"""Exact rational calculations for the convex redraw and orthogonal L theorems."""

from dataclasses import dataclass
from fractions import Fraction as Q


def convex_root_data(r, R, error):
    """Return (linear coefficient, discriminant) for the large half-diagonal root."""
    r, R, error = map(Q, (r, R, error))
    if not 0 < r <= R or error < 0:
        raise ValueError("require 0 < r <= R and nonnegative Euclidean error")
    linear = r - error
    return linear, linear**2 - 8 * error * (R + error)


def convex_certifies(r, R, error, q):
    """Does the convex-layer test force half-diagonal strictly above q?"""
    linear, _ = convex_root_data(r, R, error)
    error, R, q = map(Q, (error, R, q))
    return q > 0 and q * (linear - q) > 2 * error * (R + error)


@dataclass(frozen=True)
class Certificate:
    method: str
    side_squared: Q
    # The two anchor coordinates are (u,0) and (0,span-u) on the elbow family.
    span: Q | None = None
    u_interval: tuple[Q, Q] | None = None

    def elbow_square(self, u):
        u = Q(u)
        if self.u_interval is None or not self.u_interval[0] <= u <= self.u_interval[1]:
            raise ValueError("u must lie in the marked elbow interval")
        t = self.span
        v = t - u
        return ((u, 0), (t, u), (v, t), (0, v))


@dataclass(frozen=True)
class OrthogonalL:
    width: Q
    height: Q
    vertical_arm_width: Q
    horizontal_arm_height: Q

    def __post_init__(self):
        for name in ("width", "height", "vertical_arm_width", "horizontal_arm_height"):
            object.__setattr__(self, name, Q(getattr(self, name)))
        if not (0 < self.vertical_arm_width < self.width and
                0 < self.horizontal_arm_height < self.height):
            raise ValueError("each arm must have positive thickness below its outer span")

    def polygon(self):
        W, H = self.width, self.height
        a, b = self.vertical_arm_width, self.horizontal_arm_height
        return ((0, 0), (W, 0), (W, b), (a, b), (a, H), (0, H))

    def fits_disk(self, radius):
        """Exact feasibility of a closed disk with a given rational radius."""
        r = Q(radius)
        if r < 0:
            raise ValueError("negative radius")
        a, b = self.vertical_arm_width, self.horizontal_arm_height
        return (2 * r <= min(self.width, self.height) and
                max(a - r, 0)**2 + max(b - r, 0)**2 >= r**2)

    def reaches_disk_target(self, side_squared):
        """Exact comparison of a squared side with twice the squared inradius."""
        target = Q(side_squared) / 2
        if target < 0:
            return False
        if min(self.width, self.height)**2 <= 4 * target:
            return True
        A, B = sorted((self.vertical_arm_width, self.horizontal_arm_height))
        if B >= 2 * A:
            return B**2 <= 4 * target
        # Compare (A+B-sqrt(2AB))^2 <= target with one sign-safe squaring.
        c, d = A + B, 2 * A * B
        excess = c**2 + d - target
        return excess <= 0 or excess**2 <= 4 * c**2 * d

    def certificate(self):
        W, H = self.width, self.height
        a, b = self.vertical_arm_width, self.horizontal_arm_height
        m = min(W, H)
        if max(a, b) > m:
            return Certificate("tall walls", m**2)
        if max(a, b) == m:
            plateau = W - a if b == m else H - b
            d = min(plateau, m / 2)
            return Certificate("end plateau", d**2 + (m - d)**2)

        t = min(W, H, a + b)
        lo, hi = t - a, b
        if t < W:
            lo = hi = b
        if t < H:
            lo = hi = t - a
        u = max(lo, min(t / 2, hi))
        return Certificate("elbow", u**2 + (t - u)**2, t, (lo, hi))
