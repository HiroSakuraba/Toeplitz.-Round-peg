"""Exact diagnostics for the v276 rational diamond spirals.

The finite computations check the block formula and the areas of truncations.
The proof of the infinite claims is in the v276 research note.
"""

from fractions import Fraction
import json


Q = Fraction
ORIGIN = (Q(0), Q(0))
DIRECTIONS = ((1, 0), (0, 1), (-1, 0), (0, -1))


def block_sizes(scale):
    """Return the number of edges and the length of a one-sign burst."""
    if scale < 0:
        raise ValueError("scale must be nonnegative")
    return 4 * 4**scale, 4 * 2**scale


def action_weight(scale, start, stop):
    """Unsigned action sum for edges in the half-open index interval."""
    count, _ = block_sizes(scale)
    if not 0 <= start <= stop <= count:
        raise ValueError("edge indices outside the block")
    outer = Q(1, 2**scale)
    step = outer / (2 * count)
    size = stop - start
    sum_j = (start + stop - 1) * size // 2

    def square_prefix(n):
        return n * (n - 1) * (2*n - 1) // 6

    sum_j2 = square_prefix(stop) - square_prefix(start)
    return (size*outer**2 - outer*step*(2*sum_j+size)
            + step**2*(sum_j2+sum_j)) / 2


def block_action(scale, counterwind):
    """Exact endpoint and extremes of the primitive relative to block entry."""
    count, burst = block_sizes(scale)
    if not counterwind:
        total = action_weight(scale, 0, count)
        return total, Q(0), total
    values = [Q(0)]
    for number, start in enumerate(range(0, count, burst)):
        weight = action_weight(scale, start, start+burst)
        values.append(values[-1] + (-1)**number * weight)
    # Between burst endpoints the integral changes monotonically.
    return values[-1], min(values), max(values)


def block_vertices(scale, counterwind):
    """List exact vertices for small diagnostic blocks."""
    count, burst = block_sizes(scale)
    outer = Q(1, 2**scale)
    step = outer / (2 * count)
    for j in range(count + 1):
        sign = -1 if counterwind and (j//burst) % 2 else 1
        x, y = DIRECTIONS[(sign*j) % 4]
        radius = outer - j*step
        yield x*radius, y*radius


def signed_action(points):
    """Integral of (x dy - y dx)/2 along an oriented polygonal path."""
    return sum(((a[0]*b[1]-a[1]*b[0])/2 for a, b in zip(points, points[1:])), Q(0))


def closed_truncation(blocks, counterwind):
    """Include blocks [0, blocks), then join the two inner ends by a diameter."""
    if blocks < 1:
        raise ValueError("at least one block is required")
    arm = list(block_vertices(0, counterwind))
    for scale in range(1, blocks):
        arm.extend(list(block_vertices(scale, counterwind))[1:])
    opposite = [(-x, -y) for x, y in reversed(arm)]
    return arm + [ORIGIN] + opposite + [(Q(0), Q(-1)), arm[0]]


def report():
    cases = {}
    for name, counterwind in (("monotone", False), ("counterwinding", True)):
        values = [block_action(k, counterwind) for k in range(7)]
        cases[name] = {
            "block_endpoint_actions": [str(value[0]) for value in values],
            "block_action_oscillations": [str(high-low) for _, low, high in values],
            "closed_truncation_areas": [
                str(signed_action(closed_truncation(k, counterwind)))
                for k in (1, 2, 3, 4)
            ],
        }
    return cases


if __name__ == "__main__":
    print(json.dumps(report(), indent=2, sort_keys=True))
