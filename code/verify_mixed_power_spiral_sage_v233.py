"""Independent interval bounds for the three v233 source boxes."""

from sage.all import QQ, RealIntervalField, matrix, vector

from mixed_power_spiral_v233 import FLOORS, GRID, LAWS, WITNESSES, certify

def verify(name, law):
    field = RealIntervalField(160)
    sign, endpoints = WITNESSES[name]
    side_floor, radius_floor = FLOORS[name]
    starts, ends, errors = [], [], []
    for lower, upper in zip(endpoints[::2], endpoints[1::2]):
        a, b = QQ(lower) / GRID, QQ(upper) / GRID
        enclosed = []
        for u in (a, b):
            angle = field(sum(QQ(c) * u**(-n) for n, c in law.terms))
            signal = (u**2 * angle.cos(), u**2 * angle.sin())
            center = tuple(QQ(x.center().exact_rational()) for x in signal)
            enclosed.append((center, max(QQ(abs(x - y).upper().exact_rational())
                                          for x, y in zip(signal, center))))
        starts.append(vector(QQ, enclosed[0][0]))
        ends.append(vector(QQ, enclosed[1][0]))
        # Compute a rigorous interval for each component of the second
        # derivative before the unit complex rotation.
        parameter = field(a, b)
        slope = -sum(n * QQ(c) * parameter**(-n-1) for n, c in law.terms)
        second = sum(n * (n+1) * QQ(c) * parameter**(-n-2) for n, c in law.terms)
        real = 2 - (parameter*slope)**2
        imaginary = 4*parameter*slope + parameter**2*second
        curvature = QQ((abs(real)**2 + abs(imaginary)**2).sqrt().upper().exact_rational())
        errors.append(max(e[1] for e in enclosed) + (b-a)**2 * curvature / 8)
        assert a*a > QQ(radius_floor)

    def residual(points):
        a, b, c, d = points
        edge = b - a
        turn = sign * vector(QQ, (-edge[1], edge[0]))
        return vector(QQ, tuple(c-b-turn) + tuple(d-a-turn))

    base = residual(starts)
    cols = []
    for i, end in enumerate(ends):
        displaced = starts.copy()
        displaced[i] = end
        cols.append(residual(displaced) - base)
    system = matrix(QQ, cols).transpose()
    root = system.solve_right(-base)
    inverse_norm = max(sum(map(abs, row)) for row in system.inverse().rows())
    margin = min(min(u, 1-u) for u in root) / inverse_norm
    error = max(errors)
    assert 4*error < margin, (name, float(4*error), float(margin))
    chord = tuple(a + u*(b-a) for a, b, u in zip(starts, ends, root))
    side_allowance = sum(2*errors[i] + max(root[i], 1-root[i]) *
                         sum(map(abs, ends[i]-starts[i])) for i in (0, 1))
    side = chord[1] - chord[0]
    assert side.dot_product(side) > (QQ(side_floor) + side_allowance)**2
    assert certify(law, sign, endpoints, side_floor, radius_floor).certified


if __name__ == "__main__":
    for name, law in LAWS.items():
        verify(name, law)
        print(f"{name}: verified with SageMath interval arithmetic and rational matrices")
