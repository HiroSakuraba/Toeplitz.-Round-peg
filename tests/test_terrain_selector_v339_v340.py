from collections import Counter
from fractions import Fraction as Q
from io import StringIO
from importlib.util import find_spec
from itertools import product
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "code"))
from terrain_selector_v339 import contact_crossings, end_component, sample, write_report
from terrain_component_degree_v337 import chart_cells, shear_witnesses
from polygon_orbit_census_v223 import _affine_system, _side_sq, pt


@unittest.skipUnless(find_spec("numpy"), "v338 sample replay requires NumPy")
class SampleTests(unittest.TestCase):
    def test_full_seed_report(self):
        buffer = StringIO()
        write_report(buffer)
        self.assertEqual(buffer.getvalue(), (ROOT / "results/v340/end-wall.csv").read_text())
        rows = list(sample())
        counts = Counter(c.method if c else "unresolved" for _, _, _, _, c in rows)
        self.assertEqual(counts, {"strict":131, "tall walls":42, "left plateau":10,
                                  "right plateau":11, "one corner":34, "unresolved":24})
        def smallest(sub):
            return min((sum(h), len(h), h, seed) for seed,h,i,b,c in sub)[2]
        self.assertEqual(smallest(r for r in rows if r[3] is None), (3,3))
        self.assertEqual(smallest(r for r in rows if r[3] is None and len(set(r[1])) > 1), (2,5))
        self.assertEqual(smallest(r for r in rows if r[-1] is None), (1,1,5,2))
        for _, _, i, _, c in rows:
            if c:
                self.assertGreaterEqual(2*c.side_squared, i*i)

    @unittest.skipUnless(find_spec("sympy"), "exact determinants require SymPy")
    def test_each_crossing_by_independent_incident_chart_determinants(self):
        import sympy as sp
        from robust_squares_v336 import terrain_polygon
        for seed, h, _, _, _ in sample():
            polygon = [pt(*p) for p in terrain_polygon(h)]
            edges = list(zip(polygon, polygon[1:]+polygon[:1]))
            for c in contact_crossings(h, permit_one_corner=True):
                a,k,b = c.contacts
                points = [(a+b-k,0), (b,b-k), (k,b-a), (a,k-a)]
                charts = []
                for x,y in points:
                    incident = [(p,q) for p,q in edges
                                if min(p.x,q.x) <= x <= max(p.x,q.x)
                                and min(p.y,q.y) <= y <= max(p.y,q.y)]
                    charts.append(incident)
                determinants = [sp.Matrix(_affine_system(es,1)[0]).det()
                                for es in product(*charts)]
                with self.subTest(seed=seed, contacts=c.contacts):
                    self.assertEqual(len(determinants), 1 if c.method == "strict" else 2)
                    self.assertIn({sp.sign(d) for d in determinants}, ({1},{-1}))


class GeometryTests(unittest.TestCase):
    def test_plateau_orientation_and_bounds(self):
        self.assertEqual(end_component((4,8,1,7)).side_squared, 10)
        self.assertEqual(end_component((7,1,8,4)).contacts, (3,))
        self.assertEqual(end_component((4,4,4,7)).side_squared, 8)
        self.assertIsNone(end_component((4,3,8,7)))  # first step descends
        self.assertIsNone(end_component((4,4,4,4)))  # no rising step
        self.assertEqual(end_component((7,1,1,5)).side_squared, 16)

    def test_first_corner_and_parity_exceptions(self):
        for h, contacts, expected in (
            ((2,5,2),(0,2,3),5), ((8,4,3,1,1,6),(0,1,6),26),
            ((5,1,7,4,6,4),(0,3,6),18), ((4,3,1,6,3,8),(2,4,6),8)
        ):
            found = {c.contacts:c for c in contact_crossings(h, permit_one_corner=True)}
            self.assertEqual(found[contacts].side_squared, expected)

    def test_tapered_rectangle(self):
        t = Q(1,100)
        polygon = [(t*y+(1-2*t*y)*x,y) for x,y in ((0,0),(1,0),(1,2),(0,2))]
        # Independent full rational census of a finite taper.
        cells = dict(chart_cells(polygon))
        self.assertEqual(len(cells), 3)
        sizes = []
        for corners in cells.values():
            square, = corners
            sizes.append(_side_sq(square))
        self.assertCountEqual(sizes, [1/(1+2*t)**2] +
                              [(1+t*t)/(1+t)**4]*2)

    @unittest.skipUnless(find_spec("sympy"), "one-sided shear replay requires SymPy")
    def test_all_plateau_shear_regimes(self):
        import sympy as sp
        t = sp.Symbol('t', positive=True)
        candidates = [tuple(map(int, s)) for s in
                      ('0014','0015','0114','0115','0135','0145','0155','1145','1155')]
        for a in (Q(1,4),Q(1,2),Q(3,4)):
            polygon = ((0,0),(1,0),(1,3),(a,3),(a,1),(0,1))
            found = shear_witnesses(polygon, candidates, 1)
            self.assertEqual(set(found), {(0,1,3,5) if a < Q(1,2) else (0,1,4,5)})
            (_, side_squared), = found.values()
            u = min(a,Q(1,2))
            self.assertEqual(sp.limit(side_squared,t,0), u*u+(1-u)**2)


if __name__ == "__main__":
    unittest.main()
