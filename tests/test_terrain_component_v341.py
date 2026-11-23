import csv
from fractions import Fraction as Q
from importlib.util import find_spec
from io import StringIO
from pathlib import Path
import unittest
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT / "code"))
from polygon_component_v341 import PolygonComponents, degree, determinant, segment_floor
from polygon_orbit_census_v223 import _affine_system, _side_sq, _solve, mul, pt, sub
from terrain_component_degree_v337 import TERRAIN, GLUED, chart_cells
from terrain_ledge_v341 import ledge_squares

SMALL = (
    (((0,0),(4,0),(4,2),(3,2),(3,5),(2,5),(2,1),(0,1)),
     ((1,1),(3,0),(4,2),(2,3))),
    (((0,0),(4,0),(4,2),(2,2),(2,4),(1,4),(1,1),(0,1)),
     ((0,1),(2,0),(3,2),(1,3))),
)


class ComponentTests(unittest.TestCase):
    def test_tied_minima_have_sharp_bounds_and_room_for_detail(self):
        for polygon,coords in SMALL:
            seed = tuple(pt(*v) for v in coords)
            # Full unperturbed orbit census gives the independent upper bound.
            self.assertEqual(max(_side_sq(q) for _,qs in chart_cells(polygon) for q in qs),5)
            engine = PolygonComponents(polygon)
            cells = engine.trace(seed)
            self.assertEqual(len(cells),2)
            self.assertEqual(min(segment_floor(c.squares) for c in cells.values()),5)
            self.assertEqual(sum(s for _,s in degree(cells)[1]),1)
            self.assertTrue(any(set(q)==set(coords) for q,*_ in ledge_squares(polygon)))
            # Recompute matrices in unit-speed coordinates and check the local box.
            for ids in cells:
                tangents = []
                for j,p in zip(ids,seed):
                    a,b = engine.edges[j]
                    vec = sub(b,a)
                    length = max(abs(vec.x),abs(vec.y))
                    for end in (a,b):
                        if end != p:
                            self.assertGreaterEqual(max(abs(p.x-end.x),abs(p.y-end.y)),1)
                    tangents.append((pt(0,0),mul(vec,1/length)))
                mat,_ = _affine_system(tangents,1)
                self.assertEqual(determinant(mat),1)
                inverse_cols = [_solve(mat,[int(j==k) for j in range(4)]) for k in range(4)]
                self.assertEqual(max(sum(map(abs,row)) for row in zip(*inverse_cols)),3)

    def test_v337_independent_shear_fixtures(self):
        for poly,floor in ((TERRAIN,Q(2)),(GLUED,Q(49,2))):
            engine = PolygonComponents(poly)
            seed = engine.lattice_seeds()[0]
            for shift in (0,1):
                labeled = seed[shift:]+seed[:shift]
                cells = engine.trace(labeled)
                self.assertEqual(min(segment_floor(c.squares) for c in cells.values()),floor)
                self.assertEqual(abs(sum(s for _,s in degree(cells)[1])),1)

    def test_single_square_boundary_has_zero_labeled_degree(self):
        square = ((0,0),(2,0),(2,2),(0,2))
        cells = PolygonComponents(square).trace(tuple(pt(*v) for v in square))
        self.assertEqual(len(cells),60)
        self.assertEqual(min(segment_floor(c.squares) for c in cells.values()),2)
        _,hits = degree(cells)
        self.assertTrue(hits)
        self.assertEqual(sum(s for _,s in hits),0)

    def test_bad_seed_rejected(self):
        engine = PolygonComponents(TERRAIN)
        for seed in (((0,0),)*4, ((3,0),(4,0),(4,1),(3,1))):
            with self.assertRaises(ValueError):
                engine.trace(tuple(pt(*v) for v in seed))

    @unittest.skipUnless(find_spec("numpy"), "archived terrain constructor imports NumPy")
    def test_ledge_hits_by_separate_component_trace(self):
        from robust_squares_v336 import terrain_polygon
        with (ROOT / "results/v340/end-wall.csv").open() as f:
            previous = list(csv.DictReader(f))
        for row in previous:
            if row["method"] != "unresolved":
                continue
            polygon = terrain_polygon(tuple(map(int,row["heights"].split())))
            engine = PolygonComponents(polygon)
            for square,*_ in ledge_squares(polygon):
                seed = tuple(pt(*p) for p in square)
                if 2*_side_sq(seed) < int(row["interior_side"])**2:
                    continue
                cells = engine.trace(seed)
                self.assertEqual(min(segment_floor(c.squares) for c in cells.values()),_side_sq(seed))
                self.assertEqual(abs(sum(d for _,d in degree(cells)[1])),1)

    @unittest.skipUnless(find_spec("numpy"), "v336 corpus generator imports NumPy")
    def test_all_certificates_replay(self):
        from terrain_replay_v341 import write_report
        stream = StringIO()
        write_report(stream)
        self.assertEqual(stream.getvalue(),(ROOT/"results/v341/component-certificates.csv").read_text())
        stream.seek(0)
        rows = list(csv.DictReader(stream))
        self.assertEqual(len({(r['kind'],r['seed']) for r in rows}),476)
        self.assertEqual(sum(r['kind']=='terrain' for r in rows),282)
        self.assertEqual(sum(r['v340']=='unresolved' for r in rows),24)
        self.assertEqual(sum(bool(r['ledge_floor2']) for r in rows),11)
        for row in rows:
            self.assertIn(int(row['degree']),(-1,1))
            self.assertGreaterEqual(2*Q(row['floor2']),int(row['i'])**2)


if __name__ == '__main__':
    unittest.main()
