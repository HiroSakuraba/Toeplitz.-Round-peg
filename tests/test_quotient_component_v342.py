import csv
from fractions import Fraction as Q
from importlib.util import find_spec
from io import StringIO
from pathlib import Path
import sys
import unittest
from collections import Counter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"code"))
from polygon_component_v341 import PolygonComponents, degree
from polygon_orbit_census_v223 import census_polygon, pt
from quotient_component_v342 import orbit_replay, recorded_stabilizers
from terrain_component_degree_v337 import TERRAIN

SQUARE = ((0,0),(2,0),(2,2),(0,2))


class OrbitTests(unittest.TestCase):
    @unittest.skipUnless(find_spec("sympy"), "symbolic germ replay needs SymPy")
    def test_square_monodromy_and_unmarked_count(self):
        cells = PolygonComponents(SQUARE).trace(tuple(pt(*p) for p in SQUARE))
        self.assertEqual(sum(sign for _,sign in degree(cells)[1]),0)
        for direction in (1,-1):
            r = orbit_replay(SQUARE,SQUARE,[(direction*y,direction*x) for x,y in SQUARE])
            self.assertEqual((r.stabilizer,r.marked_charts,r.chart_orbits),
                             ((0,1,2,3),60,15))
            self.assertEqual((r.parity,r.floor_squared,r.regular_charts),
                             (1,2,((0,1,2,3),)))

    def test_finite_rhombus_upper_witnesses(self):
        for t in (Q(1,4),Q(1,25),Q(-1,10)):
            polygon = [pt(x+t*y,y+t*x) for x,y in SQUARE]
            census = census_polygon(polygon)
            self.assertEqual(census.count,1)
            self.assertEqual(census.squares[0].side_squared,2*(1-t*t)**2)

    def test_sharp_redraws_stay_inside_error_budget(self):
        centered = ((-1,-1),(1,-1),(1,1),(-1,1))
        for eps in (Q(1,100),Q(1,2)):
            limit_squared = 2*(1-eps)**2
            for eta in (eps/2,eps/4):
                radius = 1-eps+eta/2
                strain = eta/(4*radius)
                polygon = [pt(radius*(x+strain*y),radius*(y+strain*x)) for x,y in centered]
                error = max(max(abs(p.x-x),abs(p.y-y)) for p,(x,y) in zip(polygon,centered))
                self.assertLessEqual(error,eps)
                census = census_polygon(polygon)
                self.assertEqual(census.count,1)
                squared = census.squares[0].side_squared
                self.assertGreaterEqual(squared,limit_squared)
                self.assertLess(squared,limit_squared+3*eta)

    @unittest.skipUnless(find_spec("sympy"), "symbolic germ replay needs SymPy")
    def test_zero_parity_and_nonregular_deformation(self):
        s = ((1,2),(2,2),(2,3),(1,3))
        v = ((1,2),(3,-1),(-2,1),(-1,3),(2,-1),(-4,-2))
        r = orbit_replay(TERRAIN,s,v)
        self.assertEqual((r.stabilizer,r.parity,r.regular_charts),((0,),0,()))
        with self.assertRaises(ValueError):
            orbit_replay(SQUARE,SQUARE,[(0,0)]*4)

    @unittest.skipUnless(find_spec("numpy"), "v336 corpus generator needs NumPy")
    def test_archived_covering_components(self):
        from terrain_replay_v341 import fixtures
        with (ROOT/"results/v341/component-certificates.csv").open() as f:
            records = {(r["kind"],int(r["seed"])):r for r in csv.DictReader(f)}
        self.assertEqual(recorded_stabilizers(fixtures(),records),
                         {("terrain",1):282,("union",1):194})


@unittest.skipUnless(find_spec("numpy"), "archived column geometry imports NumPy")
class FrameTests(unittest.TestCase):
    def test_registry_replay(self):
        from terrain_square_annulus_v342 import write_report
        s = StringIO()
        write_report(s)
        self.assertEqual(s.getvalue(),(ROOT/"results/v342/square-frames.csv").read_text())
        s.seek(0)
        self.assertEqual(Counter((r["kind"],r["at_target"]) for r in csv.DictReader(s)),
                         {("terrain","strict"):67,("terrain","equality"):65,
                          ("terrain","none"):150,("union","strict"):7,
                          ("union","equality"):6,("union","none"):181})

    def test_threshold_and_midpoint_clipping(self):
        from terrain_square_annulus_v342 import best_frame
        for width,height,outer,margin in ((2,2,1,Q(1,2)),(6,2,3,0)):
            frame = best_frame(0,width,dict.fromkeys(range(width),height),
                               dict.fromkeys(range(width),0))
            self.assertEqual((frame.radius,frame.outer_radius,frame.registration_margin),
                             (1,outer,margin))
        self.assertIsNone(best_frame(0,7,dict.fromkeys(range(7),2),dict.fromkeys(range(7),0)))
        h = [1,2,2,2]
        f = best_frame(0,4,dict(enumerate(h)),dict.fromkeys(range(4),0))
        self.assertEqual((f.radius,f.center,f.outer_radius),(1,(2,1),2))


if __name__ == "__main__":
    unittest.main()
