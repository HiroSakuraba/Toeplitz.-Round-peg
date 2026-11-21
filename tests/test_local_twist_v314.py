from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"code"))
from local_twist_v314 import cap_clearance, center_charts, local_certificate, pt, radial_polygon, twist
from half_turn_degree_v313 import cut_crossings, polygon as ear_polygon, tube_certificate


class LocalParityTests(unittest.TestCase):
    def test_circle_seam_inside_and_outside(self):
        box=(pt(-3,-2),pt(3,-2),pt(3,2),pt(-3,2))
        for p in (box,box[::-1],box[1:]+box[:1]):
            self.assertEqual(twist(p,Q(1))["twist"],1)
            self.assertEqual(twist(p,Q(5))["twist"],0)

    def test_radial_order_and_transversality(self):
        data=local_certificate()
        self.assertEqual(data["charts"],12)
        self.assertEqual(data["radii_squared"],["457/25","2088/49","2169/49"])
        self.assertEqual([r["twist"] for r in data["circle_scan"]],[1,0,1,0])
        for (lo,hi),(_,_,r2) in zip(data["radius_squared_bands"],data["sample_centers_and_radii_squared"]):
            self.assertLess(Q(lo),Q(r2))
            self.assertLess(Q(r2),Q(hi))

    def test_cutoff_contact_and_asymmetry(self):
        diamond=(pt(5,0),pt(0,5),pt(-5,0),pt(0,-5))
        with self.assertRaisesRegex(ValueError,"unresolved"):
            twist(diamond,Q(5),max_depth=8)
        with self.assertRaisesRegex(ValueError,"half-turn"):
            twist((pt(5,0),pt(0,4),pt(-5,0),pt(0,-5)),Q(1))
        with self.assertRaises(ValueError):
            center_charts(radial_polygon(),0,Q(1,10))

    def test_split_model_and_clearance(self):
        self.assertEqual(cut_crossings()["one_ear_moved_degrees"],[1,0])
        result=tube_certificate(ear_polygon(Q(1),Q(9,10)),epsilon=Q(1,8192))
        self.assertEqual(result["side_interval"],["2/5","2/5"])
        self.assertEqual(sum(result["excluded_leaves"].values()),2164)

    def test_cap_budget(self):
        certificate=cap_clearance()
        self.assertEqual(certificate["side_floor_squared"],"961/512")
        self.assertGreater(min(map(Q,certificate["margins"].values())),0)
        with self.assertRaises(ValueError):
            cap_clearance(Q(1,8))

    def test_registration_budget_fails_closed(self):
        with self.assertRaises(ArithmeticError):
            local_certificate(Q(1,8),Q(1,64))
