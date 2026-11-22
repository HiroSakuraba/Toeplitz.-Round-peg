import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from spectral_noncollapse_v318_v322 import (
    audit, plateau, angular_floor, central_mean, subdivision_vertices
)


class SpectralNegativeControls(unittest.TestCase):
    def test_all_affine_cells(self):
        for den in range(2, 25):
            for num in range(1, den+1):
                for mixing in (Q(0), Q(1,1000), Q(1,3), Q(1)):
                    row = audit(Q(num,den), mixing)
                    self.assertEqual(row["min_subadditivity_slack"], 0)

    def test_saturation_all_windows(self):
        for b in (Q(1,100), Q(1,2), Q(9,10), Q(999,1000), Q(1)):
            for h in (Q(1,10000), Q(1,100), Q(1,10), Q(1,4), Q(1,2)):
                self.assertEqual(central_mean(b,h), angular_floor(b,h))
                self.assertGreaterEqual(central_mean(b,h,Q(1,9)), angular_floor(b,h))

    def test_strict_approximants(self):
        b = Q(9,10)
        for d in (10,100,1000,10000):
            row = audit(b,Q(1,d))
            self.assertTrue(all(s > 0 for s in row["slopes"]))
            self.assertEqual(row["central_slope"], b/d)
            self.assertEqual(plateau(b,Q(1,2),Q(1,d)), b/2)

    def test_plateau_even_near_isodiametric_equality(self):
        b = Q(999999,1000000)
        self.assertEqual(central_mean(b,(1-b)/2), 0)
        self.assertGreater(central_mean(b,1-b), 0)

    def test_exact_equality_case(self):
        self.assertEqual(audit(Q(1))["slopes"], [Q(1),Q(1)])
        self.assertEqual(central_mean(1,Q(1,100)), 1)

    def test_triangle_vertices_present(self):
        self.assertTrue({(0,0),(0,1),(1,0)} <= set(subdivision_vertices(Q(2,3))))

    def test_domain_guards(self):
        for b,x in ((0,0),(2,0),(1,-1),(1,2)):
            with self.assertRaises(ValueError):
                plateau(b,x)
        with self.assertRaises(ValueError):
            angular_floor(1,0)
        with self.assertRaises(ValueError):
            central_mean(1,1)


if __name__ == "__main__":
    unittest.main()
