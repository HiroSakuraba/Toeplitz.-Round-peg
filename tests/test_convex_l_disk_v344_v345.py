from fractions import Fraction as Q
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from convex_l_disk_v344_v345 import OrthogonalL, convex_certifies, convex_root_data
from polygon_component_v341 import PolygonComponents, degree, segment_floor
from polygon_orbit_census_v223 import pt


class DiskBoundTests(unittest.TestCase):
    def test_convex_finite_cutoff(self):
        self.assertEqual(convex_root_data(1, 2, Q(1, 100)), (Q(99, 100), Q(8193, 10000)))
        self.assertTrue(convex_certifies(1, 2, Q(1, 100), Q(9, 10)))
        self.assertFalse(convex_certifies(1, 2, Q(1, 100), Q(99, 100)))
        self.assertEqual(convex_root_data(3, 3, 0), (3, 9))
        self.assertFalse(convex_certifies(3, 3, 0, 3))
        # q is exactly a root. The cutoff test needs strict inequality.
        self.assertFalse(convex_certifies(13, 15, 1, 4))
        self.assertFalse(convex_certifies(1, 1, 2, Q(1, 2)))
        with self.assertRaises(ValueError):
            convex_certifies(2, 1, 0, 1)

    def test_marked_component_against_independent_polygon_degree(self):
        # Include each active face, two-way ties in both orientations,
        # a sliding component, and the triple tie. Use nonintegral contacts too.
        for dims in ((9, 7, 1, 1), (4, 7, 2, 3), (7, 4, 3, 2),
                     (5, 7, 2, 3), (7, 5, 3, 2), (5, 5, 3, 4),
                     (5, 5, 2, 3), (5, 5, 2, 4), (7, 7, 2, 2),
                     (Q(9, 2), 6, Q(7, 4), Q(13, 4))):
            with self.subTest(dims=dims):
                shape = OrthogonalL(*dims)
                cert = shape.certificate()
                u = sum(cert.u_interval) / 2
                seed = tuple(pt(*p) for p in cert.elbow_square(u))
                cells = PolygonComponents(shape.polygon()).trace(seed)
                self.assertEqual(min(segment_floor(c.squares) for c in cells.values()),
                                 cert.side_squared)
                self.assertEqual(abs(sum(s for _, s in degree(cells)[1])), 1)

    def test_arms_and_inradius_transition(self):
        for shape, method, floor in (
            (OrthogonalL(2, 7, 1, 2), "end plateau", 2),
            (OrthogonalL(7, 2, 2, 1), "end plateau", 2),
            (OrthogonalL(4, 7, 3, 4), "end plateau", 10),
            (OrthogonalL(4, 7, 3, 5), "tall walls", 16),
            (OrthogonalL(5, 5, 3, 4), "elbow", Q(25, 2)),
        ):
            cert = shape.certificate()
            self.assertEqual((cert.method, cert.side_squared), (method, floor))
        equal_arms = OrthogonalL(8, 9, 1, 1)
        self.assertTrue(equal_arms.fits_disk(Q(5857, 10000)))
        self.assertFalse(equal_arms.fits_disk(Q(5858, 10000)))
        for b, r in ((2, 1), (3, Q(3, 2))):
            shape = OrthogonalL(8, 9, 1, b)
            self.assertTrue(shape.fits_disk(r))
            self.assertFalse(shape.fits_disk(r + Q(1, 10000)))
            self.assertTrue(shape.reaches_disk_target(2 * r**2))
            self.assertFalse(shape.reaches_disk_target(2 * r**2 - Q(1, 10000)))

    def test_exact_disk_comparison_and_symmetries(self):
        for W in range(2, 8):
            for H in range(2, 8):
                for a in range(1, W):
                    for b in range(1, H):
                        shape = OrthogonalL(W, H, a, b)
                        floor = shape.certificate().side_squared
                        self.assertTrue(shape.reaches_disk_target(floor))
                        self.assertEqual(OrthogonalL(H, W, b, a).certificate().side_squared,
                                         floor)
                        for numerator in range(1, 13):
                            r = Q(numerator, 4)
                            # Compare the algebraic formula to disk feasibility.
                            if not shape.fits_disk(r):
                                self.assertTrue(shape.reaches_disk_target(2*r*r))
                            if shape.fits_disk(r + Q(1, 100000)):
                                self.assertFalse(shape.reaches_disk_target(2*r*r))
        scale = Q(5, 11)
        dims = (5, 5, 3, 4)
        self.assertEqual(OrthogonalL(*(scale*x for x in dims)).certificate().side_squared,
                         scale**2 * OrthogonalL(*dims).certificate().side_squared)


if __name__ == "__main__":
    unittest.main()
