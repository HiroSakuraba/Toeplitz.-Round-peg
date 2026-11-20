"""Independent finite regression suite for the restricted v269-v273 models."""
from fractions import Fraction as Q
from math import pi
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from relative_rectangle_v269_v273 import (
    _det, square_residual, ellipse_square_jacobian, ellipse_relative_certificate,
    relative_interval_parity, cyclic_dot_character, ellipse_angle_certificate,
    ellipse_angle_band_floors, permutation_sector_ledger,
    cyclic_square_residual_equivariance, dihedral_sector_ledger, rectangle_boundary_stratum,
    is_simple_polygon, notch_square_certificate, spectral_axiom_countermodel,
)


class V269RelativeIndex(unittest.TestCase):
    def test_ellipse_relative_end_signs(self):
        d=ellipse_relative_certificate(3,2)
        self.assertEqual(d["flat_end_dot_signs"],(-36,16))
        self.assertEqual(relative_interval_parity(*d["flat_end_dot_signs"]),1)

    def test_ellipse_square_exact(self):
        d=ellipse_relative_certificate(3,4)
        self.assertEqual(d["square_side_squared"],Q(576,25))
        self.assertEqual(d["square_circumradius_squared"],Q(288,25))

    def test_ellipse_four_parameter_jacobian_exact(self):
        self.assertEqual(ellipse_square_jacobian(3,4), Q(4032,25))
        self.assertEqual(ellipse_square_jacobian(4,3), -Q(4032,25))
        self.assertEqual(ellipse_square_jacobian(3,2), -Q(4*9*4*(9-4),13))

    def test_ellipse_jacobian_nonzero_for_several_nonpythagorean_axes(self):
        for a,b in ((2,3),(1,7),(Q(3,2), Q(7,3)),(5,12)):
            with self.subTest(a=a,b=b):
                self.assertNotEqual(ellipse_square_jacobian(a,b),0)

    def test_circle_rank_loss_rejected(self):
        with self.assertRaises(ValueError): ellipse_square_jacobian(3,3)
        with self.assertRaises(ValueError): ellipse_relative_certificate(-3,2)

    def test_undefined_boundary_parity_rejected(self):
        for a,b in ((0,1),(1,0),(0,0)):
            with self.assertRaises(ValueError): relative_interval_parity(a,b)

    def test_same_sign_even_parity(self):
        self.assertEqual(relative_interval_parity(-2,-3),0)
        self.assertEqual(relative_interval_parity(2,3),0)
        self.assertEqual(relative_interval_parity(2,-3),1)

    def test_tensor_character(self):
        self.assertEqual([cyclic_dot_character(k) for k in range(4)],[1,-1,1,-1])


class V270SpectralAxiomObstruction(unittest.TestCase):
    def test_dimensionless_analytical_boundary_and_derivative(self):
        r=spectral_axiom_countermodel(1)
        self.assertEqual(r["midpoint_derivative_over_radius_squared"],0)
        self.assertLess(r["maximum_derivative_over_radius_squared"],1)
        self.assertAlmostEqual(r["normalized_action"](0),0)
        self.assertAlmostEqual(r["normalized_action"](pi),1)
        self.assertAlmostEqual(r["derivative_normalized"](pi/2),0,places=30)

    def test_triangle_gap_nonnegative_regressions(self):
        r=spectral_axiom_countermodel(1)
        for i in range(0,25):
            for j in range(0,25-i):
                x,y=i*pi/24,j*pi/24
                self.assertGreaterEqual(r["normalized_triangle_gap"](x,y),-1e-15)

    def test_lipschitz_sharp_cutoff(self):
        r=spectral_axiom_countermodel(pi/2)
        self.assertAlmostEqual(r["maximum_derivative_over_radius_squared"],1)
        with self.assertRaises(ValueError): spectral_axiom_countermodel(pi/2+Q(1,100))

    def test_action_strictly_rises_across_90_degrees(self):
        r=spectral_axiom_countermodel(1)
        a=r["normalized_action"]
        self.assertLess(a(pi/2-0.01),a(pi/2))
        self.assertLess(a(pi/2),a(pi/2+0.01))


class V271SixSectorAudit(unittest.TestCase):
    def test_all_24_labellings_six_cyclic_sectors(self):
        d=permutation_sector_ledger()
        self.assertEqual(d["six_sectors"],6)
        self.assertEqual(d["labellings"],24)
        self.assertEqual(d["orbit_sizes"],(4,)*6)
        self.assertEqual(d["pairing_counts"],{0:8,1:8,2:8})

    def test_unequal_rectangle_diagonal_pairing_is_unique(self):
        d=permutation_sector_ledger(3,4)
        self.assertEqual(len(set(x["ordering"] for x in d["rows"])),24)
        self.assertEqual(len(set(x["cyclic_representative"] for x in d["rows"])),6)

    def test_square_does_not_create_extra_diagonal_pairings(self):
        self.assertEqual(permutation_sector_ledger(5,5)["pairing_counts"],{0:8,1:8,2:8})

    def test_six_cyclic_sectors_become_three_dihedral_classes(self):
        self.assertEqual(dihedral_sector_ledger()["eightfold_orbit_sizes"],(8,8,8))

    def test_target_representation_on_exact_eight_basis(self):
        d=cyclic_square_residual_equivariance()
        self.assertEqual(d["basis_checks"],8)
        self.assertEqual(d["target_determinant"],1)
        self.assertEqual(d["source_determinant"],-1)

    def test_simple_permutation_determinants(self):
        self.assertEqual(_det(((0,1),(1,0))),-1)
        self.assertEqual(_det(((1,0),(0,1))),1)
        with self.assertRaises(ValueError): _det(((1,0,2),(2,3,4)))


class V272DegenerationAudit(unittest.TestCase):
    def test_exact_rectangle_limit_strata(self):
        self.assertEqual(rectangle_boundary_stratum(1,1),"genuine_noncollapsed_rectangle")
        self.assertEqual(rectangle_boundary_stratum(0,1),"two_pair_flattening_nonzero_diagonal")
        self.assertEqual(rectangle_boundary_stratum(1,0),"two_pair_flattening_nonzero_diagonal")
        self.assertEqual(rectangle_boundary_stratum(0,0),"total_four_corner_collision")

    def test_inadmissible_length_rejected(self):
        with self.assertRaises(ValueError): rectangle_boundary_stratum(-1,1)

    def test_bump_polygons_are_actual_jordan_curves(self):
        for n in (4,5,10,100,1000):
            with self.subTest(n=n):
                d=notch_square_certificate(n)
                self.assertTrue(d["simple_Jordan_polygon"])
                self.assertTrue(is_simple_polygon(d["polygon"]))
                self.assertEqual(square_residual(d["selected_square"]),(0,0,0,0))
                self.assertEqual(d["side_squared"],Q(1,n*n))
                self.assertEqual(d["circumradius_squared"],Q(1,2*n*n))
                self.assertEqual(d["hausdorff_distance_to_unit_square_upper_bound"],Q(1,n))
                self.assertEqual(d["normalized_square"],((1,1),(0,1),(0,0),(1,0)))

    def test_selected_shrinkout_monotonic(self):
        self.assertLess(notch_square_certificate(100)["side_squared"],
                        notch_square_certificate(10)["side_squared"])

    def test_adversarial_non_Jordan_polygon_rejected(self):
        self.assertFalse(is_simple_polygon(((0,0),(2,2),(0,2),(2,0))))

    def test_invalid_bump_parameters_rejected(self):
        with self.assertRaises(ValueError):notch_square_certificate(3)
        with self.assertRaises(ValueError):notch_square_certificate(4,Q(9,10))


class V273AngleContinuation(unittest.TestCase):
    def test_all_ellipse_angle_cosines_are_realizable(self):
        for c in (Q(-9,10),Q(-1,2),Q(0),Q(1,2),Q(9,10)):
            with self.subTest(c=c):
                d=ellipse_angle_certificate(3,4,c)
                self.assertEqual(d["cos_diagonal_angle"],c)
                self.assertGreater(d["x_squared"],0)
                self.assertGreater(d["y_squared"],0)
                self.assertEqual((d["y_squared"]-d["x_squared"])/
                                 d["circumradius_squared"],c)

    def test_square_angle_and_positive_size(self):
        d=ellipse_angle_certificate(3,4,0)
        self.assertEqual(d["x_squared"],Q(144,25))
        self.assertEqual(d["y_squared"],Q(144,25))
        self.assertEqual(d["side_squared"],(Q(576,25),Q(576,25)))

    def test_exact_radius_and_side_floor_on_angle_band(self):
        d=ellipse_angle_band_floors(3,4,Q(1,2))
        self.assertEqual(d["circumradius_squared_floor"],Q(192,19))
        for c in (Q(-1,2),Q(-1,4),Q(0),Q(1,4),Q(1,2)):
            v=ellipse_angle_certificate(3,4,c)
            self.assertGreaterEqual(v["circumradius_squared"],d["circumradius_squared_floor"])
            self.assertGreaterEqual(min(v["side_squared"]),d["minimum_side_squared_floor"])

    def test_angle_floor_sharp_endpoint(self):
        d=ellipse_angle_band_floors(3,4,Q(1,2))
        r=ellipse_angle_certificate(3,4,Q(-1,2))
        self.assertEqual(r["circumradius_squared"],d["circumradius_squared_floor"])
        self.assertEqual(min(r["side_squared"]),d["minimum_side_squared_floor"])

    def test_invalid_angles_are_rejected(self):
        for c in (-1,1,2):
            with self.assertRaises(ValueError):ellipse_angle_certificate(3,4,c)
        for eta in (-Q(1,2),Q(1)):
            with self.assertRaises(ValueError):ellipse_angle_band_floors(3,4,eta)
        with self.assertRaises(ValueError):ellipse_angle_certificate(3,3,0)

if __name__ == "__main__": unittest.main()