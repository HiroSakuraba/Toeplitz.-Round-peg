"""Exact-rational v217 tests of boundary-chain additivity and failure modes."""
import random
import sys
import unittest
from fractions import Fraction as Q
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'code'))

from multiscale_degree_v217 import (
    Point, Arc, pt, arc, _at, subarc, partition_arcs,
    planted_example, cancellation_example, multiscale_ledger,
    opposite_gap_squared, boundary_has_square, ray_flux,
    _rank, _det, _feasible_box, Inconclusive, BoundaryRoot,
)


class MultiscaleDegreeTests(unittest.TestCase):
    def test_exact_box_feasibility_including_singular(self):
        self.assertEqual(_rank(((Q(1),Q(2)),(Q(2),Q(4)))),1)
        self.assertEqual(_det(((Q(1),Q(2)),(Q(3),Q(4)))),Q(-2))
        self.assertTrue(_feasible_box(((Q(1),Q(1)),(Q(2),Q(2))), (Q(1),Q(2)),2))
        self.assertFalse(_feasible_box(((Q(1),Q(1)),(Q(2),Q(2))), (Q(1),Q(3)),2))

    def test_exact_restriction_across_internal_vertices(self):
        a=arc((0,0),(Q(1,2),1),(Q(3,4),Q(1,2)),(1,0))
        self.assertEqual(_at(a,Q(1,3)),pt(Q(1,2),1))
        self.assertEqual(_at(a,Q(1,2)),pt(Q(5,8),Q(3,4)))
        x=subarc(a,Q(1,6),Q(5,6))
        self.assertEqual(x.vertices,(pt(Q(1,4),Q(1,2)),pt(Q(1,2),1),
                                     pt(Q(3,4),Q(1,2)),pt(Q(7,8),Q(1,4))))

    def test_restriction_of_full_arc_is_identical(self):
        a=cancellation_example()[2]
        self.assertEqual(subarc(a,0,1),a)

    def test_clean_split_tracks_nonzero_degree(self):
        l=multiscale_ledger(planted_example(),((Q(1,3),),(),(),()))
        self.assertEqual(l.parent_degree,-1)
        self.assertEqual(l.sum_child_flux,-1)
        self.assertEqual([v for _,v,_ in l.children],[0,-1])
        self.assertTrue(all(ok for _,_,ok in l.children))
        self.assertEqual(l.internal_wall_pairs,1)
        self.assertGreater(l.opposite_distance_squared,0)

    def test_internal_root_at_cut_does_not_change_total(self):
        l=multiscale_ledger(planted_example(),((Q(1,2),),(),(),()))
        self.assertEqual((l.parent_degree,l.sum_child_flux),(-1,-1))
        self.assertEqual(l.boundary_root_children,2)
        self.assertTrue(all(not ok for _,_,ok in l.children))
        self.assertEqual(l.internal_wall_pairs,1)

    def test_internal_wall_flux_can_depend_on_ray(self):
        # An individual child with a boundary root has no intrinsic degree.
        # Its algebraic positive-ray flux can change with ray direction,
        # while the sum over both children remains the parent degree.
        parts,_=partition_arcs(planted_example(),((Q(1,2),),(),(),()))
        children=[tuple(parts[i][j] for i,j in enumerate(idx))
                  for idx in ((0,0,0,0),(1,0,0,0))]
        a=(-2,-2,-2,-2)
        b=(-2,-2,-2,1)
        self.assertEqual([ray_flux(child,a).total for child in children],[-1,0])
        self.assertEqual([ray_flux(child,b).total for child in children],[0,-1])
        self.assertEqual(ray_flux(planted_example(),a).total,-1)
        self.assertEqual(ray_flux(planted_example(),b).total,-1)

    def test_opposite_sign_squares_split_into_two_certified_degrees(self):
        l=multiscale_ledger(cancellation_example(),((),(),(Q(1,2),),()))
        self.assertEqual((l.parent_degree,l.sum_child_flux),(0,0))
        self.assertEqual([v for _,v,_ in l.children],[-1,1])
        self.assertTrue(all(ok for _,_,ok in l.children))

    def test_two_dimension_partition_telescopes(self):
        l=multiscale_ledger(planted_example(),((Q(1,3),),(Q(2,5),),(),()))
        self.assertEqual(len(l.children),4)
        self.assertEqual(l.internal_wall_pairs,4)
        self.assertEqual(l.parent_degree,-1)
        self.assertEqual(l.sum_child_flux,-1)
        self.assertTrue(all(ok for _,_,ok in l.children))
        self.assertEqual(sum(v!=0 for _,v,_ in l.children),1)

    def test_refinement_of_all_four_axes(self):
        l=multiscale_ledger(planted_example(),((Q(1,3),),(Q(1,3),),
                                                   (Q(1,3),),(Q(1,3),)))
        self.assertEqual(len(l.children),16)
        self.assertEqual(l.internal_wall_pairs,32)
        self.assertEqual(l.parent_degree,l.sum_child_flux)
        self.assertEqual(sum(v!=0 for _,v,_ in l.children),1)
        self.assertTrue(all(ok for _,_,ok in l.children))

    def test_no_refinement_matches_parent(self):
        l=multiscale_ledger(planted_example())
        self.assertEqual(l.children,(((0,0,0,0),-1,True),))
        self.assertEqual(l.internal_wall_pairs,0)

    def test_opposite_orientation_has_zero_degree(self):
        l=multiscale_ledger(planted_example(),((Q(1,3),),(),(),()),orientation=-1)
        self.assertEqual(l.parent_degree,0)
        self.assertEqual(l.sum_child_flux,0)

    def test_reversal_changes_degree_and_preserves_additivity(self):
        a=list(planted_example())
        a[0]=Arc(tuple(reversed(a[0].vertices)))
        l=multiscale_ledger(tuple(a),((Q(1,3),),(),(),()))
        self.assertEqual(l.parent_degree,1)
        self.assertEqual(l.sum_child_flux,1)

    def test_global_translation_does_not_change_degree(self):
        a=planted_example()
        translated=tuple(Arc(tuple(pt(p.x+Q(11,7),p.y-Q(3,11))
                                     for p in z.vertices)) for z in a)
        l=multiscale_ledger(translated,((Q(1,3),),(),(),()))
        self.assertEqual(l.parent_degree,-1)

    def test_fixed_deck_shifts_do_not_alias(self):
        a=planted_example()
        wrong=multiscale_ledger(a,((Q(1,3),),(),(),()),shifts=(0,2,0,0))
        self.assertEqual(wrong.parent_degree,0)
        self.assertEqual(wrong.sum_child_flux,0)

    def test_parent_external_boundary_square_refused(self):
        a=list(planted_example())
        a[0]=subarc(a[0],Q(1,2),Q(1))
        self.assertTrue(boundary_has_square(a))
        with self.assertRaises(BoundaryRoot):
            multiscale_ledger(tuple(a))

    def test_zero_side_guard(self):
        a=list(planted_example())
        a[2]=a[0]
        self.assertEqual(opposite_gap_squared(a),0)
        with self.assertRaises(ValueError):
            multiscale_ledger(tuple(a))

    def test_invalid_cuts(self):
        with self.assertRaises(ValueError):
            partition_arcs(planted_example(),((Q(1,2),Q(1,3)),(),(),()))
        with self.assertRaises(ValueError):
            subarc(planted_example()[0],Q(1,2),Q(1,2))
        with self.assertRaises(ValueError):
            partition_arcs(planted_example(),((Q(-1,3),),(),(),()))

    def test_no_ray_search_reports_inconclusive(self):
        with self.assertRaises(Inconclusive):
            multiscale_ledger(planted_example(),max_rays=0)

    def test_random_rational_clean_cuts_and_perturbations(self):
        rng=random.Random(217217)
        original=planted_example()
        for _ in range(10):
            perturbed=[]
            for a in original:
                dx=Q(rng.randint(-1,1),4096)
                dy=Q(rng.randint(-1,1),4096)
                perturbed.append(Arc(tuple(pt(p.x+dx,p.y+dy) for p in a.vertices)))
            i=rng.randrange(4)
            t=Q(rng.choice((3,4,5,6,7,8,9,10,11,12,13,14)),17)
            cuts=[(),(),(),()];cuts[i]=(t,)
            l=multiscale_ledger(tuple(perturbed),tuple(cuts))
            self.assertEqual((l.parent_degree,l.sum_child_flux),(-1,-1))
            self.assertEqual(l.boundary_root_children,0)
            self.assertEqual(l.internal_wall_pairs,1)

    def test_crosscheck_v216_degree_when_importable(self):
        try:
            from boundary_degree_v216 import Arc as OldArc, pt as old_pt, boundary_degree
        except ImportError:
            self.skipTest('v216 module not present in local standalone staging')
        for a in (planted_example(),cancellation_example()):
            converted=tuple(OldArc(tuple(old_pt(p.x,p.y) for p in x.vertices)) for x in a)
            old=boundary_degree(converted)
            new=multiscale_ledger(a,((Q(1,3),),(),(),()))
            self.assertEqual(new.parent_degree,old.degree)


if __name__=='__main__':
    unittest.main()