"""v232: multi-candidate fail-closed search and exact C0 phase-neighborhood proofs."""
import sys
import unittest
from pathlib import Path
from fractions import Fraction as Q
from unittest.mock import patch
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))

from spiral_multicandidate_v232 import (
    CandidateScan, PhasePerturbationCertificate, _source_signature,
    scan_candidates, certify_candidates, adaptive_multicandidate,
    certify_phase_neighborhood, report,
)
from spiral_pruned_v231 import STORED_WITNESSES, STORED_STRONG_FLOORS, replay_stored_witness


class PhaseRobustnessTests(unittest.TestCase):
    def test_three_phase_neighborhoods_are_exactly_valid(self):
        for p in (Q(2),Q(5,2),Q(3)):
            cert=certify_phase_neighborhood(p,Q(1,10**6))
            self.assertTrue(cert.exact_degree_valid,(p,cert))
            self.assertLess(cert.four_E_total,cert.margin)
            self.assertEqual((cert.side_strict_floor,cert.radius_strict_floor),
                             STORED_STRONG_FLOORS[p])

    def test_unperturbed_certificates_retain_strength(self):
        for p in STORED_WITNESSES:
            cert=certify_phase_neighborhood(p,Q(0))
            self.assertTrue(cert.exact_degree_valid)
            self.assertTrue(replay_stored_witness(p)['certified'])

    def test_positive_eta_reduces_margin_monotonically(self):
        for p in STORED_WITNESSES:
            a=certify_phase_neighborhood(p,0)
            b=certify_phase_neighborhood(p,Q(1,10**7))
            c=certify_phase_neighborhood(p,Q(1,10**6))
            self.assertEqual(a.margin,b.margin)
            self.assertEqual(b.margin,c.margin)
            self.assertEqual(b.four_E_total-a.four_E_total,Q(4,10**7))
            self.assertEqual(c.four_E_total-a.four_E_total,Q(4,10**6))
            self.assertTrue(c.exact_degree_valid)

    def test_too_large_phase_change_fails_closed(self):
        for p in STORED_WITNESSES:
            c=certify_phase_neighborhood(p,Q(1,100))
            self.assertFalse(c.exact_degree_valid)

    def test_invalid_phase_and_unproved_exponent_are_rejected(self):
        with self.assertRaises(ValueError):certify_phase_neighborhood(Q(2),Q(-1,1000))
        with self.assertRaises(ValueError):certify_phase_neighborhood(Q(7,3))

    def test_replay_is_independent_of_numerical_discovery(self):
        with patch('spiral_multicandidate_v232.scan_candidates',side_effect=AssertionError('search forbidden')):
            for p in STORED_WITNESSES:
                self.assertTrue(certify_phase_neighborhood(p).exact_degree_valid)

    def test_machine_report_has_three_certified_entries(self):
        data=report()['phase_neighborhoods']
        self.assertEqual(len(data),3)
        self.assertTrue(all(e['degree_valid'] and e['strict_margin'] for e in data))


class SearchAccountingTests(unittest.TestCase):
    def test_signature_ignores_dihedral_labeling(self):
        pts=((0.,0.),(1.,0.),(1.,1.),(0.,1.))
        self.assertEqual(_source_signature(pts),_source_signature(tuple(reversed(pts))))
        self.assertEqual(_source_signature(pts),_source_signature(pts[1:]+pts[:1]))

    def test_bad_budgets_fail_before_numerical_search(self):
        for kwargs in ({'max_pairs':0},{'topk':0},{'max_candidates':0},{'min_margin':1.0}):
            with self.assertRaises(ValueError):scan_candidates(2,.55,**kwargs)

    def test_small_budget_reports_truncation(self):
        s=scan_candidates(2,inner=.55,step=.22,max_pairs=3,max_candidates=4)
        self.assertEqual(s.pairs_examined,3)
        self.assertTrue(s.pair_budget_exhausted)
        self.assertLessEqual(len(s.distinct_candidates),s.cap)

    def test_scan_returns_diverse_proposals_with_bounded_memory(self):
        scan=scan_candidates(2,inner=.22,step=.22,max_pairs=2500,max_candidates=8)
        self.assertEqual(scan.samples>10,True)
        self.assertLessEqual(scan.pairs_examined,2500)
        self.assertLessEqual(len(scan.distinct_candidates),8)
        self.assertGreater(len(scan.distinct_candidates),0)
        keys=[_source_signature(c['points']) for c in scan.distinct_candidates]
        self.assertEqual(len(set(keys)),len(keys))
        self.assertTrue(all(c['status']=='found' for c in scan.distinct_candidates))

    def test_empty_candidate_set_is_inconclusive(self):
        r=certify_candidates([])
        self.assertFalse(r['certified'])
        self.assertEqual(r['status'],'candidate-exhausted-inconclusive')
        self.assertEqual(r['attempted'],0)

    def test_bad_certification_budget_rejected(self):
        with self.assertRaises(ValueError):certify_candidates([],max_checks=0)

    def test_rejected_first_candidate_does_not_prevent_second(self):
        first={'status':'found','ids':(0,1,2,3)}
        second={'status':'found','ids':(4,5,6,7)}
        with patch('spiral_multicandidate_v232.local_arcs_at_newton',side_effect=lambda c,half_angle:c):
            with patch('spiral_multicandidate_v232.certify_on_smooth_spiral',
                       side_effect=[{'certified':False,'reason':'uncertified'},
                                    {'certified':True,'witness':'second'}]):
                r=certify_candidates([first,second],half_angles=(.01,),max_checks=2)
        self.assertTrue(r['certified'])
        self.assertEqual(r['candidate_rank'],2)
        self.assertEqual(r['attempted'],2)
        self.assertEqual(r['rejections'],{'uncertified':1})

    def test_exception_in_one_candidate_does_not_stop_search(self):
        with patch('spiral_multicandidate_v232.local_arcs_at_newton',
                   side_effect=[ValueError('bad Newton'), {'status':'found'}]):
            with patch('spiral_multicandidate_v232.certify_on_smooth_spiral',
                       return_value={'certified':True}):
                r=certify_candidates([{},{}],half_angles=(.01,),max_checks=2)
        self.assertTrue(r['certified'])
        self.assertEqual(r['rejections'],{'ValueError':1})

    def test_hard_limit_reports_budget_inconclusive(self):
        with patch('spiral_multicandidate_v232.local_arcs_at_newton',return_value={}):
            with patch('spiral_multicandidate_v232.certify_on_smooth_spiral',
                       return_value={'certified':False,'reason':'insufficient-margin'}):
                r=certify_candidates([{},{}],half_angles=(.01,.02),max_checks=1)
        self.assertEqual(r['status'],'budget-inconclusive')
        self.assertFalse(r['certified'])
        self.assertEqual(r['attempted'],1)

    def test_fully_empty_adaptive_ladder_is_not_false_disproof(self):
        s=CandidateScan(10,2,False,0,(),4)
        with patch('spiral_multicandidate_v232.scan_candidates',return_value=s):
            r=adaptive_multicandidate(2,radius_ladder=(.55,.45),max_pairs=2)
        self.assertFalse(r['certified'])
        self.assertEqual(r['status'],'search-inconclusive')
        self.assertEqual(len(r['search_ladder']),2)

    def test_single_radius_multi_proposal_can_succeed_at_second_candidate(self):
        def cand(i):return {'status':'found','ids':(i,i+1,i+2,i+3)}
        s=CandidateScan(100,500,True,12,(cand(0),cand(4)),2)
        with patch('spiral_multicandidate_v232.scan_candidates',return_value=s):
            with patch('spiral_multicandidate_v232.local_arcs_at_newton',side_effect=lambda c,half_angle:c):
                with patch('spiral_multicandidate_v232.certify_on_smooth_spiral',
                           side_effect=[{'certified':False,'reason':'margin'},
                                        {'certified':True,'witness':'later'}]):
                    r=adaptive_multicandidate(2,radius_ladder=(.55,),
                          max_checks_per_radius=2,half_angles=(.01,))
        self.assertTrue(r['certified'])
        self.assertEqual(r['candidate_rank'],2)
        self.assertEqual(r['rejections'],{'margin':1})


if __name__=='__main__':unittest.main()
