"""v231 exact a-posteriori degree certificates for three irregular spiral curves."""
import sys
import unittest
from pathlib import Path
from fractions import Fraction as Q
from math import isqrt,pi,sin,cos,hypot
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))

from spiral_pruned_v231 import (
    _atan_interval,_pi_interval,_trig_rational_enclosure,_round_fraction,
    _R_residual,_R_matrix,_inverse_and_root,_certified_point,_chord_error,
    spiral_samples,SpatialEdges,Sample,solve_chart,discover,newton_square,
    local_arcs_at_newton,certify_on_smooth_spiral,certify_discovered_square,
    blind_discover_and_certify,replay_stored_witness,STORED_WITNESSES,
    smooth_spiral,
)

class TrigEnclosures(unittest.TestCase):
    def test_pi_interval_machin_positivity(self):
        lo,hi=_pi_interval()
        self.assertLess(lo,hi)
        self.assertLess(Q('3.141592653589793'),lo)
        self.assertLess(hi,Q('3.141592653589794'))

    def test_rational_atan_interval(self):
        for x in (Q(1,5),Q(1,239)):
            lo,hi=_atan_interval(x)
            self.assertLess(lo,hi)
            self.assertLess(hi-lo,Q(1,10**25))

    def test_trigonometric_enclosures_independent_float_sanity(self):
        # Float sanity checks are not used as evidence of interval correctness.
        for t in (Q(1),Q(3,2),Q(7),Q(37,3),Q(145),Q(299)):
            s,c,e=_trig_rational_enclosure(t)
            self.assertGreater(e,0)
            self.assertLess(abs(float(s)-sin(float(t))),float(e)+2e-14)
            self.assertLess(abs(float(c)-cos(float(t))),float(e)+2e-14)

    def test_rounding_rational(self):
        self.assertEqual(_round_fraction(Q(3,10),10),Q(3,10))
        self.assertEqual(_round_fraction(-Q(3,10),10),-Q(3,10))
        self.assertEqual(_round_fraction(Q(17,100),10),Q(1,5))

    def test_endpoint_enclosure(self):
        for p in (Q(2),Q(5,2),Q(3)):
            _,us=STORED_WITNESSES[p]
            for ui in us:
                u=Q(ui,2**30);t=u*u
                v,e=_certified_point(Sample(t,0,0),p)
                self.assertTrue(0<t<=1)
                self.assertGreater(e,0)
                th=float(u**(-int(2*p)))
                self.assertLess(abs(float(v[0])-float(t)*cos(th)),float(e)+2e-14)
                self.assertLess(abs(float(v[1])-float(t)*sin(th)),float(e)+2e-14)

class DegreeLemma(unittest.TestCase):
    def test_identity_matrix_inverse(self):
        A=tuple(tuple(Q(i==j) for j in range(4)) for i in range(4))
        u,M=_inverse_and_root(A,tuple(Q(1,2) for _ in range(4)))
        self.assertEqual(u,(Q(1,2),)*4)
        self.assertEqual(M,1)

    def test_affine_four_vertex_residual(self):
        from spiral_pruned_v231 import _R_matrix
        pts=((Q(0),Q(0)),(Q(1),Q(0)),(Q(1),Q(1)),(Q(0),Q(1)))
        self.assertEqual(_R_residual(pts,1),(0,0,0,0))
        A,b=_R_matrix(pts,((0,1),)*4,1)
        self.assertEqual(len(A),4)
        self.assertEqual(len(b),4)

    def test_chord_bound_exact_and_positive(self):
        for p in (Q(2),Q(5,2),Q(3)):
            lo,hi=Q(2,5),Q(21,50)
            self.assertGreater(_chord_error(lo,hi,p),0)
            self.assertEqual(_chord_error(lo,lo,p),0)
            self.assertEqual(_chord_error(lo,hi,p),_chord_error(hi,lo,p))

    def test_chord_second_derivative_bound_numerical_sanity(self):
        for p in (2,2.5,3):
            lo,hi=Q(2,5),Q(21,50)
            B=float(_chord_error(lo,hi,p))
            a=smooth_spiral(float(lo),p);b=smooth_spiral(float(hi),p)
            for j in range(21):
                u=j/20
                z=smooth_spiral(float(lo)*(1-u)+float(hi)*u,p)
                chord=(a[0]*(1-u)+b[0]*u,a[1]*(1-u)+b[1]*u)
                self.assertLess(hypot(z[0]-chord[0],z[1]-chord[1]),B+1e-13)

    def test_interval_errors_from_true_spiral_at_stored_endpoints(self):
        for p,(_,us) in STORED_WITNESSES.items():
            for ui in us:
                t=Q(ui,2**30)**2
                pt,e=_certified_point(Sample(t,0,0),p)
                self.assertLess(float(e),1e-8)
                z=smooth_spiral(float(t),float(p))
                self.assertLess(abs(float(pt[0])-z[0]),float(e)+5e-14)
                self.assertLess(abs(float(pt[1])-z[1]),float(e)+5e-14)

class CertificateReplay(unittest.TestCase):
    def test_machine_readable_manifest_matches_proof_code(self):
        import json
        manifest=Path(__file__).resolve().parents[1]/'research'/'v231-certificate-manifest.json'
        data=json.loads(manifest.read_text())
        from spiral_pruned_v231 import STORED_STRONG_FLOORS
        self.assertEqual(len(data['witnesses']),3)
        for entry in data['witnesses']:
            p=Q(entry['p'])
            sig,us=STORED_WITNESSES[p]
            self.assertEqual(sig,entry['sigma'])
            self.assertEqual(tuple(us),tuple(entry['u']))
            self.assertEqual(STORED_STRONG_FLOORS[p],
                             (Q(entry['side_strict_floor']),Q(entry['radius_strict_floor'])))
            self.assertTrue(replay_stored_witness(p)['exact_separator_4E_lt_1over20000_lt_boundary'])

    def test_stored_sources_have_distinct_parameter_intervals(self):
        for p,(_,us) in STORED_WITNESSES.items():
            blocks=[(Q(us[i],2**30)**2,Q(us[i+1],2**30)**2) for i in (0,2,4,6)]
            for lo,hi in blocks:self.assertLess(lo,hi)
            ordered=sorted(blocks)
            for x,y in zip(ordered,ordered[1:]):self.assertLess(x[1],y[0])

    def test_replay_p2(self):
        x=replay_stored_witness(2)
        self.assertTrue(x['certified'],x)
        self.assertLess(x['max_residual_perturbation'],x['degree_boundary_margin'])
        self.assertEqual(x['certified_side_strictly_greater_than'],'1/4')
        self.assertEqual(x['certified_min_radius_strictly_greater_than'],'6/25')
        self.assertTrue(x['exact_separator_4E_lt_1over20000_lt_boundary'])

    def test_replay_p5_2(self):
        x=replay_stored_witness(Q(5,2))
        self.assertTrue(x['certified'],x)
        self.assertEqual(x['certified_side_strictly_greater_than'],'3/5')
        self.assertEqual(x['certified_min_radius_strictly_greater_than'],'3/10')
        self.assertTrue(x['exact_separator_4E_lt_1over20000_lt_boundary'])
        self.assertLess(x['max_residual_perturbation'],x['degree_boundary_margin'])

    def test_replay_p3(self):
        x=replay_stored_witness(3)
        self.assertTrue(x['certified'],x)
        self.assertEqual(x['certified_side_strictly_greater_than'],'7/10')
        self.assertEqual(x['certified_min_radius_strictly_greater_than'],'3/8')
        self.assertTrue(x['exact_separator_4E_lt_1over20000_lt_boundary'])
        self.assertLess(x['max_residual_perturbation'],x['degree_boundary_margin'])

    def test_sign_reversal_is_not_a_certified_square(self):
        for p,(sig,us) in STORED_WITNESSES.items():
            samps=tuple(Sample(Q(ui,2**30)**2,0.,0.) for ui in us)
            x={'status':'found','p':p,'inner':None,'samples_data':samps,'ids':(0,2,4,6),'sigma':-sig}
            r=certify_on_smooth_spiral(x)
            self.assertFalse(r['certified'],(p,r))

    def test_impossible_side_floor_refused(self):
        p=Q(2);sig,us=STORED_WITNESSES[p]
        samps=tuple(Sample(Q(ui,2**30)**2,0.,0.) for ui in us)
        x={'status':'found','p':p,'inner':None,'samples_data':samps,'ids':(0,2,4,6),'sigma':sig}
        r=certify_on_smooth_spiral(x,side_floor=Q(2))
        self.assertFalse(r['certified'])
        self.assertEqual(r['reason'],'side-floor-insufficient')

    def test_powers_beyond_supported_half_integer_rejected(self):
        with self.assertRaises(ValueError):spiral_samples(Q(7,3),.3)
        with self.assertRaises(ValueError):spiral_samples(0,.3)

class DiscoveryTests(unittest.TestCase):
    def test_grid_nearest_from_segment_boxes(self):
        samps=(Sample(Q(1,4),0,0),Sample(Q(1,2),1,0),Sample(Q(3,4),2,0))
        s=SpatialEdges(samps,cell=.4)
        self.assertIn(0,[i for _,i in s.nearest_edges((.5,0.),max_distance=.1)])
        self.assertEqual(s.nearest_edges((400.,400.)),())

    def test_blind_discovery_p2(self):
        r=blind_discover_and_certify(2)
        self.assertTrue(r['certified'],r)
        self.assertLess(r['search_pairs'],30000)

    def test_blind_discovery_p5_2(self):
        r=blind_discover_and_certify(Q(5,2))
        self.assertTrue(r['certified'],r)
        self.assertLess(r['search_pairs'],30000)

    def test_blind_discovery_p3(self):
        r=blind_discover_and_certify(3)
        self.assertTrue(r['certified'],r)
        self.assertLess(r['search_pairs'],30000)

    def test_no_preloaded_candidate_location(self):
        p=2; candidate=discover(p,.22,.22,max_pairs=10000,topk=4,min_margin=.14)
        self.assertEqual(candidate['status'],'found')
        self.assertEqual(len(set(candidate['ids'])),4)
        t,res,it=newton_square(candidate)
        self.assertLess(res,5e-14)
        self.assertEqual(len(t),4)

if __name__=='__main__':unittest.main()
