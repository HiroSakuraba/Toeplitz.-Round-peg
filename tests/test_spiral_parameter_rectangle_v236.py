"""Check mixed derivatives and failure controls under SageMath."""

import importlib.util
from pathlib import Path
import sys
import unittest

if importlib.util.find_spec("sage") is None:
    raise unittest.SkipTest("The v236 interval tests run in the SageMath workflow")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "code"))
from sage.all import QQ, SR, cos, sin, var
from spiral_parameter_rectangle_v236 import (
    DOMAIN, FIELD, MANIFEST, _ray_acceleration, _split, load_tree, replay, verify_tile,
)


class ParameterRectangleTests(unittest.TestCase):
    def test_mixed_ray_acceleration_encloses_symbolic_second_derivative(self):
        t = var("t")
        u0, lam0, amp0 = QQ(4)/7, QQ(2)/3, -QQ(1)/250
        u_speed, lam_speed, amp_speed = QQ(3)/170, -QQ(1)/23, QQ(1)/41
        slope, intercept = QQ(7)/5, -QQ(2)/7
        u, lam, amp = u0+u_speed*t, lam0+lam_speed*t, amp0+amp_speed*t
        phase = u**(-4) + lam*u**(-2) + amp*(slope*u+intercept)
        result = _ray_acceleration(FIELD(u0), FIELD(lam0), FIELD(amp0),
                                   FIELD(u_speed), (FIELD(lam_speed), FIELD(amp_speed)),
                                   (slope, intercept))
        for actual, trig in zip(result, (cos, sin)):
            expr = SR(u*u*trig(phase)).diff(t, 2).subs(t=0)
            self.assertIn(expr.n(prec=160), actual)

    def test_cover_rejects_missing_child(self):
        with self.assertRaisesRegex(ValueError, "invalid certificate"):
            replay([0, [], []])

    def test_first_leaf_and_orientation(self):
        bounds = DOMAIN
        node = load_tree(MANIFEST)
        while len(node) == 3:
            axis, node, _ = node
            bounds, _ = _split(bounds, axis)
        orientation, *box = node
        self.assertEqual(verify_tile(bounds, orientation, *box), "certified")
        self.assertNotEqual(verify_tile(bounds, -orientation, *box), "certified")
