import importlib.util
import math
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
SCRIPT = ROOT / "paper" / "generate_xy_plot.py"


def load_plot_module():
    specification = importlib.util.spec_from_file_location("generate_xy_plot", SCRIPT)
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


class XYPlotTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.plot = load_plot_module()

    def test_curve_endpoints(self):
        self.assertEqual(self.plot.signalling_defect(0), (0.0, "weak"))
        defect, regime = self.plot.signalling_defect(math.pi / 4)
        self.assertEqual(regime, "strong")
        self.assertAlmostEqual(defect, 1.5)

    def test_middle_root_satisfies_quartic(self):
        theta = 0.5
        defect, regime = self.plot.signalling_defect(theta)
        coefficients = self.plot.quartic_coefficients(math.tan(theta))
        self.assertEqual(regime, "middle")
        self.assertAlmostEqual(
            self.plot.evaluate_polynomial(coefficients, defect),
            0,
            places=10,
        )

    def test_generated_points_cover_all_regimes(self):
        regimes = {regime for _, _, regime in self.plot.points()}
        self.assertEqual(regimes, {"weak", "middle", "strong"})

    def test_partial_swap_optimizer_endpoints(self):
        self.assertEqual(self.plot.partial_swap_optimizer(0), (1.0, 0.5))
        shrinkage, symmetric_weight = self.plot.partial_swap_optimizer(
            math.pi / 2
        )
        self.assertEqual(shrinkage, 0.0)
        self.assertAlmostEqual(symmetric_weight, 0.75)

    def test_partial_swap_interior_optimizer_matches_scalar_value(self):
        phi = math.asin(0.8)
        shrinkage, _ = self.plot.partial_swap_optimizer(phi)
        self.assertAlmostEqual(shrinkage, 0.339610, places=5)
        self.assertAlmostEqual(
            self.plot.partial_swap_fixed_defect(0.8, shrinkage),
            1.415252,
            places=5,
        )

    def test_xy_optimizer_tracks_witness_transition(self):
        self.assertEqual(
            self.plot.xy_optimizer(0),
            (1.0, 1.0, 0.0, "weak"),
        )
        transverse, longitudinal, even_weight, regime = (
            self.plot.xy_optimizer(0.45)
        )
        self.assertEqual(regime, "middle")
        self.assertAlmostEqual(longitudinal, 2 * transverse - 1)
        self.assertGreater(even_weight, 0)
        self.assertLess(even_weight, 0.5)
        defect, _ = self.plot.signalling_defect(0.45)
        self.assertAlmostEqual(
            4
            * self.plot.xy_largest_eigenvalue(
                0.45,
                transverse,
                even_weight,
            ),
            defect,
            places=6,
        )
        transverse, longitudinal, even_weight, regime = (
            self.plot.xy_optimizer(math.pi / 4)
        )
        self.assertAlmostEqual(transverse, 0)
        self.assertAlmostEqual(longitudinal, 0)
        self.assertEqual(even_weight, 0.25)
        self.assertEqual(regime, "strong")


if __name__ == "__main__":
    unittest.main()
