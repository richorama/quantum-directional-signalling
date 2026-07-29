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


if __name__ == "__main__":
    unittest.main()
