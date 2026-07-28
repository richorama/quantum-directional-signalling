import unittest
from fractions import Fraction

from quantum_coarse_graining.exact import Gaussian, is_unitary, matrix_unit
from quantum_coarse_graining.xy import (
    xy_block_characteristic,
    xy_boundary_certificate,
    xy_first_threshold_tangent_polynomial,
    xy_middle_quartic_coefficients,
    xy_second_threshold_tangent_polynomial,
    xy_pauli_channel,
    xy_strong_defect,
    xy_strong_probabilities,
    xy_strong_spectrum_certificate,
    xy_strong_threshold_polynomial,
    xy_transfer_parameters,
    xy_unitary,
    xy_weak_defect,
    xy_weak_probabilities,
    xy_weak_threshold_polynomial,
)


class XYSignallingTests(unittest.TestCase):
    def test_exact_xy_unitary_and_transfer_angle(self):
        factor = (Fraction(4, 5), Fraction(3, 5))
        self.assertTrue(is_unitary(xy_unitary(factor)))
        self.assertEqual(
            xy_transfer_parameters(factor),
            (Fraction(7, 25), Fraction(24, 25)),
        )

    def test_weak_branch_exact_control(self):
        cosine, sine = Fraction(12, 13), Fraction(5, 13)
        self.assertLess(xy_weak_threshold_polynomial(cosine, sine), 0)
        self.assertEqual(xy_weak_defect(cosine, sine), Fraction(120, 169))
        self.assertEqual(sum(xy_weak_probabilities(cosine, sine)), 1)

    def test_strong_branch_exact_control(self):
        cosine, sine = Fraction(9, 41), Fraction(40, 41)
        self.assertGreater(xy_strong_threshold_polynomial(cosine, sine), 0)
        self.assertEqual(xy_strong_defect(cosine, sine), Fraction(2440, 1681))
        self.assertEqual(sum(xy_strong_probabilities(cosine, sine)), 1)
        self.assertTrue(xy_strong_spectrum_certificate(cosine, sine))

    def test_iswap_endpoint(self):
        self.assertEqual(
            xy_strong_probabilities(Fraction(0), Fraction(1)),
            (Fraction(1, 4),) * 4,
        )
        self.assertEqual(
            xy_strong_defect(Fraction(0), Fraction(1)),
            Fraction(3, 2),
        )

    def test_pauli_channel_preserves_trace(self):
        probabilities = xy_strong_probabilities(
            Fraction(9, 41),
            Fraction(40, 41),
        )
        image = xy_pauli_channel(matrix_unit(2, 0, 0), probabilities)
        self.assertEqual(image[0][0] + image[1][1], Gaussian(1))

    def test_block_characteristic_is_monic(self):
        coefficients = xy_block_characteristic(
            Fraction(9, 41),
            Fraction(40, 41),
            Fraction(729, 1681),
            Fraction(81, 1681),
            Fraction(1, 4),
        )
        self.assertEqual(coefficients[0], 1)

    def test_branch_guards(self):
        with self.assertRaises(ValueError):
            xy_weak_defect(Fraction(3, 5), Fraction(4, 5))
        with self.assertRaises(ValueError):
            xy_strong_defect(Fraction(3, 5), Fraction(4, 5))
        with self.assertRaises(ValueError):
            xy_unitary((Fraction(1, 2), Fraction(1, 2)))

    def test_combined_boundary_certificate(self):
        self.assertTrue(xy_boundary_certificate())

    def test_conditional_middle_quartic_control(self):
        self.assertEqual(
            xy_middle_quartic_coefficients(Fraction(1, 2)),
            (
                Fraction(175, 64),
                Fraction(-3321, 64),
                Fraction(-977, 32),
                Fraction(411, 4),
                Fraction(5),
            ),
        )
        self.assertEqual(
            xy_first_threshold_tangent_polynomial(Fraction(1, 2)),
            Fraction(97, 64),
        )
        self.assertEqual(
            xy_second_threshold_tangent_polynomial(Fraction(1, 2)),
            Fraction(-5, 8),
        )


if __name__ == "__main__":
    unittest.main()
