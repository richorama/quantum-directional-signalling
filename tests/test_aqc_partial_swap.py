import unittest
from fractions import Fraction

from quantum_coarse_graining.exact import is_unitary, matrix_unit
from quantum_coarse_graining.partial_swap import (
    depolarizing_channel,
    equal_cartan_is_partial_swap,
    fixed_channel_defect,
    fixed_channel_terms,
    partial_swap_boundary_certificate,
    partial_swap_unitary,
    qudit_fixed_channel_defect,
    qudit_fixed_channel_terms,
    qudit_partial_swap_certificate,
    qudit_weak_partial_swap_defect,
    qudit_weak_threshold,
    strong_stationarity_holds,
    weak_partial_swap_defect,
)


class PartialSwapTests(unittest.TestCase):
    def test_partial_swap_is_unitary(self):
        self.assertTrue(
            is_unitary(
                partial_swap_unitary(Fraction(40, 41), Fraction(9, 41))
            )
        )

    def test_equal_cartan_angles_reduce_to_partial_swap(self):
        self.assertTrue(
            equal_cartan_is_partial_swap((Fraction(4, 5), Fraction(3, 5)))
        )

    def test_depolarizing_channel_range(self):
        operator = matrix_unit(2, 0, 1)
        self.assertEqual(depolarizing_channel(operator, 1), operator)
        with self.assertRaises(ValueError):
            depolarizing_channel(operator, Fraction(-1, 2))

    def test_weak_branch_is_exact(self):
        cosine = Fraction(40, 41)
        sine = Fraction(9, 41)
        self.assertEqual(
            fixed_channel_terms(cosine, sine, 1),
            (Fraction(0), 4 * sine * sine),
        )
        self.assertEqual(
            weak_partial_swap_defect(cosine, sine),
            Fraction(18, 41),
        )
        self.assertEqual(
            fixed_channel_defect(cosine, sine, 1, 2 * sine),
            Fraction(18, 41),
        )

    def test_weak_branch_rejects_strong_coupling(self):
        with self.assertRaises(ValueError):
            weak_partial_swap_defect(Fraction(4, 5), Fraction(3, 5))

    def test_swap_endpoint_recovers_three_halves(self):
        self.assertEqual(
            fixed_channel_defect(0, 1, 0, 1),
            Fraction(3, 2),
        )
        self.assertTrue(strong_stationarity_holds(0, 1, 0, 1))
        self.assertTrue(partial_swap_boundary_certificate())

    def test_strong_stationarity_rejects_the_first_norm_branch(self):
        self.assertFalse(
            strong_stationarity_holds(
                Fraction(0),
                Fraction(1),
                Fraction(-1, 3),
                Fraction(2, 3),
            )
        )

    def test_rejects_incorrect_radical(self):
        with self.assertRaises(ValueError):
            fixed_channel_defect(
                Fraction(40, 41),
                Fraction(9, 41),
                1,
                Fraction(1),
            )

    def test_qutrit_weak_threshold_and_defect(self):
        self.assertEqual(qudit_weak_threshold(3), Fraction(3, 4))
        self.assertEqual(
            qudit_weak_partial_swap_defect(
                3,
                Fraction(4, 5),
                Fraction(3, 5),
            ),
            Fraction(6, 5),
        )
        self.assertEqual(
            qudit_fixed_channel_defect(
                3,
                Fraction(4, 5),
                Fraction(3, 5),
                1,
                Fraction(18, 5),
            ),
            Fraction(6, 5),
        )

    def test_qutrit_swap_endpoint(self):
        terms = qudit_fixed_channel_terms(3, 0, 1, 0)
        self.assertEqual(terms[4], 9)
        self.assertEqual(
            qudit_fixed_channel_defect(3, 0, 1, 0, 3),
            Fraction(16, 9),
        )
        self.assertTrue(qudit_partial_swap_certificate())

    def test_qutrit_negative_cptp_endpoint_uses_boundary_maximum(self):
        self.assertEqual(
            qudit_fixed_channel_defect(
                3,
                0,
                1,
                Fraction(-1, 8),
                Fraction(21, 8),
            ),
            Fraction(15, 8),
        )

    def test_dimension_four_exact_controls(self):
        self.assertEqual(qudit_weak_threshold(4), Fraction(13, 15))
        self.assertEqual(
            qudit_fixed_channel_defect(
                4,
                Fraction(3, 5),
                Fraction(4, 5),
                1,
                Fraction(32, 5),
            ),
            Fraction(8, 5),
        )
        self.assertEqual(
            qudit_fixed_channel_defect(4, 0, 1, 0, 4),
            Fraction(15, 8),
        )

    def test_qudit_rejects_invalid_dimension(self):
        with self.assertRaises(ValueError):
            qudit_weak_threshold(1)


if __name__ == "__main__":
    unittest.main()
