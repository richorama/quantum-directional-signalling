import unittest
from fractions import Fraction

from quantum_coarse_graining.cartan import (
    cartan_effective_channel,
    cartan_effective_channel_is_cptp,
    cartan_joint_pauli_covariance_certificate,
    cartan_pauli_coefficients,
    cartan_pauli_weights,
    cartan_symmetry_certificate,
    cartan_unitary,
    cartan_visible_sector_certificate,
)
from quantum_coarse_graining.exact import Gaussian, is_unitary, matrix_unit
from quantum_coarse_graining.ising import effective_dephasing, ising_unitary

IDENTITY = (Fraction(1), Fraction(0))
THREE_FOUR = (Fraction(3, 5), Fraction(4, 5))
FOUR_THREE = (Fraction(4, 5), Fraction(3, 5))


class CartanSymmetryTests(unittest.TestCase):
    def test_exact_cartan_unitaries(self):
        for parameters in (
            (IDENTITY, IDENTITY, THREE_FOUR),
            (THREE_FOUR, FOUR_THREE, IDENTITY),
            (THREE_FOUR, FOUR_THREE, THREE_FOUR),
        ):
            self.assertTrue(is_unitary(cartan_unitary(parameters)))

    def test_single_z_factor_reduces_to_ising(self):
        parameters = (IDENTITY, IDENTITY, THREE_FOUR)
        self.assertEqual(
            cartan_unitary(parameters),
            ising_unitary(*THREE_FOUR),
        )
        coefficients = cartan_pauli_coefficients(parameters)
        self.assertEqual(
            coefficients,
            (
                Gaussian(Fraction(3, 5)),
                Gaussian(),
                Gaussian(),
                Gaussian(0, Fraction(-4, 5)),
            ),
        )
        operator = matrix_unit(2, 0, 1)
        self.assertEqual(
            cartan_effective_channel(operator, parameters),
            effective_dephasing(operator, *THREE_FOUR),
        )

    def test_pauli_weights_certify_a_cptp_channel(self):
        parameters = (THREE_FOUR, FOUR_THREE, THREE_FOUR)
        weights = cartan_pauli_weights(parameters)
        self.assertEqual(sum(weights), 1)
        self.assertTrue(all(weight >= 0 for weight in weights))
        self.assertTrue(cartan_effective_channel_is_cptp(parameters))

    def test_candidate_is_exact_on_the_visible_sector(self):
        parameters = (THREE_FOUR, FOUR_THREE, THREE_FOUR)
        self.assertTrue(cartan_visible_sector_certificate(parameters))

    def test_joint_pauli_covariance(self):
        parameters = (THREE_FOUR, FOUR_THREE, THREE_FOUR)
        self.assertTrue(cartan_joint_pauli_covariance_certificate(parameters))
        self.assertTrue(cartan_symmetry_certificate(parameters))

    def test_rejects_invalid_factor(self):
        with self.assertRaises(ValueError):
            cartan_unitary(
                (
                    IDENTITY,
                    (Fraction(1, 2), Fraction(1, 2)),
                    IDENTITY,
                )
            )


if __name__ == "__main__":
    unittest.main()
