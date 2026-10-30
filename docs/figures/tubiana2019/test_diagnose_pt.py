"""Sanity checks for the convergence and conditional-frequency calculations."""
import unittest

import arviz as az
import numpy as np

from diagnose_pt import diagnostic, sequence_covariance, project_sequences


class DiagnosticsTests(unittest.TestCase):
    def test_memory_bounded_sequence_statistics(self):
        rng = np.random.default_rng(91)
        sequences = rng.integers(0, 20, size=(2, 103, 27), dtype=np.int8)
        x = np.eye(20)[sequences].reshape(2, 103, 540)
        center, covariance = sequence_covariance(sequences.reshape(-1, 27))
        expected_center = x.mean(axis=(0, 1))
        centered = (x - expected_center).reshape(-1, 540)
        np.testing.assert_allclose(center, expected_center, atol=1e-15)
        np.testing.assert_allclose(covariance, centered.T @ centered / len(centered), atol=1e-14)
        vectors = rng.normal(size=(540, 5))
        np.testing.assert_allclose(project_sequences(sequences, center, vectors),
                                   (x-center) @ vectors, atol=1e-13)

    def test_independent_draws(self):
        values = np.random.default_rng(41).normal(size=(4, 4000))
        result = diagnostic(values)
        self.assertLess(result["rhat"], 1.01)
        self.assertGreater(result["bulk_ess"], 8000)

    def test_separated_chains(self):
        values = np.random.default_rng(41).normal(size=(4, 4000))
        values += 3 * np.arange(4)[:, None]
        result = diagnostic(values)
        self.assertGreater(result["rhat"], 1.5)
        self.assertLess(result["bulk_ess"], 20)

    def test_constant_is_not_evidence_of_mixing(self):
        result = diagnostic(np.ones((4, 1000)))
        self.assertIsNone(result["rhat"])
        self.assertIsNone(result["bulk_ess"])

    def test_chains_frozen_at_different_values(self):
        values = np.repeat(np.arange(4)[:, None], 1000, axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            result = diagnostic(values)
        self.assertTrue(np.isinf(result["rhat"]))
        self.assertLess(result["bulk_ess"], 10)

    def test_conditional_frequency_mcse(self):
        rng = np.random.default_rng(12)
        keep = rng.random((4, 10000)) < .6
        residue = rng.random(keep.shape) < .3
        conditional = residue[keep].mean()
        influence = keep * (residue - conditional) / keep.mean()
        actual = float(az.mcse(influence, method="mean"))
        expected = np.sqrt(.3 * .7 / (.6 * keep.size))
        self.assertAlmostEqual(actual / expected, 1, delta=.1)


if __name__ == "__main__":
    unittest.main()
