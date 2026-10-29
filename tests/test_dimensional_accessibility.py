"""Fast regression tests; full numerical suites are executable under code/."""
import sys
import unittest
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
import reproduce
import verify_extensions


class DimensionalAccessibilityTests(unittest.TestCase):
    def test_exact_kernel(self):
        for s in (0.01, 1., 100.):
            self.assertAlmostEqual(reproduce.kernel(1.5, s)/(1/s**2+1/s), 1., places=12)

    def test_aqual_eigenvalues(self):
        a = verify_extensions.aqual_tensor(np.array([0., 0., 1.]))
        np.testing.assert_allclose(np.linalg.eigvalsh(a),
                                   [1/np.sqrt(2), 1/np.sqrt(2), 3/(2*np.sqrt(2))])

    def test_invalid_arguments(self):
        with self.assertRaises(ValueError):
            reproduce.kernel(1., 1.)
        with self.assertRaises(ValueError):
            reproduce.detuning_root(2., 0.01)
        with self.assertRaises(ValueError):
            verify_extensions.aqual_tensor(np.zeros(3))


if __name__ == '__main__':
    unittest.main()
