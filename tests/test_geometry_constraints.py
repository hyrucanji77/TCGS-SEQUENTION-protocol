"""Fast regression checks for the geometric and power-law-conversion supplement."""
import math
import sys
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'code'))
import verify_geometry_constraints as geometry


class GeometryConstraintTests(unittest.TestCase):
    def test_hyperbolic_and_flat_sections(self):
        self.assertEqual(geometry.sectional_curvatures(0., 2., 3.), (0., 0.))
        self.assertEqual(geometry.sectional_curvatures(1., 2., 3.), (-.25, -.25))

    def test_published_conversion_confidence(self):
        limits = geometry.published_limits()['sources']
        self.assertEqual([r['published_confidence_percent'] for r in limits], [68, 68, 68])
        self.assertAlmostEqual(limits[1]['ell_limit_micrometres'], .7068583470577035)

    def test_hoskins_historical_ordering(self):
        rows = {r['id']: r for r in geometry.published_limits()['sources']}
        self.assertLess(rows['hoskins1985']['beta2_abs_limit'],
                        rows['hoyle2004']['beta2_abs_limit'])
        self.assertLess(rows['adelberger2007']['beta2_abs_limit'],
                        rows['hoskins1985']['beta2_abs_limit'])
        self.assertAlmostEqual(rows['hoskins1985']['ell_limit_micrometres'],
                               2.0420352248333655)

    def test_published_provenance(self):
        rows = {r['id']: r for r in geometry.published_limits()['sources']}
        p = rows['hoskins1985']['coefficient_and_confidence_provenance']
        self.assertEqual(p['doi'], '10.1103/PhysRevLett.98.131104')
        self.assertEqual(p['pdf_page_1based'], 3)
        self.assertIn('reference 15', p['table'])
        for row in geometry.published_limits()['later_literature']:
            self.assertFalse(row['quantitative_beta2_input_used'])

    def test_conditional_scale_hierarchy(self):
        data = geometry.cross_sector_hierarchy(geometry.convert_beta2(4.5e-4), 1.2e-10)
        self.assertTrue(1.27e23 < data['comparison_acceleration_m_s2'] < 1.28e23)
        self.assertTrue(9.4e-34 < data['xi_times_eta_H'] < 9.5e-34)
        self.assertAlmostEqual(data['scale_ratio']*data['xi_times_eta_H'], 1.)
        self.assertAlmostEqual(data['unsuppressed_curvature_radius_m']/data['ell_m'],
                               data['scale_ratio'], delta=1e-14*data['scale_ratio'])
        for ell, acceleration in ((0., 1.), (1., -1.), (math.inf, 1.)):
            with self.assertRaises(ValueError):
                geometry.cross_sector_hierarchy(ell, acceleration)

    def test_invalid_parameters(self):
        for beta, ref in ((-1., 1.), (1., 0.), (math.inf, 1.)):
            with self.assertRaises(ValueError):
                geometry.convert_beta2(beta, ref)
        with self.assertRaises(ValueError):
            geometry.sectional_curvatures(1., 0., 1.)


if __name__ == '__main__':
    unittest.main()
