#!/usr/bin/env python3
"""
Unit Tests for TCGS-SEQUENTION Protocol
=======================================

Run with: python -m pytest tests/test_protocol.py
Or simply: python tests/test_protocol.py
"""

import sys
import os
import numpy as np

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from protocol import ThreeGateProtocol, ProtocolConfig
from statistics import compute_synchrony_statistic, permutation_test_synchrony


def test_synchrony_statistic_perfect():
    """Test S* = 1.0 for perfectly synchronized times."""
    times = np.array([1000, 1000, 1000, 1000, 1000])
    delta_t = 100
    
    s_stat, ref = compute_synchrony_statistic(times, delta_t)
    
    assert s_stat == 1.0, f"Expected S* = 1.0 for perfect synchrony, got {s_stat}"
    assert ref == 1000, f"Expected reference = 1000, got {ref}"
    print("✓ test_synchrony_statistic_perfect PASSED")


def test_synchrony_statistic_partial():
    """Test S* for partially synchronized times."""
    times = np.array([1000, 1050, 1100, 5000, 6000])  # 3/5 within window
    delta_t = 200
    
    s_stat, ref = compute_synchrony_statistic(times, delta_t)
    
    # Median is 1100, so 1000, 1050, 1100 are within ±200
    assert 0.5 <= s_stat <= 0.7, f"Expected S* around 0.6, got {s_stat}"
    print("✓ test_synchrony_statistic_partial PASSED")


def test_synchrony_statistic_none():
    """Test S* for completely asynchronous times."""
    times = np.array([1000, 5000, 10000, 20000, 50000])
    delta_t = 100
    
    s_stat, ref = compute_synchrony_statistic(times, delta_t)
    
    # Only 1 time (median) should be within window
    assert s_stat <= 0.4, f"Expected low S* for async times, got {s_stat}"
    print("✓ test_synchrony_statistic_none PASSED")


def test_permutation_pvalue_perfect():
    """Test p-value is low for perfect synchrony."""
    times = np.array([1000, 1001, 1002, 1003, 1004])
    delta_t = 100
    
    p_value, null = permutation_test_synchrony(times, delta_t, n_permutations=1000, seed=42)
    
    assert p_value < 0.05, f"Expected p < 0.05 for near-perfect synchrony, got {p_value}"
    print("✓ test_permutation_pvalue_perfect PASSED")


def test_protocol_config():
    """Test configuration hash generation."""
    config1 = ProtocolConfig(delta_t=2000, alpha=0.05)
    config2 = ProtocolConfig(delta_t=2000, alpha=0.05)
    config3 = ProtocolConfig(delta_t=1000, alpha=0.05)
    
    # Same parameters should give same hash
    assert config1.config_hash == config2.config_hash, "Same config should give same hash"
    
    # Different parameters should give different hash
    assert config1.config_hash != config3.config_hash, "Different config should give different hash"
    print("✓ test_protocol_config PASSED")


def test_protocol_g1_analysis():
    """Test G1 gate analysis."""
    config = ProtocolConfig(delta_t=500, s_critical=0.6, alpha=0.05, n_permutations=1000)
    protocol = ThreeGateProtocol(config)
    
    # Synchronized gene
    sync_data = {'pop_A': 1000, 'pop_B': 1100, 'pop_C': 1050, 'pop_D': 1150, 'pop_E': 950}
    result = protocol.analyze_gene_g1('test_gene', sync_data)
    
    assert result.s_statistic > 0.8, f"Expected high S* for synchronized data, got {result.s_statistic}"
    assert result.passes_g1, "Expected G1 pass for synchronized data"
    print("✓ test_protocol_g1_analysis PASSED")


def test_protocol_g2_laboratory():
    """Test G2 gate for laboratory-isolated populations."""
    config = ProtocolConfig()
    protocol = ThreeGateProtocol(config)
    
    populations = ['pop_A', 'pop_B', 'pop_C', 'pop_D']
    result = protocol.analyze_independence_g2(populations, is_laboratory_isolated=True)
    
    assert result.passes_g2, "Laboratory isolation should auto-pass G2"
    print("✓ test_protocol_g2_laboratory PASSED")


def test_full_protocol():
    """Test full protocol execution."""
    config = ProtocolConfig(delta_t=500, s_critical=0.6, alpha=0.05, n_permutations=1000)
    protocol = ThreeGateProtocol(config)
    
    data = {
        'gene1': {'pop_A': 1000, 'pop_B': 1100, 'pop_C': 1050, 'pop_D': 1150},
        'gene2': {'pop_A': 5000, 'pop_B': 8000, 'pop_C': 3000, 'pop_D': 9000},
    }
    
    results = protocol.analyze(data, is_laboratory_isolated=True)
    
    assert 'summary' in results, "Results should contain summary"
    assert 'decision' in results['summary'], "Summary should contain decision"
    assert results['summary']['decision'] in ['SUPPORT', 'FALSIFY', 'AMBIGUOUS'], \
        f"Invalid decision: {results['summary']['decision']}"
    print("✓ test_full_protocol PASSED")


def run_all_tests():
    """Run all tests."""
    print("\n" + "=" * 50)
    print("Running TCGS-SEQUENTION Protocol Tests")
    print("=" * 50 + "\n")
    
    tests = [
        test_synchrony_statistic_perfect,
        test_synchrony_statistic_partial,
        test_synchrony_statistic_none,
        test_permutation_pvalue_perfect,
        test_protocol_config,
        test_protocol_g1_analysis,
        test_protocol_g2_laboratory,
        test_full_protocol,
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except AssertionError as e:
            print(f"✗ {test.__name__} FAILED: {e}")
            failed += 1
        except Exception as e:
            print(f"✗ {test.__name__} ERROR: {e}")
            failed += 1
    
    print("\n" + "=" * 50)
    print(f"Results: {passed} passed, {failed} failed")
    print("=" * 50 + "\n")
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)
