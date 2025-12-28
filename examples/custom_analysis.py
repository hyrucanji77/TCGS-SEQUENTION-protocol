#!/usr/bin/env python3
"""
Example: Running the Protocol on Custom Data
=============================================

This script demonstrates how to run the three-gate protocol
on your own dataset.

Usage:
    python examples/custom_analysis.py
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from protocol import ThreeGateProtocol, ProtocolConfig

def main():
    print("=" * 60)
    print("TCGS-SEQUENTION Protocol: Custom Data Analysis")
    print("=" * 60)
    
    # =================================================================
    # STEP 1: Prepare your data
    # =================================================================
    # 
    # Your data should be a dictionary with structure:
    #   {gene_name: {population_name: first_detection_time}}
    #
    # Example:
    
    my_data = {
        'gene_A': {
            'population_1': 1000,
            'population_2': 1200,
            'population_3': 1100,
            'population_4': 1300,
            'population_5': 950,
        },
        'gene_B': {
            'population_1': 5000,
            'population_2': 7500,
            'population_3': 6000,
            'population_4': 8000,
            'population_5': 5500,
        },
        'gene_C': {
            'population_1': 3000,
            'population_2': 3100,
            'population_3': 2900,
            'population_4': 3050,
            # population_5 has no mutation in this gene
        },
    }
    
    print("\nInput data:")
    for gene, pops in my_data.items():
        print(f"  {gene}: {len(pops)} populations")
    
    # =================================================================
    # STEP 2: Configure the protocol
    # =================================================================
    #
    # Adjust parameters based on your system:
    #   - delta_t: Synchrony window (in your time units)
    #   - s_critical: Minimum synchrony fraction (0.6 = 60%)
    #   - alpha: Significance level (0.05 = 5%)
    #   - min_populations: Minimum populations required
    
    config = ProtocolConfig(
        delta_t=500.0,         # ±500 time units for synchrony
        s_critical=0.6,        # Require 60% populations synchronized
        alpha=0.05,            # 5% significance level
        n_permutations=10000,  # 10,000 permutations
        min_populations=4      # Need at least 4 populations
    )
    
    print(f"\nConfiguration:")
    print(f"  Synchrony window: ±{config.delta_t}")
    print(f"  S* threshold: {config.s_critical}")
    print(f"  Significance level: {config.alpha}")
    print(f"  Config hash: {config.config_hash}")
    
    # =================================================================
    # STEP 3: Run the protocol
    # =================================================================
    
    protocol = ThreeGateProtocol(config)
    
    # For laboratory-isolated populations (like LTEE), set:
    #   is_laboratory_isolated=True
    # This automatically passes Gate G2.
    
    results = protocol.analyze(
        my_data,
        is_laboratory_isolated=True  # Set to False for natural populations
    )
    
    # =================================================================
    # STEP 4: Interpret results
    # =================================================================
    
    print("\n" + "=" * 60)
    print("RESULTS")
    print("=" * 60)
    
    # Gate G1 results
    print("\n--- Gate G1: Synchrony ---")
    for gene, g1 in results['g1_results'].items():
        status = "✓ PASS" if g1.passes_g1 else "✗ FAIL"
        print(f"  {gene}: S*={g1.s_statistic:.3f}, p={g1.p_value:.4f} [{status}]")
    
    # Gate G2 result
    print("\n--- Gate G2: Independence ---")
    g2 = results['g2_result']
    status = "✓ PASS" if g2.passes_g2 else "✗ FAIL"
    print(f"  {g2.notes} [{status}]")
    
    # Gate G3 results
    print("\n--- Gate G3: Common-Cause Exclusion ---")
    for gene, g3 in results['g3_results'].items():
        status = "✓ PASS" if g3.passes_g3 else "✗ FAIL"
        print(f"  {gene}: retained={g3.synchrony_retained:.1%}, ΔAIC={g3.delta_aic:.1f} [{status}]")
    
    # Summary
    summary = results['summary']
    print("\n" + "=" * 60)
    print("FINAL DECISION")
    print("=" * 60)
    print(f"\n  Decision: {summary['decision']}")
    print(f"  {summary['interpretation']}")
    
    if summary['fully_passing']:
        print(f"\n  Genes passing all gates: {summary['fully_passing']}")
    
    print("\n" + "=" * 60)

if __name__ == "__main__":
    main()
