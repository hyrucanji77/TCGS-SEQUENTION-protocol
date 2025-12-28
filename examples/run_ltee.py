#!/usr/bin/env python3
"""
Example: Running the LTEE Analysis
==================================

This script demonstrates how to run the three-gate protocol
on the E. coli LTEE dataset.

Usage:
    python examples/run_ltee.py
"""

import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

from ltee_analysis import run_full_protocol, ProtocolConfig, export_results

def main():
    print("=" * 60)
    print("TCGS-SEQUENTION Protocol: LTEE Analysis")
    print("=" * 60)
    
    # Create configuration with preregistered parameters
    config = ProtocolConfig(
        delta_t=2000.0,       # ±2000 generation synchrony window
        s_critical=0.6,       # Minimum 60% populations synchronized
        alpha=0.05,           # 5% significance level
        n_permutations=10000  # 10,000 permutation tests
    )
    
    print(f"\nConfiguration hash: {config.config_hash}")
    print(f"Timestamp: {config.timestamp}")
    
    # Run the full three-gate protocol
    print("\nRunning protocol...")
    result = run_full_protocol(config)
    
    # Print results
    print("\n" + "=" * 60)
    print("RESULTS SUMMARY")
    print("=" * 60)
    
    print(f"\nOverall Decision: {result.overall_decision}")
    print(f"Interpretation: {result.interpretation}")
    
    print(f"\nGate G1 (Synchrony): {len(result.g1_passing_genes)}/{len(result.g1_results)} genes pass")
    print(f"Gate G2 (Independence): {'PASS' if result.g2_result.passes_g2 else 'FAIL'}")
    print(f"Gate G3 (Common-Cause): {len(result.g3_passing_genes)}/{len(result.g1_passing_genes)} genes pass")
    
    print(f"\nGenes passing ALL gates: {result.fully_passing_genes}")
    
    # Export results
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'output')
    os.makedirs(output_dir, exist_ok=True)
    
    export_results(result, output_dir)
    print(f"\nResults exported to: {output_dir}")
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)

if __name__ == "__main__":
    main()
