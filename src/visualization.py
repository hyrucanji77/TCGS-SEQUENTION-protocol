#!/usr/bin/env python3
"""
TCGS-SEQUENTION: Visualization Functions
=========================================

Plotting functions for protocol results.

Author: Henry Arellano-Peña
Framework: TCGS-SEQUENTION
License: CC BY 4.0
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from typing import Dict, List, Optional, Any


def plot_synchrony_overview(g1_results: Dict[str, Any],
                            delta_t: float,
                            output_path: Optional[str] = None,
                            figsize: Tuple[int, int] = (12, 8)) -> plt.Figure:
    """
    Create overview plot of G1 synchrony results.
    
    Parameters
    ----------
    g1_results : dict
        Dictionary of G1 results by gene
    delta_t : float
        Synchrony window size
    output_path : str, optional
        Path to save figure
    figsize : tuple
        Figure size
        
    Returns
    -------
    fig : Figure
        Matplotlib figure
    """
    genes = list(g1_results.keys())
    n_genes = len(genes)
    
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    
    # Panel 1: S* statistics
    ax1 = axes[0, 0]
    s_values = [g1_results[g].s_statistic for g in genes]
    passes = [g1_results[g].passes_g1 for g in genes]
    colors = ['green' if p else 'red' for p in passes]
    
    bars = ax1.barh(range(n_genes), s_values, color=colors, alpha=0.7, edgecolor='black')
    ax1.axvline(x=0.6, color='red', linestyle='--', linewidth=2, label='S* threshold')
    ax1.set_yticks(range(n_genes))
    ax1.set_yticklabels(genes)
    ax1.set_xlabel('S* Statistic')
    ax1.set_title('Gate G1: Synchrony Statistics')
    ax1.legend()
    ax1.set_xlim(0, 1.1)
    
    # Panel 2: p-values
    ax2 = axes[0, 1]
    p_values = [g1_results[g].p_value for g in genes]
    colors = ['green' if p <= 0.05 else 'red' for p in p_values]
    
    ax2.barh(range(n_genes), [-np.log10(p + 1e-10) for p in p_values], 
             color=colors, alpha=0.7, edgecolor='black')
    ax2.axvline(x=-np.log10(0.05), color='red', linestyle='--', linewidth=2, label='α = 0.05')
    ax2.set_yticks(range(n_genes))
    ax2.set_yticklabels(genes)
    ax2.set_xlabel('-log₁₀(p-value)')
    ax2.set_title('Gate G1: Statistical Significance')
    ax2.legend()
    
    # Panel 3: Population counts
    ax3 = axes[1, 0]
    n_pops = [g1_results[g].n_populations for g in genes]
    ax3.barh(range(n_genes), n_pops, color='steelblue', alpha=0.7, edgecolor='black')
    ax3.set_yticks(range(n_genes))
    ax3.set_yticklabels(genes)
    ax3.set_xlabel('Number of Populations')
    ax3.set_title('Populations with Parallel Evolution')
    
    # Panel 4: Summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    n_passing = sum(passes)
    summary_text = f"""
    GATE G1 SUMMARY
    ━━━━━━━━━━━━━━━━━━━━
    Genes analyzed: {n_genes}
    Genes passing: {n_passing}
    Pass rate: {100*n_passing/n_genes:.1f}%
    
    Criteria:
    • S* ≥ 0.6
    • p ≤ 0.05
    • Δt = {delta_t:,.0f} generations
    """
    ax4.text(0.1, 0.5, summary_text, fontsize=12, fontfamily='monospace',
             verticalalignment='center', transform=ax4.transAxes)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_timing_heatmap(gene_data: Dict[str, Dict[str, float]],
                        populations: List[str],
                        output_path: Optional[str] = None,
                        figsize: Tuple[int, int] = (14, 10)) -> plt.Figure:
    """
    Create heatmap of mutation timing across genes and populations.
    
    Parameters
    ----------
    gene_data : dict
        Data structure: {gene: {population: time}}
    populations : list
        List of population names
    output_path : str, optional
        Path to save figure
    figsize : tuple
        Figure size
        
    Returns
    -------
    fig : Figure
        Matplotlib figure
    """
    genes = list(gene_data.keys())
    n_genes = len(genes)
    n_pops = len(populations)
    
    # Build timing matrix
    timing_matrix = np.full((n_genes, n_pops), np.nan)
    for i, gene in enumerate(genes):
        for j, pop in enumerate(populations):
            if pop in gene_data[gene]:
                timing_matrix[i, j] = gene_data[gene][pop]
    
    fig, ax = plt.subplots(figsize=figsize)
    
    im = ax.imshow(timing_matrix / 1000, aspect='auto', cmap='viridis')
    
    ax.set_xticks(range(n_pops))
    ax.set_xticklabels(populations, rotation=45, ha='right')
    ax.set_yticks(range(n_genes))
    ax.set_yticklabels(genes)
    
    ax.set_xlabel('Population', fontsize=12)
    ax.set_ylabel('Gene', fontsize=12)
    ax.set_title('Mutation Timing Heatmap (×1000 generations)', fontsize=14)
    
    cbar = plt.colorbar(im, ax=ax, label='Generations (×1000)')
    
    # Add text annotations
    for i in range(n_genes):
        for j in range(n_pops):
            if not np.isnan(timing_matrix[i, j]):
                text = ax.text(j, i, f'{timing_matrix[i,j]/1000:.1f}',
                              ha='center', va='center', color='white', fontsize=8)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_decision_summary(results: Dict[str, Any],
                          output_path: Optional[str] = None,
                          figsize: Tuple[int, int] = (10, 8)) -> plt.Figure:
    """
    Create decision summary visualization.
    
    Parameters
    ----------
    results : dict
        Complete protocol results
    output_path : str, optional
        Path to save figure
    figsize : tuple
        Figure size
        
    Returns
    -------
    fig : Figure
        Matplotlib figure
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.axis('off')
    
    summary = results.get('summary', {})
    
    # Gate boxes
    gate_data = [
        ('G1\nSynchrony', len(summary.get('g1_passing', [])) > 0,
         f"{len(summary.get('g1_passing', []))} genes"),
        ('G2\nIndependence', summary.get('g2_passes', False),
         'PASS' if summary.get('g2_passes', False) else 'FAIL'),
        ('G3\nCommon-Cause', len(summary.get('g3_passing', [])) > 0,
         f"{len(summary.get('g3_passing', []))} genes"),
    ]
    
    for i, (label, passes, detail) in enumerate(gate_data):
        color = 'lightgreen' if passes else 'lightcoral'
        rect = mpatches.FancyBboxPatch(
            (0.1 + i * 0.3, 0.55), 0.25, 0.35,
            boxstyle="round,pad=0.02",
            facecolor=color, edgecolor='black', linewidth=2
        )
        ax.add_patch(rect)
        ax.text(0.225 + i * 0.3, 0.82, label, ha='center', va='center',
                fontsize=14, fontweight='bold')
        ax.text(0.225 + i * 0.3, 0.68, 'PASS' if passes else 'FAIL',
                ha='center', va='center', fontsize=16, fontweight='bold')
        ax.text(0.225 + i * 0.3, 0.58, detail, ha='center', va='center', fontsize=10)
    
    # Decision box
    decision = summary.get('decision', 'UNKNOWN')
    decision_colors = {
        'SUPPORT': 'gold',
        'FALSIFY': 'lightcoral',
        'AMBIGUOUS': 'lightyellow',
        'UNKNOWN': 'lightgray'
    }
    
    rect = mpatches.FancyBboxPatch(
        (0.2, 0.1), 0.6, 0.35,
        boxstyle="round,pad=0.02",
        facecolor=decision_colors.get(decision, 'lightgray'),
        edgecolor='black', linewidth=3
    )
    ax.add_patch(rect)
    ax.text(0.5, 0.35, 'PROTOCOL DECISION', ha='center', va='center',
            fontsize=14, fontweight='bold')
    ax.text(0.5, 0.22, decision, ha='center', va='center',
            fontsize=24, fontweight='bold')
    
    # Passing genes
    fully_passing = summary.get('fully_passing', [])
    if fully_passing:
        ax.text(0.5, 0.02, f'Genes passing all gates: {", ".join(fully_passing)}',
                ha='center', va='center', fontsize=10, style='italic')
    
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    plt.title('TCGS-SEQUENTION Protocol Results', fontsize=16, fontweight='bold', y=0.98)
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


def plot_gene_detail(gene: str,
                     times: Dict[str, float],
                     delta_t: float,
                     output_path: Optional[str] = None,
                     figsize: Tuple[int, int] = (10, 6)) -> plt.Figure:
    """
    Create detailed visualization for a single gene.
    
    Parameters
    ----------
    gene : str
        Gene name
    times : dict
        Dictionary mapping populations to emergence times
    delta_t : float
        Synchrony window size
    output_path : str, optional
        Path to save figure
    figsize : tuple
        Figure size
        
    Returns
    -------
    fig : Figure
        Matplotlib figure
    """
    pops = list(times.keys())
    values = list(times.values())
    n_pops = len(pops)
    
    fig, ax = plt.subplots(figsize=figsize)
    
    # Horizontal bar plot
    y_pos = np.arange(n_pops)
    bars = ax.barh(y_pos, values, color='steelblue', alpha=0.8, edgecolor='black')
    
    # Reference line (median)
    median_time = np.median(values)
    ax.axvline(x=median_time, color='red', linestyle='--', linewidth=2,
               label=f'Median: {median_time:,.0f}')
    
    # Synchrony window
    ax.axvspan(median_time - delta_t, median_time + delta_t,
               alpha=0.2, color='green', label=f'±{delta_t:,.0f} window')
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(pops)
    ax.set_xlabel('Generation of First Detection', fontsize=12)
    ax.set_ylabel('Population', fontsize=12)
    ax.set_title(f'{gene}: Mutation Timing Across {n_pops} Populations',
                 fontsize=14, fontweight='bold')
    ax.legend(loc='lower right')
    
    # Add value labels
    for i, (p, v) in enumerate(zip(pops, values)):
        ax.text(v + max(values) * 0.01, i, f'{v:,.0f}', va='center', fontsize=9)
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
    
    return fig


# Type hint import
from typing import Tuple
