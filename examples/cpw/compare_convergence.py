#!/usr/bin/env python3
"""
Compare AMR-only vs AMR+TMOP convergence for Palace simulations.

This script reads data from Palace postpro directories and plots
error indicator norm vs DOFs for convergence comparison.
"""

import json
import os
import glob
import matplotlib.pyplot as plt


def read_palace_data(postpro_dir):
    """
    Read Palace simulation data from a postpro directory.
    
    Returns list of dicts with 'iteration', 'dofs', 'error' for each AMR iteration.
    """
    data = []
    
    # Find all iteration directories
    iter_dirs = sorted(glob.glob(os.path.join(postpro_dir, 'iteration*')))
    
    for iter_dir in iter_dirs:
        # Extract iteration number from directory name
        iter_name = os.path.basename(iter_dir)
        iter_num = int(iter_name.replace('iteration', ''))
        
        # Read palace.json for DOFs
        palace_json = os.path.join(iter_dir, 'palace.json')
        if os.path.exists(palace_json):
            with open(palace_json, 'r') as f:
                palace_data = json.load(f)
                dofs = palace_data['Problem']['DegreesOfFreedom']
        else:
            continue
        
        # Read error-indicators.csv for error norm
        # Each iteration directory should have its own error indicators
        error_csv = os.path.join(iter_dir, 'error-indicators.csv')
        if os.path.exists(error_csv):
            with open(error_csv, 'r') as f:
                lines = f.readlines()
                if len(lines) >= 2:
                    # Parse the data line (second line)
                    values = lines[1].strip().split(',')
                    error_norm = float(values[0])
                else:
                    error_norm = None
        else:
            error_norm = None
        
        if error_norm is not None:
            data.append({
                'iteration': iter_num,
                'dofs': dofs,
                'error': error_norm
            })
    
    # Also check the final error-indicators.csv in the root postpro dir
    final_error_csv = os.path.join(postpro_dir, 'error-indicators.csv')
    final_palace_json = os.path.join(postpro_dir, 'palace.json')
    
    if os.path.exists(final_error_csv) and os.path.exists(final_palace_json):
        with open(final_palace_json, 'r') as f:
            final_data = json.load(f)
            final_dofs = final_data['Problem']['DegreesOfFreedom']
        
        with open(final_error_csv, 'r') as f:
            lines = f.readlines()
            if len(lines) >= 2:
                values = lines[1].strip().split(',')
                final_error = float(values[0])
                
                # Add final state if different from last iteration
                if not data or data[-1]['dofs'] != final_dofs:
                    data.append({
                        'iteration': len(data) + 1,
                        'dofs': final_dofs,
                        'error': final_error,
                        'final': True
                    })
    
    return data


def plot_comparison(amr_data, tmop_data, output_file='convergence_comparison.png'):
    """Plot error vs DOFs comparison."""
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if amr_data:
        amr_dofs = [d['dofs'] for d in amr_data]
        amr_error = [d['error'] for d in amr_data]
        ax.loglog(amr_dofs, amr_error, 'b-o', linewidth=2, markersize=8, label='AMR only')
    
    if tmop_data:
        tmop_dofs = [d['dofs'] for d in tmop_data]
        tmop_error = [d['error'] for d in tmop_data]
        ax.loglog(tmop_dofs, tmop_error, 'r-s', linewidth=2, markersize=8, label='AMR + TMOP')
    
    ax.set_xlabel('Degrees of Freedom', fontsize=12)
    ax.set_ylabel('Error Indicator Norm', fontsize=12)
    ax.set_title('AMR Convergence: With and Without TMOP Mesh Optimization', fontsize=14)
    ax.legend(fontsize=11)
    ax.grid(True, which='both', linestyle='--', alpha=0.7)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Saved plot to {output_file}")
    plt.close()


def print_data_table(name, data):
    """Print data as a formatted table."""
    print(f"\n{name}:")
    print("-" * 50)
    print(f"{'Iteration':>10} {'DOFs':>12} {'Error Norm':>15}")
    print("-" * 50)
    for d in data:
        print(f"{d['iteration']:>10} {d['dofs']:>12} {d['error']:>15.4e}")


if __name__ == '__main__':
    # Read data from postpro directories
    script_dir = os.path.dirname(os.path.abspath(__file__))
    
    amr_only_dir = os.path.join(script_dir, 'postpro', 'amr_only')
    tmop_dir = os.path.join(script_dir, 'postpro', 'tmop_test')
    
    amr_data = []
    tmop_data = []
    
    if os.path.exists(amr_only_dir):
        amr_data = read_palace_data(amr_only_dir)
        print_data_table("AMR Only", amr_data)
    else:
        print(f"Warning: AMR-only directory not found: {amr_only_dir}")
    
    if os.path.exists(tmop_dir):
        tmop_data = read_palace_data(tmop_dir)
        print_data_table("AMR + TMOP", tmop_data)
    else:
        print(f"Warning: TMOP directory not found: {tmop_dir}")
    
    if amr_data or tmop_data:
        output_file = os.path.join(script_dir, 'convergence_comparison.png')
        plot_comparison(amr_data, tmop_data, output_file)
        print("\nComparison complete!")
    else:
        print("\nNo data found to plot.")
