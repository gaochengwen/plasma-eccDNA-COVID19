#!/usr/bin/env python3
"""
Reproduce Figure 3 Panel D from /home/gao/eccDNA/revise/Figure 3.pdf
Sorted strictly by recurrence frequency from high to low (top to bottom).
Outputs Figure_3D.pdf, Figure_3D.png, and Figure_3D.svg to /home/gao/eccDNA/revise/.
"""

import os
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Wedge, Rectangle

def parse_pos(ecc_str: str):
    chrom, pos = ecc_str.split(':')
    c_num = int(chrom.replace('chr', '')) if chrom.replace('chr', '').isdigit() else 99
    start_pos = int(pos.split('_')[0]) if '_' in pos else int(pos.split('-')[0])
    return (c_num, start_pos)

def main():
    # Set up matplotlib style for Nature publication standards
    mpl.rcParams['font.sans-serif'] = 'Arial'
    mpl.rcParams['font.family'] = 'sans-serif'
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['axes.edgecolor'] = '#000000'
    mpl.rcParams['axes.linewidth'] = 0.8

    COVID_RED = '#C82829'  # Nature Red
    GREY_BG = '#E5E5E5'    # Light grey for undetected

    source_dir = '/home/gao/eccDNA/revise/source_data'
    output_dir = '/home/gao/eccDNA/revise'

    f_matrix = os.path.join(source_dir, 'covid_specific_exact.matrix.min10.tsv')
    df_matrix = pd.read_csv(f_matrix, sep='\t', index_col=0)

    f_sorted = os.path.join(source_dir, 'Figure_3D_sorted_recurrent_circles.tsv')
    df_sorted = pd.read_csv(f_sorted, sep='\t')

    # Primary sort: recurrence frequency descending
    # Secondary sort: chromosome number and start coordinate ascending
    df_sorted['chrom_num'] = df_sorted['eccDNA'].apply(lambda x: parse_pos(x)[0])
    df_sorted['start_pos'] = df_sorted['eccDNA'].apply(lambda x: parse_pos(x)[1])
    df_sorted = df_sorted.sort_values(
        by=['recurrence', 'chrom_num', 'start_pos'],
        ascending=[False, True, True]
    ).reset_index(drop=True)

    loci_ordered = df_sorted['eccDNA'].tolist()
    recurrences = df_sorted['recurrence'].tolist()
    samples = df_matrix.columns.tolist()

    n_loci = len(loci_ordered)
    n_samples = len(samples)

    fig, ax = plt.subplots(figsize=(7.5, 7.5), dpi=300)
    ax.set_aspect('equal')
    ax.axis('off')

    cx, cy = 0, 0
    r_inner_base = 2.2
    r_outer_base = 5.2
    track_width = (r_outer_base - r_inner_base) / n_loci
    r_gap = 0.012

    total_angle_deg = 270.0
    angle_per_sample = total_angle_deg / n_samples
    a_gap = 0.5  # angular gap in degrees

    # 1. Concentric rings (track heatmaps)
    # Locus 0 (highest frequency = 14) at outermost ring (r_idx = 26)
    # Locus 26 (lowest frequency = 10) at innermost ring (r_idx = 0)
    for i_locus in range(n_loci):
        r_idx = n_loci - 1 - i_locus
        r1 = r_inner_base + r_idx * track_width + r_gap / 2
        r2 = r_inner_base + (r_idx + 1) * track_width - r_gap / 2
        
        locus_name = loci_ordered[i_locus]
        
        for j_sample in range(n_samples):
            sample_name = samples[j_sample]
            val = df_matrix.loc[locus_name, sample_name] if locus_name in df_matrix.index else 0
            
            th_start = - j_sample * angle_per_sample - a_gap / 2
            th_end = - (j_sample + 1) * angle_per_sample + a_gap / 2
            
            color = COVID_RED if val == 1 else GREY_BG
            w = Wedge((cx, cy), r2, th_end, th_start, width=(r2-r1), facecolor=color, edgecolor='none')
            ax.add_patch(w)

    # 2. Sample labels around perimeter
    r_label = r_outer_base + 0.18
    for j_sample in range(n_samples):
        sample_name = samples[j_sample]
        th_mid_deg = - (j_sample + 0.5) * angle_per_sample
        th_mid_rad = np.radians(th_mid_deg)
        
        lx = cx + r_label * np.cos(th_mid_rad)
        ly = cy + r_label * np.sin(th_mid_rad)
        
        rot_deg = th_mid_deg
        if -90 < th_mid_deg <= 0:
            rot_deg = th_mid_deg
            ha = 'left'
            va = 'center'
        elif -180 < th_mid_deg <= -90:
            rot_deg = th_mid_deg + 180
            ha = 'right'
            va = 'center'
        else:  # -270 <= th_mid_deg <= -180
            rot_deg = th_mid_deg + 180
            ha = 'right'
            va = 'center'
            
        ax.text(lx, ly, sample_name, fontsize=5.5, rotation=rot_deg, ha=ha, va=va, rotation_mode='anchor')

    # 3. Bar chart in top-right quadrant
    bar_x_start = 0.25
    bar_scale = 0.165
    bar_h = track_width * 0.72

    for i_locus in range(n_loci):
        r_idx = n_loci - 1 - i_locus
        r_mid = r_inner_base + (r_idx + 0.5) * track_width
        
        rec = recurrences[i_locus]
        locus_name = loci_ordered[i_locus].replace('_', '–')
        
        bar_width = rec * bar_scale
        y_bottom = r_mid - bar_h / 2
        
        rect = Rectangle((bar_x_start, y_bottom), bar_width, bar_h, facecolor=COVID_RED, edgecolor='none')
        ax.add_patch(rect)
        
        ax.text(bar_x_start + bar_width + 0.08, r_mid, locus_name, fontsize=4.8, va='center', ha='left')

    # 4. Scale axis at bottom of bar chart
    r_min_ring = r_inner_base + 0.5 * track_width
    axis_y = r_min_ring - 0.22

    ax.plot([bar_x_start, bar_x_start + 16 * bar_scale], [axis_y, axis_y], color='black', lw=0.8)
    for tick_val in [0, 8, 16]:
        tx = bar_x_start + tick_val * bar_scale
        ax.plot([tx, tx], [axis_y, axis_y - 0.07], color='black', lw=0.8)
        ax.text(tx, axis_y - 0.12, str(tick_val), fontsize=5.5, ha='center', va='top')

    # Panel label 'd'
    ax.text(-6.6, 6.2, 'd', fontsize=12, fontweight='bold', ha='left', va='top')

    ax.set_xlim(-7.0, 7.0)
    ax.set_ylim(-7.0, 7.0)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Figure_3D.pdf'), bbox_inches='tight')
    plt.savefig(os.path.join(output_dir, 'Figure_3D.png'), bbox_inches='tight', dpi=300)
    plt.savefig(os.path.join(output_dir, 'Figure_3D.svg'), bbox_inches='tight')
    print('Figure_3D successfully generated and saved to', output_dir)

if __name__ == '__main__':
    main()
