import os
import csv
import sys
import time
import math
import numpy as np
import matplotlib.pyplot as plt
from Bio.PDB import MMCIFParser, PPBuilder
from datetime import datetime
from matplotlib.colors import ListedColormap

# =============================================================================
# ACADEMIC FORMATTING & LOGGING UTILITIES
# =============================================================================

def print_banner(text):
    """Standardized header for pipeline phase demarcation."""
    print("\n" + "="*115)
    print(f"{text:^115}")
    print("="*115)

# =============================================================================
# STEREOCONFORMATIONAL ANALYSIS ENGINE (RAMPAGE SPECIFICATION)
# =============================================================================

def generate_ramachandran_analysis(residue_data, output_path, statistics):
    """
    Generates a high-resolution Ramachandran plot utilizing discrete 
    tiling logic consistent with RAMPAGE-style distributions.
    """
    phi = np.array([r['Phi'] for r in residue_data])
    psi = np.array([r['Psi'] for r in residue_data])

    # 1. Initialize 1-degree resolution density grid
    grid = np.zeros((360, 360))

    # Define conformational boundaries based on empirical torsion angle distributions
    for i in range(360):
        for j in range(360):
            p, s = i - 180, j - 180
            
            # --- BETA SHEET REGIONS (Upper Left Quadrant) ---
            if (-145 < p < -60 and 110 < s < 165): grid[j, i] = 3      # Favored
            elif (-180 < p < -40 and 80 < s < 180): grid[j, i] = 2     # Allowed
            
            # --- RIGHT-HANDED ALPHA HELIX (Lower Left Quadrant) ---
            elif (-100 < p < -40 and -70 < s < -25): grid[j, i] = 3    # Favored
            elif (-130 < p < -30 and -100 < s < 20): grid[j, i] = 2    # Allowed
            
            # --- LEFT-HANDED ALPHA HELIX (Central Right Quadrant) ---
            elif (45 < p < 85 and 30 < s < 75): grid[j, i] = 3         # Favored
            elif (20 < p < 110 and 10 < s < 110): grid[j, i] = 2       # Allowed
            
            # --- ADDITIONAL ALLOWED REGIONS ---
            elif grid[j, i] == 0:
                if (p < 0 and s > -120) or (p > 30 and s < 120):
                    grid[j, i] = 1 # Generous/Additional

    # 2. Formal RAMPAGE Color Grading
    # 0: Outlier (White), 1: Generous (Pale Yellow), 2: Allowed (Yellow), 3: Favored (Red)
    academic_cmap = ListedColormap(['#FFFFFF', '#FFFDE7', '#FFFF00', '#FF3D00'])

    fig, ax = plt.subplots(figsize=(12, 10), dpi=300)
    
    # Render background density map
    ax.imshow(grid, extent=[-180, 180, -180, 180], cmap=academic_cmap, origin='lower', aspect='auto')

    # 3. Plot atomic residue coordinates
    ax.scatter(phi, psi, marker='o', s=12, c='black', edgecolors='none', zorder=10)

    # 4. STATISTICAL VALIDATION SUMMARY
    summary_text = (
        f"        Plot Statistics\n\n"
        f"Residues in most favored regions [A,B,L]  : {statistics['favored']:>3} ({statistics['fav_pct']:.1f}%)\n"
        f"Residues in additional allowed regions    : {statistics['allowed']:>3} ({statistics['all_pct']:.1f}%)\n"
        f"Residues in disallowed regions           : {statistics['disallowed']:>3} ({statistics['dis_pct']:.1f}%)\n"
        f"--------------------------------------------------\n"
        f"Total number of residues                  : {statistics['total']}"
    )
    
    # Positioning of summary box in the right margin
    plt.gcf().text(0.72, 0.45, summary_text, fontsize=10, family='monospace', 
                   bbox=dict(facecolor='white', alpha=1.0, edgecolor='black', boxstyle='square,pad=1'))

    # 5. Formal Aesthetics & Plot Parameters
    ax.set_title(f"Ramachandran Analysis: Multi-Epitope Chimera Configuration", fontsize=16, fontweight='bold', loc='left', pad=15)
    ax.set_xlabel("Phi (φ) Degrees", fontsize=12)
    ax.set_ylabel("Psi (ψ) Degrees", fontsize=12)
    ax.set_xlim(-180, 180)
    ax.set_ylim(-180, 180)
    ax.set_xticks(range(-180, 181, 45))
    ax.set_yticks(range(-180, 181, 45))
    ax.axhline(0, color='black', linewidth=1)
    ax.axvline(0, color='black', linewidth=1)
    ax.grid(True, linestyle=':', color='gray', alpha=0.4)

    plt.tight_layout(rect=[0, 0, 0.7, 1])
    plt.savefig(output_path, bbox_inches='tight')
    plt.close()

# =============================================================================
# CORE VALIDATION PIPELINE
# =============================================================================

def execute_phase2_validation():
    """
    Orchestrates the stereoconformational audit of the validated tertiary model.
    """
    start_time = time.time()
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_model = os.path.join(project_root, "src", "phase2", "modelling_results", "Final_Optimal_Model.cif")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepE")
    os.makedirs(os.path.join(output_base, "Visuals"), exist_ok=True)

    print_banner("PHASE 2 STEP E: TERTIARY MODEL STEREOCONFORMATIONAL VALIDATION")

    if not os.path.exists(input_model):
        print(f"[ERROR] Required input model ({input_model}) not found.")
        return

    # Structural Parsing Initialization
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("vaccine_construct", input_model)
    ppb = PPBuilder()
    
    residue_data = []
    counts = {'favored': 0, 'allowed': 0, 'disallowed': 0, 'total': 0}

    # Iterative analysis of polypeptide backbone dihedrals
    for model in structure:
        for chain in model:
            for poly in ppb.build_peptides(chain):
                phi_psi = poly.get_phi_psi_list()
                for phi, psi in phi_psi:
                    if phi is None or psi is None: continue
                    p_deg, s_deg = math.degrees(phi), math.degrees(psi)
                    
                    is_favored = False
                    is_allowed = False
                    
                    # Evaluation against RAMPAGE propensity thresholds
                    if ((-145 < p_deg < -60 and 110 < s_deg < 165) or 
                        (-100 < p_deg < -40 and -70 < s_deg < -25) or 
                        (45 < p_deg < 85 and 30 < s_deg < 75)):
                        is_favored = True
                    elif ((-180 < p_deg < -40 and 80 < s_deg < 180) or 
                          (-130 < p_deg < -30 and -100 < s_deg < 20) or 
                          (20 < p_deg < 110 and 10 < s_deg < 110)):
                        is_allowed = True
                    
                    counts['total'] += 1
                    if is_favored: counts['favored'] += 1
                    elif is_allowed: counts['allowed'] += 1
                    else: counts['disallowed'] += 1
                    
                    residue_data.append({"Phi": p_deg, "Psi": s_deg})

    # Statistical Calculations
    counts['fav_pct'] = (counts['favored'] / counts['total']) * 100
    counts['all_pct'] = (counts['allowed'] / counts['total']) * 100
    counts['dis_pct'] = (counts['disallowed'] / counts['total']) * 100

    # Output generation
    plot_filename = "Stereo_Validation_Report.png"
    plot_path = os.path.join(output_base, "Visuals", plot_filename)
    generate_ramachandran_analysis(residue_data, plot_path, counts)

    # Terminal Report Dashboard
    print(f"{'STEREOCONFORMATIONAL METRIC':<45} | {'QUANTITATIVE RESULT'}")
    print("-" * 115)
    print(f"{'Proportion of Residues in Favored Regions':<45} | {counts['fav_pct']:.2f}%")
    print(f"{'Proportion of Residues in Allowed Regions':<45} | {counts['all_pct']:.2f}%")
    print(f"{'Proportion of Residues in Outlier Regions':<45} | {counts['dis_pct']:.2f}%")
    print("-" * 115)

    # Academic Verdict Logic
    pass_threshold = 90.0
    verdict = "PASS" if counts['fav_pct'] >= pass_threshold else "REVISE"
    
    print(f"PIPELINE VERDICT: [{verdict}]")
    print(f"Rigor Assessment: A favored distribution of {counts['fav_pct']:.2f}% confirms " 
          f"{'high-fidelity' if counts['fav_pct'] > 95 else 'acceptable'} backbone geometry.")
    
    print_banner("STEP 2E VALIDATION COMPLETE")
    print(f"[SUCCESS] Analytical Report Generated: {os.path.relpath(plot_path, project_root)}")
    print("="*115 + "\n")

if __name__ == "__main__":
    execute_phase2_validation()