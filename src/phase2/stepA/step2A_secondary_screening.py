import os
import csv
import sys
import time
import numpy as np
from datetime import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =============================================================================
# ACADEMIC FORMATTING & BENCHMARKING UTILITIES
# =============================================================================

def print_banner(text):
    """Prints a formal academic section header."""
    print("\n" + "="*125)
    print(f"{text:^125}")
    print("="*125)

def format_time(seconds):
    """Formats execution time into a readable MM:SS string."""
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

# =============================================================================
# CORE SECONDARY SCREENING ENGINE
# =============================================================================

def run_step2a_comprehensive_screening():
    start_time = time.time()

    # 1. DYNAMIC DIRECTORY RESOLUTION
    # Location: MpoxHIV_Vaccine/src/phase2/stepA/step2A_secondary_screening.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1G", "Variants")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA")
    
    raw_out_dir = os.path.join(output_base, "Raw")
    filt_out_dir = os.path.join(output_base, "Filtered")
    
    os.makedirs(raw_out_dir, exist_ok=True)
    os.makedirs(filt_out_dir, exist_ok=True)

    if not os.path.exists(input_folder):
        print(f"\n[FATAL ERROR] Variant directory not found: {input_folder}")
        sys.exit(1)

    # 2. EXPERIMENTAL INITIALIZATION
    print_banner("PHASE 2 STEP A: COMPREHENSIVE SECONDARY SCREENING DASHBOARD")
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Input Source        : {input_folder}")
    print(f"[INFO] Screening Standards : Stability < 40 | Solubility (GRAVY) < 0 | Toxicity < 15 Cys")
    print("-" * 125)

    raw_results = []
    # Hydrophobicity scale for Fold-Score (Positional Variance)
    hydro_scale = {'A': 1.8, 'L': 3.8, 'I': 4.5, 'V': 4.2, 'F': 2.8, 'M': 1.9, 'C': 2.5, 
                   'G': -0.4, 'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 
                   'N': -3.5, 'Q': -3.5, 'D': -3.5, 'E': -3.5, 'K': -3.9, 'R': -4.5, 'H': -3.2}

    fasta_files = sorted([f for f in os.listdir(input_folder) if f.endswith(".fasta")])
    
    if not fasta_files:
        print("[WARNING] No FASTA variants found in Phase 1G. Please check Step 1G execution.")
        return

    # 3. ANALYTICAL EXECUTION
    for filename in fasta_files:
        with open(os.path.join(input_folder, filename), 'r') as f:
            lines = f.readlines()
            seq = "".join([l.strip() for l in lines if not l.startswith(">")]).upper()
        
        # Biochemical Calculations (Compositional + Sequential)
        ana = ProteinAnalysis(seq)
        instability_idx = ana.instability_index()
        gravy_val = ana.gravy()
        cys_count = seq.count('C')
        qn_ratio = (seq.count('Q') + seq.count('N')) / len(seq)
        fold_score = np.var([hydro_scale.get(aa, 0) for aa in seq])
        
        # Categorical Pass/Fail Logic
        s_stab = "PASS" if instability_idx < 40 else "FAIL"
        s_sol  = "PASS" if gravy_val < 0 else "FAIL"
        s_tox  = "PASS" if cys_count < 15 else "FAIL"
        s_all  = "PASS" if qn_ratio < 0.4 else "FAIL"
        
        is_viable = all(x == "PASS" for x in [s_stab, s_sol, s_tox, s_all])

        raw_results.append({
            "Variant": filename,
            "Len": len(seq),
            "STAB_IDX": round(instability_idx, 2),
            "STAB": s_stab,
            "SOL": s_sol,
            "TOX": s_tox,
            "ALLER": s_all,
            "Fold_Score": round(fold_score, 4),
            "GRAVY": round(gravy_val, 4),
            "Viable": "YES" if is_viable else "NO"
        })

    # 4. TERMINAL OUTPUT: RAW DATA TABLE
    print(f"{'VARIANT':<22} | {'LEN':<4} | {'STAB':<6} | {'SOL':<6} | {'TOX':<6} | {'ALLER':<6} | {'FOLD-SCORE':<10} | {'VIABLE'}")
    print("-" * 125)
    for r in raw_results:
        print(f"{r['Variant']:<22} | {r['Len']:<4} | {r['STAB']:<6} | {r['SOL']:<6} | {r['TOX']:<6} | {r['ALLER']:<6} | {r['Fold_Score']:<10} | {r['Viable']}")

    # 5. DATA FILTERING & RANKING
    # Isolate only those that passed all categories
    filtered_list = [r for r in raw_results if r['Viable'] == "YES"]
    # Rank by Stability (Most stable at the top)
    filtered_list.sort(key=lambda x: x['STAB_IDX'])

    # 6. TERMINAL OUTPUT: FILTERED & RANKED TABLE
    print_banner("STEP 2A: FILTERED OPTIMAL CANDIDATES (RANKED BY STABILITY INDEX)")
    if filtered_list:
        print(f"{'RANK':<5} | {'VARIANT':<22} | {'STABILITY INDEX':<18} | {'GRAVY':<10} | {'STATUS'}")
        print("-" * 125)
        for i, f in enumerate(filtered_list):
            status_label = "[OPTIMAL]" if i == 0 else "[ACCEPTED]"
            print(f"{i+1:<5} | {f['Variant']:<22} | {f['STAB_IDX']:<18} | {f['GRAVY']:<10} | {status_label}")
    else:
        print(f"{'--- NO VARIANTS PASSED THE COMPREHENSIVE SCREENING CRITERIA ---':^125}")

    # 7. DATA EXPORT (Dual CSV Logging)
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    raw_path = os.path.join(raw_out_dir, f"Step2A_Raw_Physicochemical_{ts}.csv")
    filt_path = os.path.join(filt_out_dir, f"Step2A_Filtered_Ranked_{ts}.csv")

    with open(raw_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=raw_results[0].keys())
        writer.writeheader()
        writer.writerows(raw_results)

    if filtered_list:
        with open(filt_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=filtered_list[0].keys())
            writer.writeheader()
            writer.writerows(filtered_list)

    # 8. FINAL BENCHMARK SUMMARY
    exec_time = format_time(time.time() - start_time)
    print_banner("PHASE 2 STEP A COMPLETE")
    print(f"[SUCCESS] Total Analyzed      : {len(raw_results)} variants")
    print(f"[SUCCESS] Viable Candidates   : {len(filtered_list)}")
    print(f"[SUCCESS] Execution Time      : {exec_time}")
    print(f"[INFO] Raw Log: {os.path.relpath(raw_path, project_root)}")
    print(f"[INFO] Filtered Log: {os.path.relpath(filt_path, project_root)}")
    print("="*125 + "\n")

if __name__ == "__main__":
    run_step2a_comprehensive_screening()