import os
import json
import shutil
import csv
import time
import numpy as np
from datetime import datetime

# =============================================================================
# ACADEMIC FORMATTING & UTILITIES
# =============================================================================

def print_banner(text):
    print("\n" + "="*125)
    print(f"{text:^125}")
    print("="*125)

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

# =============================================================================
# HEURISTIC DATA FINDER (The "LMAO" Killer)
# =============================================================================

def find_confidence_list(obj):
    """
    Recursively scans the JSON object for the longest list of floats.
    In AlphaFold JSONs, the pLDDT is always the primary large numeric array.
    """
    max_list = []
    
    if isinstance(obj, dict):
        for k, v in obj.items():
            found = find_confidence_list(v)
            if len(found) > len(max_list):
                max_list = found
    elif isinstance(obj, list):
        # Check if this is a list of numbers
        if len(obj) > 0 and isinstance(obj[0], (int, float)):
            return obj
        # Otherwise keep searching nested lists
        for item in obj:
            found = find_confidence_list(item)
            if len(found) > len(max_list):
                max_list = found
                
    return max_list

# =============================================================================
# CORE PARSER ENGINE
# =============================================================================

def run_alphafold3_server_parser():
    start_time = time.time()
    
    # 1. Path Resolution
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    target_dir = os.path.join(project_root, "src", "phase2", "modelling_results")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepD")
    
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    os.makedirs(raw_out_dir := raw_dir, exist_ok=True)
    os.makedirs(filt_out_dir := filt_dir, exist_ok=True)

    print_banner("PHASE 2 STEP D.2: ALPHAFOLD3 MULTI-SEED CONFIDENCE AUDIT")
    print(f"[INFO] Audit Strategy  : Heuristic Numeric List Detection (Key-Agnostic)")
    print(f"[INFO] Target Directory : {os.path.relpath(target_dir, project_root)}")
    print("-" * 125)

    results = []
    
    # 2. Iterate through Seeds 0-4
    for i in range(5):
        summary_file = None
        full_data_file = None
        cif_file = None
        
        # Match files by seed index
        for f in os.listdir(target_dir):
            if f"summary_confidences_{i}.json" in f: summary_file = f
            if f"full_data_{i}.json" in f: full_data_file = f
            if f"model_{i}.cif" in f: cif_file = f

        if not summary_file or not full_data_file: continue

        try:
            # Get PTM and Rank Score from Summary
            with open(os.path.join(target_dir, summary_file), 'r') as f:
                s_data = json.load(f)
            
            # Get pLDDT from Full Data using the Heuristic Hunter
            with open(os.path.join(target_dir, full_data_file), 'r') as f:
                f_data = json.load(f)

            plddt_scores = find_confidence_list(f_data)
            avg_plddt = np.mean(plddt_scores) if plddt_scores else 0.0
            
            # Extract standard AF3 metrics
            ptm = s_data.get('ptm', 0.0)
            rank_score = s_data.get('ranking_score', 0.0)

            results.append({
                "Seed": i,
                "Avg_pLDDT": round(avg_plddt, 2),
                "pTM": round(ptm, 4),
                "Rank_Score": round(rank_score, 4),
                "CIF": cif_file
            })
        except Exception as e:
            print(f"[WARNING] Seed {i} failed to parse: {e}")

    if not results:
        print("[ERROR] No valid data found. Check your file naming.")
        return

    # 3. RANKING (Highest Rank Score Wins)
    results.sort(key=lambda x: x['Rank_Score'], reverse=True)
    winner = results[0]

    # 4. TERMINAL DASHBOARD
    print(f"{'MODEL SEED':<15} | {'AVG pLDDT %':<15} | {'pTM SCORE':<15} | {'RANK SCORE':<15} | {'STATUS'}")
    print("-" * 125)
    for r in results:
        status = "OPTIMAL" if r == winner else "REJECTED"
        print(f"Seed {r['Seed']:<10} | {r['Avg_pLDDT']:>12.2f}% | {r['pTM']:>15.3f} | {r['Rank_Score']:>15.3f} | {status}")

 # 5. PREPARE BEST MODEL (Corrected Pathing)
    if winner.get('CIF'):
        source_path = os.path.join(target_dir, winner['CIF'])
        # Destination should be the 'Filtered' folder in Step_Outputs
        destination_path = os.path.join(filt_out_dir, "Final_Optimal_Model.cif")
        
        shutil.copy(source_path, destination_path)
        
        print("-" * 125)
        print(f"[SUCCESS] High-Confidence Candidate: {winner['CIF']}")
        print(f"[SUCCESS] Final Model Exported to  : {os.path.relpath(destination_path, project_root)}")

        
    # 6. CSV LOGGING
    ts = datetime.now().strftime("%Y%m%d_%H%M")
    with open(os.path.join(raw_out_dir, f"AF3_Raw_Audit_{ts}.csv"), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys()); w.writeheader(); w.writerows(results)
    with open(os.path.join(filt_out_dir, f"AF3_Winner_Report_{ts}.csv"), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=results[0].keys()); w.writeheader(); w.writerow(winner)

    print_banner("STEP D.2 AUDIT COMPLETE")
    print(f"[SUCCESS] Execution Time : {format_time(time.time() - start_time)}")
    print(f"[INFO] Files saved to: Step_Outputs/Phase2/StepD/")

if __name__ == "__main__":
    run_alphafold3_server_parser()