import os
import sys
import time
import re
import csv
from datetime import datetime
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def format_time(seconds):
    """Formats time in seconds to a readable MM:SS format."""
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def run_step1b_final_dual_output():
    start_time = time.time()

    # 1. DYNAMIC DIRECTORY RESOLUTION
    # Current location: [ROOT]/src/phase1/stepB/step1B_stability.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Traverse up 3 levels to reach the project root
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1A")
    output_path = os.path.join(project_root, "Step_Outputs", "Phase1B")
    
    # Create subfolders for better organization
    raw_out = os.path.join(output_path, "Raw_Stability")
    filt_out = os.path.join(output_path, "Filtered_Stability")
    os.makedirs(raw_out, exist_ok=True)
    os.makedirs(filt_out, exist_ok=True)
    
    if not os.path.exists(input_folder):
        print("\n[ERROR] Phase 1A output directory not found. Please execute Step 1A first.")
        return

    fasta_files =[f for f in os.listdir(input_folder) if f.endswith(".fasta")]
    total_files = len(fasta_files)

    if total_files == 0:
        print("\n[ERROR] No FASTA files found in Phase 1A directory.")
        return

    # 2. EXPERIMENTAL PREVIEW & HEADER
    print("\n" + "="*80)
    print(f"{'PHASE 1B: PHYSICOCHEMICAL STABILITY ANALYSIS (PROTPARAM)':^80}")
    print("="*80)
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Input Directory     : {input_folder}")
    print(f"[INFO] Target Files        : {total_files} FASTA sequences")
    print(f"[INFO] Stability Threshold : Instability Index < 40.0")
    print("-" * 80)

    # 3. MATHEMATICAL EVALUATION MATRIX
    raw_list = []
    filtered_list =[]
    
    sys.stdout.write("[PROCESS] Initiating ExPASy analytical engine...\n\n")
    
    for i, filename in enumerate(sorted(fasta_files)):
        file_path = os.path.join(input_folder, filename)
        variant_id = filename.replace(".fasta", "")
        target_group = filename.split('_Var')[0]
        
        # Dynamic terminal preview string
        elapsed = format_time(time.time() - start_time)
        status_str = f"[ EVAL ] Record {i+1:03d}/{total_files:03d} | Target: {target_group:<12} | Elapsed: {elapsed}"
        sys.stdout.write(f"\r{status_str:<80}")
        sys.stdout.flush()
        
        try:
            for record in SeqIO.parse(file_path, "fasta"):
                clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', str(record.seq).upper())
                if not clean_seq: continue

                # Physicochemical Calculations
                analysis = ProteinAnalysis(clean_seq)
                mw = analysis.molecular_weight() / 1000
                gravy = analysis.gravy()
                idx = analysis.instability_index()
                
                status = "STABLE" if idx < 40 else "UNSTABLE"
                
                data_row = {
                    "Variant_ID": variant_id,
                    "Target": target_group,
                    "MW_kDa": round(mw, 2),
                    "GRAVY": round(gravy, 2),
                    "Instability_Index": round(idx, 2),
                    "Status": status
                }

                # ALWAYS add to Raw log
                raw_list.append(data_row)

                # ONLY add to Filtered log if it meets the scientific threshold
                if idx < 40:
                    filtered_list.append(data_row)

        except Exception as e:
            # Print newline before warning to preserve the dynamic status bar
            print(f"\n[WARNING] Mathematical evaluation failure on {filename}. Reason: {e}")
            continue

    # Drop down a line after the dynamic progress bar finishes
    print() 

    # 4. DATA EXPORT & REPORTING
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
    raw_csv = os.path.join(raw_out, f"Phase1B_Raw_Full_{timestamp}.csv")
    filt_csv = os.path.join(filt_out, f"Phase1B_Filtered_Stable_{timestamp}.csv")

    if raw_list:
        # Save Raw CSV
        with open(raw_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=raw_list[0].keys())
            writer.writeheader()
            writer.writerows(raw_list)
        
        # Save Filtered CSV
        with open(filt_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=filtered_list[0].keys())
            writer.writeheader()
            writer.writerows(filtered_list)

        # Post-Execution Statistics
        total_time = format_time(time.time() - start_time)
        survivor_rate = (len(filtered_list) / len(raw_list)) * 100 if raw_list else 0

        print("\n" + "="*80)
        print(f"{'STABILITY ANALYSIS COMPLETE':^80}")
        print("="*80)
        print(f"[SUCCESS] Total Sequences Analyzed : {len(raw_list)}")
        print(f"[SUCCESS] Stable Candidates Passed : {len(filtered_list)} ({survivor_rate:.1f}% Survivor Rate)")
        print(f"[SUCCESS] Total Execution Time     : {total_time}")
        print("-" * 80)
        print(f"[INFO] Raw Matrix Log    : {os.path.basename(raw_csv)}")
        print(f"[INFO] Filtered Matrix   : {os.path.basename(filt_csv)}")
        print("="*80 + "\n")

if __name__ == "__main__":
    run_step1b_final_dual_output()