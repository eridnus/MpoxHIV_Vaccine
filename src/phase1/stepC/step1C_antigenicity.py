import os
import sys
import time
import csv
import re
from datetime import datetime

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def format_time(seconds):
    """Formats time in seconds to a readable MM:SS format."""
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def calculate_antigenicity(sequence):
    """
    Evaluates sequence immunogenicity using the Kolaskar & Tongaonkar scale.
    Threshold for standard antigenic potential is >= 1.00.
    """
    kt_scale = {
        'A': 1.064, 'C': 1.412, 'D': 0.866, 'E': 0.851, 'F': 1.091,
        'G': 0.874, 'H': 1.105, 'I': 1.152, 'K': 0.930, 'L': 1.250,
        'M': 1.126, 'N': 0.851, 'P': 1.064, 'Q': 1.010, 'R': 0.873,
        'S': 1.012, 'T': 0.909, 'V': 1.187, 'W': 1.085, 'Y': 1.255
    }
    
    clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', sequence.upper())
    if not clean_seq: return 0
    
    values =[kt_scale.get(aa, 1.0) for aa in clean_seq]
    return sum(values) / len(values)

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def run_step1c_unified_antigenicity():
    start_time = time.time()

    # 1. DYNAMIC DIRECTORY RESOLUTION
    # Current location: [ROOT]/src/phase1/stepC/step1C_antigenicity.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Traverse up 3 levels to reach the project root
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1A")
    output_path = os.path.join(project_root, "Step_Outputs", "Phase1C")
    
    # Unified Folder Structure
    raw_out = os.path.join(output_path, "Raw_Antigenicity")
    filt_out = os.path.join(output_path, "Filtered_Antigenicity")
    os.makedirs(raw_out, exist_ok=True)
    os.makedirs(filt_out, exist_ok=True)

    if not os.path.exists(input_folder):
        print("\n[ERROR] Phase 1A output directory not found. Please execute Step 1A first.")
        return

    fasta_files = sorted([f for f in os.listdir(input_folder) if f.endswith(".fasta")])
    total_files = len(fasta_files)

    if total_files == 0:
        print("\n[ERROR] No FASTA files found in Phase 1A directory.")
        return

    # 2. EXPERIMENTAL PREVIEW & HEADER
    print("\n" + "="*80)
    print(f"{'PHASE 1C: KOLASKAR & TONGAONKAR ANTIGENICITY SCREENING':^80}")
    print("="*80)
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Input Directory     : {input_folder}")
    print(f"[INFO] Target Files        : {total_files} FASTA sequences")
    print(f"[INFO] Antigenic Threshold : Score >= 1.00")
    print("-" * 80)

    # 3. MATHEMATICAL EVALUATION MATRIX
    raw_results = []
    filtered_results = []
    
    sys.stdout.write("[PROCESS] Initiating Antigenicity algorithm...\n\n")

    for i, filename in enumerate(fasta_files):
        file_path = os.path.join(input_folder, filename)
        variant_id = filename.replace(".fasta", "")
        target_group = filename.split('_Var')[0]
        
        try:
            with open(file_path, "r") as f:
                sequence = "".join([line.strip() for line in f if not line.startswith(">")])

            score = calculate_antigenicity(sequence)
            status = "ANTIGENIC" if score >= 1.0 else "LOW-ANTIGEN"
            
            # Original Qualitative Logic
            if score > 1.05: 
                note = "High probability of B-cell epitopes; strong candidate."
            elif score >= 1.0: 
                note = "Standard antigenic profile; acceptable for vaccine use."
            else: 
                note = "Lower immunogenic potential."

            data_row = {
                "Variant_ID": variant_id,
                "Target": target_group,
                "Antigenicity_Score": round(score, 4),
                "Status": status,
                "Notes": note
            }

            # DUAL OUTPUT LOGIC
            raw_results.append(data_row)
            if score >= 1.0:
                filtered_results.append(data_row)

            # Dynamic terminal preview string
            elapsed = format_time(time.time() - start_time)
            status_str = f"[ EVAL ] Record {i+1:03d}/{total_files:03d} | Score: {score:.3f} ({status}) | Elapsed: {elapsed}"
            sys.stdout.write(f"\r{status_str:<80}")
            sys.stdout.flush()

        except Exception as e:
            print(f"\n[WARNING] Evaluation failure on {filename}. Reason: {e}")
            continue

    # Drop down a line after the dynamic progress bar finishes
    print()

    # 4. DATA EXPORT & REPORTING
    timestamp = datetime.now().strftime("%Y-%m-%d_%H%M")
    raw_csv = os.path.join(raw_out, f"Phase1C_Raw_Antigen_{timestamp}.csv")
    filt_csv = os.path.join(filt_out, f"Phase1C_Filtered_Antigenic_{timestamp}.csv")

    if raw_results:
        # Save RAW
        with open(raw_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=raw_results[0].keys())
            writer.writeheader()
            writer.writerows(raw_results)
        
        # Save FILTERED
        with open(filt_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=filtered_results[0].keys())
            writer.writeheader()
            writer.writerows(filtered_results)

        # Post-Execution Statistics
        total_time = format_time(time.time() - start_time)
        survivor_rate = (len(filtered_results) / len(raw_results)) * 100 if raw_results else 0

        print("\n" + "="*80)
        print(f"{'ANTIGENICITY SCREENING COMPLETE':^80}")
        print("="*80)
        print(f"[SUCCESS] Total Sequences Analyzed : {len(raw_results)}")
        print(f"[SUCCESS] Antigenic Candidates     : {len(filtered_results)} ({survivor_rate:.1f}% Survivor Rate)")
        print(f"[SUCCESS] Total Execution Time     : {total_time}")
        print("-" * 80)
        print(f"[INFO] Raw Matrix Log    : {os.path.basename(raw_csv)}")
        print(f"[INFO] Filtered Matrix   : {os.path.basename(filt_csv)}")
        print("="*80 + "\n")

if __name__ == "__main__":
    run_step1c_unified_antigenicity()