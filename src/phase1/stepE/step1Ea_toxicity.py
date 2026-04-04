import os
import sys
import time
import csv
from datetime import datetime

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def run_step1ea_toxicity():
    start_time = time.time()
    
    # DYNAMIC PATHING
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    # Input: Latest benchmark from 1Dc
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1D", "Phase1Dc", "Filtered_Benchmarks")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase1E", "Phase1Ea")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)

    print("\n" + "="*80)
    print(f"{'PHASE 1Ea: TOXICITY SCREENING':^80}")
    print("="*80)

    # Locate latest file
    files = sorted([f for f in os.listdir(input_folder) if "Min_50pct" in f])
    if not files:
        print("[ERROR] No 50% Benchmark file found in 1Dc/Filtered_Benchmarks.")
        return
    input_csv = os.path.join(input_folder, files[-1])
    
    print(f"[INFO] Processing input: {os.path.basename(input_csv)}")

    raw_data, filtered_data = [], []
    fieldnames = ["Target", "Variant", "Peptide", "Toxicity_Score", "Cysteine_Count", "Toxicity_Status"]

    with open(input_csv, 'r') as f:
        reader = list(csv.DictReader(f))
        total_rows = len(reader)
        
        for i, row in enumerate(reader):
            pep = row['Peptide'].upper()
            c_count = pep.count('C')
            hydro_ratio = sum(1 for aa in pep if aa in "VILFMWAY") / len(pep)
            is_toxic = (c_count > 2) or (hydro_ratio > 0.8)
            
            # Progress Tracking
            elapsed = format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{total_rows:03d} | Toxicity Check | Elapsed: {elapsed}")
            sys.stdout.flush()

            clean_row = {
                "Target": row.get('Target', 'N/A'),
                "Variant": row.get('Variant', 'N/A'),
                "Peptide": pep,
                "Toxicity_Score": round(hydro_ratio, 2),
                "Cysteine_Count": c_count,
                "Toxicity_Status": "TOXIC" if is_toxic else "NON-TOXIC"
            }
            
            raw_data.append(clean_row)
            if not is_toxic:
                filtered_data.append(clean_row)

    print(f"\n[INFO] Screening complete.")

    # EXPORT
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    with open(os.path.join(raw_dir, f"Phase1Ea_Raw_{ts}.csv"), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(raw_data)
        
    with open(os.path.join(filt_dir, f"Phase1Ea_Filtered_{ts}.csv"), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(filtered_data)

    print(f"[SUCCESS] {len(filtered_data)} non-toxic candidates identified.")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_step1ea_toxicity()