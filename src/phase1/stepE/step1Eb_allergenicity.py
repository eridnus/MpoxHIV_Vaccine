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

def run_step1eb_allergenicity():
    start_time = time.time()
    
    # DYNAMIC PATHING
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1E", "Phase1Ea", "Filtered")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase1E", "Phase1Eb")
    raw_dir = os.path.join(output_base, "Raw")
    filt_dir = os.path.join(output_base, "Filtered")
    os.makedirs(raw_dir, exist_ok=True)
    os.makedirs(filt_dir, exist_ok=True)

    print("\n" + "="*80)
    print(f"{'PHASE 1Eb: ALLERGENICITY SCREENING':^80}")
    print("="*80)

    files = sorted([f for f in os.listdir(input_folder) if f.endswith(".csv")])
    if not files:
        print("[ERROR] No input CSV found in 1Ea/Filtered.")
        return
    input_csv = os.path.join(input_folder, files[-1])
    
    print(f"[INFO] Processing input: {os.path.basename(input_csv)}")

    raw_data, filtered_data = [], []
    fieldnames = ["Target", "Variant", "Peptide", "QN_Ratio", "Surface_Charge", "Allergen_Status"]

    with open(input_csv, 'r') as f:
        reader = list(csv.DictReader(f))
        total_rows = len(reader)
        
        for i, row in enumerate(reader):
            pep = row['Peptide'].upper()
            qn_ratio = (pep.count('Q') + pep.count('N')) / len(pep)
            charge_count = sum(1 for aa in pep if aa in "DEHK")
            is_allergen = (qn_ratio > 0.3) or (charge_count > 4)
            
            # Progress Tracking
            elapsed = format_time(time.time() - start_time)
            sys.stdout.write(f"\r[ PROCESS ] {i+1:03d}/{total_rows:03d} | Allergen Check | Elapsed: {elapsed}")
            sys.stdout.flush()

            clean_row = {
                "Target": row['Target'],
                "Variant": row['Variant'],
                "Peptide": pep,
                "QN_Ratio": round(qn_ratio, 2),
                "Surface_Charge": charge_count,
                "Allergen_Status": "ALLERGEN" if is_allergen else "NON-ALLERGEN"
            }
            
            raw_data.append(clean_row)
            if not is_allergen:
                filtered_data.append(clean_row)

    print(f"\n[INFO] Screening complete.")

    # EXPORT
    ts = datetime.now().strftime("%Y-%m-%d_%H%M")
    with open(os.path.join(raw_dir, f"Phase1Eb_Raw_{ts}.csv"), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(raw_data)
        
    with open(os.path.join(filt_dir, f"Phase1Eb_Filtered_{ts}.csv"), 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(filtered_data)

    print(f"[SUCCESS] {len(filtered_data)} non-allergenic candidates identified.")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_step1eb_allergenicity()