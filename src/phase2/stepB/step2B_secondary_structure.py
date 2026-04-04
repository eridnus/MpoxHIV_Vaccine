import os
import csv
import sys
import time
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# =============================================================================
# ACADEMIC FORMATTING UTILITIES
# =============================================================================

def print_banner(text):
    print("\n" + "="*110)
    print(f"{text:^110}")
    print("="*110)

def format_time(seconds):
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

# =============================================================================
# CORE SECONDARY STRUCTURE ENGINE
# =============================================================================

def run_step2b_secondary_structure():
    start_time = time.time()

    # 1. DYNAMIC PATHING
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    
    # Input: The Filtered Optimal Candidates from Step 2A
    # Note: We specifically want to process the Rank 1 winner from your last run
    input_csv = os.path.join(project_root, "Step_Outputs", "Phase2", "StepA", "Filtered")
    variant_folder = os.path.join(project_root, "Step_Outputs", "Phase1G", "Variants")
    output_base = os.path.join(project_root, "Step_Outputs", "Phase2", "StepB")
    
    os.makedirs(output_base, exist_ok=True)

    # Find the latest filtered CSV to identify the Rank 1 winner
    csv_files = sorted([f for f in os.listdir(input_csv) if f.endswith(".csv")])
    if not csv_files:
        print("[ERROR] No filtered candidates found. Run Step 2A first.")
        return
    
    with open(os.path.join(input_csv, csv_files[-1]), 'r') as f:
        reader = list(csv.DictReader(f))
        winner_name = reader[0]['Variant'] # Rank 1 Candidate

    print_banner("PHASE 2 STEP B: SECONDARY STRUCTURE PROPENSITY ANALYSIS")
    print(f"[INFO] Target Variant    : {winner_name} (Rank 1)")
    print(f"[INFO] Methodology Reference: PsiPred / RaptorX Standards")
    print("-" * 110)

    # 2. LOAD WINNING SEQUENCE
    with open(os.path.join(variant_folder, winner_name), 'r') as f:
        lines = f.readlines()
        sequence = "".join([l.strip() for l in lines[1:]]).upper()

    # 3. STRUCTURAL PROPENSITY CALCULATIONS (Secondary Structure Fraction)
    analysis = ProteinAnalysis(sequence)
    helix, turn, sheet = analysis.secondary_structure_fraction()
    
    # Calculate Coil percentage (The remainder)
    coil = 1.0 - (helix + sheet + turn)

    # 4. TERMINAL DASHBOARD
    print(f"{'STRUCTURAL ELEMENT':<30} | {'FRACTION':<15} | {'PERCENTAGE':<15}")
    print("-" * 110)
    print(f"{'Alpha-Helix (H)':<30} | {helix:>15.4f} | {helix*100:>14.2f}%")
    print(f"{'Beta-Sheet (E)':<30} | {sheet:>15.4f} | {sheet*100:>14.2f}%")
    print(f"{'Turns (T)':<30} | {turn:>15.4f} | {turn*100:>14.2f}%")
    print(f"{'Random Coils (C)':<30} | {coil:>15.4f} | {coil*100:>14.2f}%")
    print("-" * 110)

    # 5. EXPORT FOR PSI-PRED / RAPTOR-X
    # These servers require a clean FASTA with a specific header
    submission_fasta = os.path.join(output_base, f"Step2B_Submission_{winner_name}")
    with open(submission_fasta, 'w') as f:
        f.write(f">Submission_PsiPred_RaptorX_{winner_name}\n{sequence}\n")

    # 6. LOGGING RESULTS
    report_path = os.path.join(output_base, "Step2B_Secondary_Structure_Report.csv")
    with open(report_path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Element", "Fraction", "Percentage"])
        writer.writerow(["Alpha-Helix", helix, f"{helix*100:.2f}%"])
        writer.writerow(["Beta-Sheet", sheet, f"{sheet*100:.2f}%"])
        writer.writerow(["Turns", turn, f"{turn*100:.2f}%"])
        writer.writerow(["Random Coils", coil, f"{coil*100:.2f}%"])

    total_time = format_time(time.time() - start_time)
    print_banner("STEP 2B COMPLETE")
    print(f"[SUCCESS] Secondary Propensity Calculated.")
    print(f"[SUCCESS] Execution Time : {total_time}")
    print(f"[INFO] Submission FASTA  : {os.path.relpath(submission_fasta, project_root)}")
    print(f"[INFO] Analysis Report   : {os.path.relpath(report_path, project_root)}")
    print("="*110 + "\n")

if __name__ == "__main__":
    run_step2b_secondary_structure()