import csv
import os

def check_input():
    # Adjust this path to your actual input csv file
    path = "Step_Outputs/Phase1F/Filtered/Phase1F_Elite_Vaccine_Candidates_2026-04-04_1809.csv"
    with open(path, 'r') as f:
        reader = csv.DictReader(f)
        peptides = [row['Peptide'] for row in reader]
    
    unique_peptides = set(peptides)
    print(f"\n[DEBUG] Total epitopes in CSV: {len(peptides)}")
    print(f"[DEBUG] Unique epitopes found: {len(unique_peptides)}")
    if len(peptides) == len(unique_peptides):
        print("[OK] All epitopes are unique.")
    else:
        print("[!!!] WARNING: You have duplicates in your input data!")

if __name__ == "__main__":
    check_input()