import os, random, csv, hashlib

def print_banner(text): print(f"\n{'='*80}\n{text:^80}\n{'='*80}")

def run_step1g_final_anchor():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1F", "Filtered")
    variant_dir = os.path.join(project_root, "Step_Outputs", "Phase1G", "Variants")
    os.makedirs(variant_dir, exist_ok=True)

    print_banner("PHASE 1G: STABILITY-ANCHORED ASSEMBLY (EAAAK)")
    
    files = sorted([f for f in os.listdir(input_folder) if f.endswith(".csv")])
    with open(os.path.join(input_folder, files[-1]), 'r') as f:
        # Get unique epitopes
        epitopes = list(set([row['Peptide'] for row in csv.DictReader(f)]))

    # THE STABILITY ANCHOR CONFIGURATION
    sol_tag = "MHHHHHH"
    adjuvant = "GIINTLQKYYCRVRGGRCAVLSCLPKEEQIGKCSTRGRKCCRRKK"
    linker = "EAAAK" # The Stability Anchor

    for v in range(1, 6):
        sample = random.sample(epitopes, 15)
        # Sequence: Tag -> Adjuvant -> Linker -> [Epitope-Linker-Epitope]
        # We use EAAAK between EVERY element to ensure max stability
        construct = f"{sol_tag}{adjuvant}{linker}{linker.join(sample)}"
        h = hashlib.md5(construct.encode()).hexdigest()[:8]
        
        out_file = os.path.join(variant_dir, f"Variant_{v}.fasta")
        with open(out_file, 'w') as f:
            f.write(f">MpoxHIV_Variant_{v}|Hash:{h}\n{construct}\n")
        print(f"[ASSEMBLY] Variant_{v} | Length: {len(construct)} AA | Hash: {h}")
        
    print_banner("ASSEMBLY COMPLETE")

if __name__ == "__main__": run_step1g_final_anchor()