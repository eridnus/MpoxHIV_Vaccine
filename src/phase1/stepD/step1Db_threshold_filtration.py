import os, sys, time, re, csv, requests
from datetime import datetime

def print_banner(text): print(f"\n{'='*80}\n{text:^80}\n{'='*80}")

def get_gravy(pep):
    """Calculates GRAVY (Grand Average of Hydropathy)."""
    hydro = {'A': 1.8, 'L': 3.8, 'I': 4.5, 'V': 4.2, 'F': 2.8, 'M': 1.9, 'C': 2.5, 'G': -0.4, 
             'P': -1.6, 'S': -0.8, 'T': -0.7, 'W': -0.9, 'Y': -1.3, 'N': -3.5, 'Q': -3.5, 
             'D': -3.5, 'E': -3.5, 'K': -3.9, 'R': -4.5, 'H': -3.2}
    return sum(hydro.get(aa, 0) for aa in pep) / len(pep)

def run_step1db_optimized():
    # 1. Path Resolution
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    input_folder = os.path.join(project_root, "Step_Outputs", "Phase1A")
    output_folder = os.path.join(project_root, "Step_Outputs", "Phase1D", "Phase1Db")
    os.makedirs(output_folder, exist_ok=True)
    
    print_banner("PHASE 1Db: IEDB AFFINITY & SOLUBILITY FILTER")
    
    fasta_files = sorted([f for f in os.listdir(input_folder) if f.endswith(".fasta")])
    all_results = []
    skipped = 0

    # 2. Main Processing
    MHCI_URL = "https://tools-cluster-interface.iedb.org/tools_api/mhci/"
    
    for i, f_name in enumerate(fasta_files):
        target = f_name.split('_Var')[0]
        sys.stdout.write(f"\r[PROCESS] {i+1:02d}/{len(fasta_files):02d} | Target: {target:<12} | Filtered: {skipped}")
        
        with open(os.path.join(input_folder, f_name), "r") as f:
            seq = "".join([l.strip() for l in f.readlines()[1:]])
            clean_seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', seq.upper())

        if len(clean_seq) < 9: continue

        # IEDB Request
        payload = {'method': 'recommended', 'sequence_text': clean_seq, 'allele': 'HLA-A*02:01', 'length': '9'}
        try:
            response = requests.post(MHCI_URL, data=payload, timeout=90)
            if response.status_code == 200:
                lines = response.text.strip().split('\n')
                if len(lines) > 1:
                    header = lines[0].split('\t')
                    # Find indices dynamically
                    r_idx = header.index('percentile_rank')
                    p_idx = header.index('peptide')
                    
                    for line in lines[1:]:
                        cols = line.split('\t')
                        pep = cols[p_idx]
                        if float(cols[r_idx]) <= 1.0: 
                            # SOLUBILITY FILTER: Only accept soluble binders
                            if get_gravy(pep) < 0.2: 
                                all_results.append({"Target": target, "Peptide": pep, "Rank": cols[r_idx]})
                            else:
                                skipped += 1
            time.sleep(1.0) # Respect rate limits
        except Exception as e:
            continue

    # 3. Export
    out_file = os.path.join(output_folder, f"Phase1Db_Elite_{datetime.now().strftime('%Y%m%d')}.csv")
    with open(out_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Target", "Peptide", "Rank"])
        writer.writeheader()
        writer.writerows(all_results)
    
    print_banner("FILTRATION COMPLETE")
    print(f"[SUCCESS] {len(all_results)} high-affinity soluble binders saved.")
    print(f"[INFO] Log file: {os.path.basename(out_file)}\n")

if __name__ == "__main__": run_step1db_optimized()