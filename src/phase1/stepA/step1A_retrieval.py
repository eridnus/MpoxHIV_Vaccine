import os
import sys
import time
from datetime import datetime
from Bio import Entrez

# =============================================================================
# EXPERIMENTAL CONFIGURATION
# =============================================================================

# NCBI Identification (Required for API compliance)
Entrez.email = "enzoleonor.3309@gmail.com" 

# Target viral proteomes for Chimeric Vaccine construct
TARGET_CLUSTERS = {
    "Mpox_L1R": "L1R monkeypox",
    "Mpox_B5R": "B5R monkeypox",
    "Mpox_A35R": "A35R monkeypox",
    "HIV_gp120": "envelope glycoprotein gp120 hiv-1",
    "HIV_gp41": "envelope glycoprotein gp41 hiv-1",
    "HIV_p24": "gag p24 hiv-1",
    "HIV_p17": "gag p17 hiv-1"
}

# Sampling parameters
VARIANTS_PER_TARGET = 30 
TOTAL_PROJECTED = len(TARGET_CLUSTERS) * VARIANTS_PER_TARGET

# =============================================================================
# CORE PROCESSING FUNCTION
# =============================================================================

def format_time(seconds):
    """Formats time in seconds to a readable MM:SS format."""
    mins, secs = divmod(int(seconds), 60)
    return f"{mins:02d}m:{secs:02d}s"

def run_high_density_retrieval():
    start_time = time.time()
    
    # 1. DYNAMIC DIRECTORY RESOLUTION
    # Current location: [ROOT]/src/phase1/stepA/step1A_retrieval.py
    script_dir = os.path.dirname(os.path.abspath(__file__))
    # Traverse up 3 levels to reach the project root
    project_root = os.path.abspath(os.path.join(script_dir, "..", "..", ".."))
    phase1a_path = os.path.join(project_root, "Step_Outputs", "Phase1A")
    
    os.makedirs(phase1a_path, exist_ok=True)

    # 2. EXPERIMENTAL PREVIEW & HEADER
    print("\n" + "="*80)
    print(f"{'PHASE 1A: VIRAL PROTEOME SEQUENCE RETRIEVAL (NCBI ENTREZ)':^80}")
    print("="*80)
    print(f"[INFO] Initialization Time : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"[INFO] Target Antigens     : {len(TARGET_CLUSTERS)} defined")
    print(f"[INFO] Variants per Target : {VARIANTS_PER_TARGET}")
    print(f"[INFO] Projected Yield     : {TOTAL_PROJECTED} FASTA sequences")
    print(f"[INFO] Output Directory    : {phase1a_path}")
    print("-" * 80)

    # 3. WORKSPACE PREPARATION (Purging old data)
    sys.stdout.write("[PROCESS] Purging prior sequence data to ensure experimental integrity...")
    sys.stdout.flush()
    purged_count = 0
    for f in os.listdir(phase1a_path):
        if f.endswith(".fasta"):
            os.remove(os.path.join(phase1a_path, f))
            purged_count += 1
    print(f" Done. ({purged_count} files removed)")
    print("-" * 80)

    # 4. SEQUENCE ACQUISITION PROTOCOL
    successful_downloads = 0

    for label, query in TARGET_CLUSTERS.items():
        print(f"\n[INFO] Establishing NCBI connection for target: {label}")
        
        try:
            # 4a. Query execution
            search_handle = Entrez.esearch(db="protein", term=query, retmax=VARIANTS_PER_TARGET)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            id_list = search_results["IdList"]
            retrieved_count = len(id_list)
            
            if retrieved_count == 0:
                print(f"[WARNING] No records found for query: '{query}'")
                continue

            # 4b. Data fetching with dynamic terminal output
            for i, ncbi_id in enumerate(id_list):
                # Calculate metrics
                elapsed = format_time(time.time() - start_time)
                
                # Dynamic terminal preview string (padded with spaces to overwrite previous lines cleanly)
                status_str = f"[ FETCH ] {label} | Record {i+1:02d}/{retrieved_count:02d} | ID: {ncbi_id:<12} | Elapsed: {elapsed}"
                sys.stdout.write(f"\r{status_str:<80}")
                sys.stdout.flush()
                
                # Download payload
                fetch_handle = Entrez.efetch(db="protein", id=ncbi_id, rettype="fasta", retmode="text")
                fasta_data = fetch_handle.read()
                fetch_handle.close()
                
                # Write to disk
                file_name = f"{label}_Var_{i+1:02d}_{ncbi_id}.fasta"
                with open(os.path.join(phase1a_path, file_name), "w") as f:
                    f.write(fasta_data)
                
                successful_downloads += 1
                
                # Regulated delay to comply with NCBI API restrictions (Max 3 req/sec)
                time.sleep(1.0) 
            
            # Print newline after completing the target batch
            print() 

        except Exception as e:
            print(f"\n[ERROR] Protocol failure during {label} acquisition. Reason: {e}")

    # 5. POST-EXECUTION REPORT
    total_time = format_time(time.time() - start_time)
    print("\n" + "="*80)
    print(f"{'ACQUISITION PROTOCOL COMPLETE':^80}")
    print("="*80)
    print(f"[SUCCESS] Total Sequences Retrieved : {successful_downloads}/{TOTAL_PROJECTED}")
    print(f"[SUCCESS] Total Execution Time      : {total_time}")
    print(f"[INFO] Data formatting complete. Proceed to Phase 1B for stability thresholding.")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_high_density_retrieval()