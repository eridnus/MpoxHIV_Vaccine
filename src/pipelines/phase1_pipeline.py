import os
import subprocess
import sys
import time
from datetime import datetime

# =============================================================================
# MASTER PIPELINE CONFIGURATION
# =============================================================================

def get_project_paths():
    """
    Resolves project root dynamically. 
    Location: src/pipelines/phase1_pipeline.py
    Root: ../..
    """
    script_dir = os.path.dirname(os.path.abspath(__file__))
    src_dir = os.path.dirname(script_dir) # Moves from pipelines/ to src/
    project_root = os.path.dirname(src_dir) # Moves from src/ to root/
    
    return {
        "root": project_root,
        "phase1_dir": os.path.join(src_dir, "phase1")
    }

def print_banner(text):
    print("\n" + "="*80)
    print(f"{text:^80}")
    print("="*80)

# =============================================================================
# ORCHESTRATION LOGIC
# =============================================================================

def run_pipeline():
    paths = get_project_paths()
    start_time = time.time()
    
    # Ordered Pipeline Steps: (Label, Folder, Script)
    pipeline_steps = [
        ("PHASE 1A: RETRIEVAL", "stepA", "step1A_retrieval.py"),
        ("PHASE 1B: STABILITY", "stepB", "step1B_stability.py"),
        ("PHASE 1C: ANTIGENICITY", "stepC", "step1C_antigenicity.py"),
        ("PHASE 1Da: EPITOPE ID", "stepD", "step1Da_epitopes_identification.py"),
        ("PHASE 1Db: IEDB FILTRATION", "stepD", "step1Db_threshold_filtration.py"),
        ("PHASE 1Dc: CONSERVANCY", "stepD", "step1Dc_conservancy_benchmark.py"),
        ("PHASE 1Ea: TOXICITY", "stepE", "step1Ea_toxicity.py"),
        ("PHASE 1Eb: ALLERGENICITY", "stepE", "step1Eb_allergenicity.py"),
        ("PHASE 1F: CONVERGENCE", "stepF", "step1F_convergence.py"),
        ("PHASE 1G: ASSEMBLY", "stepG", "step1G_construction.py")
    ]

    print_banner("MPOX-HIV VACCINE DESIGN: PHASE 1 PIPELINE")
    print(f"[INFO] Root Directory : {paths['root']}")
    print(f"[INFO] Pipeline Start  : {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("-" * 80)

    # Execute Sequence
    for label, folder, script in pipeline_steps:
        script_path = os.path.join(paths['phase1_dir'], folder, script)
        
        if not os.path.exists(script_path):
            print(f"\n[FATAL] Script missing: {script_path}")
            sys.exit(1)

        print(f"\n[EXEC] Triggering -> {label}")
        
        try:
            # Execute and stream output
            process = subprocess.Popen(
                [sys.executable, script_path],
                stdout=sys.stdout,
                stderr=sys.stderr,
                bufsize=1
            )
            process.wait()
            
            if process.returncode != 0:
                print(f"\n[FATAL] {label} terminated with code {process.returncode}.")
                sys.exit(1)
                
        except Exception as e:
            print(f"\n[FATAL] Execution crash in {label}: {e}")
            sys.exit(1)

    # FINAL SUMMARY
    duration = time.time() - start_time
    mins, secs = divmod(int(duration), 60)
    
    print_banner("PHASE 1 COMPLETE: DATA INTEGRITY VERIFIED")
    print(f"[SUCCESS] Total Pipeline Execution Time: {mins:02d}m:{secs:02d}s")
    print(f"[INFO] Final Construct saved in: {os.path.join(paths['root'], 'Step_Outputs', 'Phase1G')}")
    print("="*80 + "\n")

if __name__ == "__main__":
    run_pipeline()