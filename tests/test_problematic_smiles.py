"""
Tests that Gypsum-DL can gracefully handle SMILES strings that cause errors
during 3D coordinate generation, without crashing.
"""
import os
import shutil
from gypsum_dl import utils
from gypsum_dl.start import prepare_molecules

def test_problematic_smiles(test_dir):
    """
    Tests that Gypsum-DL can process known problematic SMILES without crashing.

    The primary goal is to ensure that `prepare_molecules` completes successfully
    and does not raise an unhandled exception, which would indicate a crash.
    The test passes if the function call finishes, regardless of whether
    individual molecules succeed or fail.
    """
    output_folder = os.path.join(test_dir, "tmp/problematic_smiles_test")
    
    # 1. Define problematic SMILES and set up input/output directories
    problematic_smiles_data = {
        "problematic_nitro": "O=C1CCCCC[N@@H+]1C(=O)Nc1cccc([N+](=O)[O-])c1",
        "problematic_zinc_complex": "CC(=O)[O-].CC(=O)[O-][Zn+2]1234S=C(N)N[N+]1=C(C)c1cccc(C(=[N+]2NC(=S3)N)C)[n+]41 601849",
        "invalid_smiles": "moosedogfacecat",
    }
    
    # Delete test output directory if it exists, then create it.
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)

    # Create the input smiles file
    input_smi_path = os.path.join(output_folder, "input.smi")
    with open(input_smi_path, "w") as f:
        for name, smi in problematic_smiles_data.items():
            f.write(f"{smi}\t{name}\n")

    # 2. Set up Gypsum-DL parameters
    params = {
        "source": input_smi_path,
        "output_folder": output_folder,
        "job_manager": "serial",
        "max_variants_per_compound": 1,
        "thoroughness": 1,
        "separate_output_files": True,
    }

    # 3. Run molecule preparation. This test passes if this call completes
    #    without raising an unhandled exception (i.e., it does not crash).
    try:
        prepare_molecules(params)
        # If the function completes, the test has passed.
        utils.log("")
        utils.log("TEST `test_problematic_smiles` PASSED")
        utils.log("=====================================")
        utils.log("`prepare_molecules` completed without crashing on problematic inputs.")
        assert True
    except Exception as e:
        # If any unhandled exception occurs, the program has crashed.
        assert False, f"CRASH DETECTED: `prepare_molecules` raised an unexpected exception: {e}"