'''
PyMOL Protein Folding Plugin

By Colby T. Ford, Ph.D.
License: GPLv3
'''

from __future__ import absolute_import
from __future__ import print_function
import os, tempfile, random, string, sys, subprocess, json, requests
from pathlib import Path

# Add a check for the new PyYAML dependency
try:
    import yaml
except ModuleNotFoundError:
    raise ImportError("The 'PyYAML' library is required for the new Boltz functionality. Please install it in your PyMOL's Python environment (e.g., 'pip install pyyaml').")

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PyMOLfold', run_plugin_gui)

## Global reference to avoid garbage collection of our dialog
dialog = None

def run_plugin_gui():
    '''
    Open the custom dialog
    '''
    global dialog
    if dialog is None:
        dialog = make_dialog()
    dialog.show()

#########################
## Folding Functions
#########################

## Helper Functions
def fix_multichain_pdb_str(pdb_string:str = "", split_res:str = "UNK") -> str:
    """
    Fixes chain identifiers in a PDB string for multiple chains.
    """
    chain_ids = list(string.ascii_uppercase + string.ascii_lowercase)
    fixed_lines = []
    chain_idx = 0
    for line in pdb_string.splitlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            residue = line[17:20]
            line = line[:21] + chain_ids[chain_idx] + line[22:]
            if residue == split_res:
                chain_idx += 1
                continue
            fixed_lines.append(line)
    return "\n".join(fixed_lines)

## ESM Folding
def fold_esm(model_name:str, aa_sequence:str, temperature:float=0.7, num_steps:int=8, token:str=""):
    """
    Protein folding using ESM models
    """
    try:
        from esm.sdk import client
        from esm.sdk.api import ESMProtein, GenerationConfig
    except ModuleNotFoundError as e:
        raise Exception(f"esm module not found: {str(e)}")

    try:
        model = client(model=model_name, url="https://forge.evolutionaryscale.ai", token=token)
    except Exception as e:
        raise Exception(f"Error getting ESM model with token: {str(e)}")

    structure_prediction_config = GenerationConfig(track="structure", num_steps=num_steps, temperature=temperature)
    structure_prediction_prompt = ESMProtein(sequence=aa_sequence)
    structure_prediction = model.generate(structure_prediction_prompt, structure_prediction_config)
    structure_prediction_chain = structure_prediction.to_protein_chain()
    pdb_string = structure_prediction_chain.to_pdb_string()

    if "multimer" in model_name:
        pdb_string = fix_multichain_pdb_str(pdb_string, split_res="UNK")

    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb:
        temp_pdb.write(pdb_string.encode())
        temp_pdb_path = temp_pdb.name

    return temp_pdb_path, None

## Chai Folding
def fold_chai(aa_sequence:str, ligand:str=None, ligand_type:str="smiles", num_trunk_recycles:int=3, num_diffn_timesteps:int=200, seed:int=1337):
    """
    Protein folding using Chai models
    """
    try:
        from chai_lab.chai1 import run_inference
        import torch
    except ModuleNotFoundError as e:
        raise Exception(f"chai_lab module not found: {str(e)}")
    
    fasta_content = f">protein|name=chain_A\n{aa_sequence}\n"
    if ligand and ligand_type:
        fasta_content += f">ligand|name=chain_B\n{ligand}\n"

    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
        temp_fasta.write(fasta_content.encode())
        temp_fasta_path = temp_fasta.name

    output_dir = tempfile.mkdtemp()
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    candidates = run_inference(fasta_file=Path(temp_fasta_path), output_dir=Path(output_dir), num_trunk_recycles=num_trunk_recycles, num_diffn_timesteps=num_diffn_timesteps, seed=seed, device=device, use_esm_embeddings=True)

    if not candidates.cif_paths:
        raise ValueError("No structure files were generated")

    return candidates.cif_paths[0], None

## Boltz Folding
def fold_boltz(
    protein_seqs_str: str, protein_ids_str: str,
    ligand_seqs_str: str, ligand_ids_str: str, ligand_type: str,
    na_seqs_str: str, na_ids_str: str, na_type: str,
    use_msa_server: bool, recycling_steps: int, sampling_steps: int, is_cyclic: bool, 
    predict_affinity: bool, binder_chain_id: str, force_cpu: bool, no_kernels: bool
    ):
    """
    Generates YAML from structured UI input and runs Boltz prediction.
    """
    try:
        import torch
    except ModuleNotFoundError as e:
        raise Exception(f"Could not import required module: {str(e)}")

    # Add version key to match Boltz examples
    yaml_data = {'version': 1, 'sequences': []}
    
    # Process Proteins
    protein_seqs = [s.strip() for s in protein_seqs_str.strip().splitlines() if s.strip()]
    protein_ids = [i.strip() for i in protein_ids_str.strip().split(',') if i.strip()]
    if protein_seqs and len(protein_seqs) != len(protein_ids):
        raise ValueError(f"Mismatch between number of protein sequences ({len(protein_seqs)}) and chain IDs ({len(protein_ids)}).")
    for i, seq in enumerate(protein_seqs):
        entity_content = {'id': protein_ids[i], 'sequence': seq}
        if is_cyclic:
            entity_content['cyclic'] = True
        if use_msa_server:
            entity_content['msa'] = 'empty'
        yaml_data['sequences'].append({'protein': entity_content})

    # Process Ligands
    ligand_seqs = [s.strip() for s in ligand_seqs_str.strip().splitlines() if s.strip()]
    ligand_ids = [i.strip() for i in ligand_ids_str.strip().split(',') if i.strip()]
    if ligand_seqs and len(ligand_seqs) != len(ligand_ids):
        raise ValueError(f"Mismatch between number of ligands ({len(ligand_seqs)}) and ligand chain IDs ({len(ligand_ids)}).")
    for i, seq in enumerate(ligand_seqs):
        yaml_data['sequences'].append({'ligand': {'id': ligand_ids[i], ligand_type: seq}})
        
    # Process Nucleic Acids
    na_seqs = [s.strip() for s in na_seqs_str.strip().splitlines() if s.strip()]
    na_ids = [i.strip() for i in na_ids_str.strip().split(',') if i.strip()]
    if na_seqs and len(na_seqs) != len(na_ids):
        raise ValueError(f"Mismatch between number of nucleic acids ({len(na_seqs)}) and nucleic acid chain IDs ({len(na_ids)}).")
    for i, seq in enumerate(na_seqs):
        yaml_data['sequences'].append({na_type: {'id': na_ids[i], 'sequence': seq}})

    if not yaml_data['sequences']:
        raise ValueError("No molecules defined. Please provide at least one protein, ligand, or nucleic acid sequence.")

    # Add properties if requested
    if predict_affinity:
        if not binder_chain_id:
            raise ValueError("Binder Chain ID must be specified for affinity prediction.")
        yaml_data['properties'] = [{'affinity': {'binder': binder_chain_id}}]

    # Write the final YAML to a temporary file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix=".yaml", encoding='utf-8') as temp_yaml:
        yaml.dump(yaml_data, temp_yaml, default_flow_style=False, sort_keys=False)
        temp_yaml.flush()
        temp_yaml_path = temp_yaml.name
        temp_yaml_filename = os.path.basename(temp_yaml_path).replace(".yaml", "")

    output_dir = tempfile.mkdtemp()
    
    if force_cpu:
        device = "cpu"
    else:
        device = "gpu" if torch.cuda.is_available() else "cpu"
    
    result = None
    try:
        cmd = [
            "boltz", "predict", temp_yaml_path,
            "--out_dir", output_dir,
            "--accelerator", device,
            "--output_format", "mmcif",
            "--recycling_steps", str(recycling_steps),
            "--sampling_steps", str(sampling_steps)
        ]
        if use_msa_server:
            cmd.append("--use_msa_server")
        if no_kernels:
            cmd.append("--no_kernels")

        print("Running Boltz with command:", " ".join(cmd))
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
             raise subprocess.CalledProcessError(result.returncode, cmd, output=result.stdout, stderr=result.stderr)

    except subprocess.CalledProcessError as e:
        raise Exception(f"Error during Boltz prediction:\nSTDOUT: {e.stdout}\nSTDERR: {e.stderr}")
    except Exception as e:
        raise Exception(f"Error running subprocess: {str(e)}")

    # --- Aggressive file searching logic ---
    structure_files = []
    for root, dirs, files in os.walk(output_dir):
        for file in files:
            if file.endswith(".cif") or file.endswith(".pdb"):
                structure_files.append(os.path.join(root, file))

    if not structure_files:
        # --- Create a directory tree for debugging if no files are found ---
        dir_listing = []
        for root, dirs, files in os.walk(output_dir):
            level = root.replace(output_dir, '').count(os.sep)
            indent = ' ' * 4 * (level)
            dir_listing.append(f'{indent}{os.path.basename(root)}/')
            subindent = ' ' * 4 * (level + 1)
            for f in files:
                dir_listing.append(f'{subindent}{f}')
        
        dir_tree_str = "\n".join(dir_listing) if dir_listing else "Directory is empty."

        raise Exception(
            f"Prediction seems to have completed, but no structure file (.cif or .pdb) was found anywhere inside the output directory.\n\n"
            f"--- Directory Tree of '{output_dir}' ---\n{dir_tree_str}\n\n"
            f"--- Boltz STDOUT ---\n{result.stdout if result else 'N/A'}\n\n"
            f"--- Boltz STDERR ---\n{result.stderr if result else 'N/A'}"
        )
    
    structure_files.sort()
    folded_path = structure_files[0]
    
    prediction_folder = os.path.dirname(folded_path)
    affinity_json_path = os.path.join(prediction_folder, f"affinity_{temp_yaml_filename}.json")
    affinity_results = None
    if os.path.exists(affinity_json_path):
        with open(affinity_json_path, 'r') as f:
            affinity_data = json.load(f)
            affinity_results = (f"Affinity Prediction Results:\n\n"
                                f"Predicted Affinity (log(IC50), µM): {affinity_data.get('affinity_pred_value', 'N/A')}\n"
                                f"Binding Probability: {affinity_data.get('affinity_probability_binary', 'N/A')}")
            
    return folded_path, affinity_results


def fold_protenix(aa_sequence:str, ligand:str=None, ligand_type:str="smiles", use_msa_server:bool=False, seed:int=1337):
    """
    Protein folding using Protenix model
    """
    try:
        import protenix
    except ModuleNotFoundError as e:
        raise Exception(f"Could not import required module: {str(e)}")
    
    json_content = [{"sequences": [{"proteinChain": {"sequence": aa_sequence, "count": 1}}], "name": "pymolfold"}]
    
    if ligand:
        if ligand_type.upper() == "CCD":
            ligand = f"CCD_{ligand[:4]}" if ligand.upper().startswith("CCD_") else f"CCD_{ligand}"
        ligand_dict = {"ligand": {"ligand": ligand, "count": 1}}
        json_content[0]["sequences"].append(ligand_dict)

    with tempfile.NamedTemporaryFile(delete=False, suffix=".json") as temp_json:
        temp_json.write(json.dumps(json_content).encode())
        temp_json_path = temp_json.name

    output_dir = tempfile.mkdtemp()
        
    try:
        cmd = ["protenix", "predict", "--input", temp_json_path, "--out_dir", output_dir, "--seeds", str(seed)]
        if use_msa_server:
            cmd.append("--use_msa_server")
        subprocess.run(cmd, capture_output=True, text=True, check=True)
    except Exception as e:
        raise Exception(f"Error during structure prediction: {str(e)}")

    cif_files = [str(file) for file in Path(os.path.join(output_dir, "pymolfold")).rglob("*") if file.is_file() and str(file).endswith("_sample_0.cif")]
    if not cif_files:
        raise Exception(f"No valid result found in {os.path.join(output_dir, 'pymolfold')}")
    
    return cif_files[0], None

## Database Functions
def get_afdb_structure(database_id):
    afdb_response = requests.get(f"https://alphafold.ebi.ac.uk/api/prediction/{database_id}", headers={"accept": "application/json"})
    if not afdb_response.ok:
        raise Exception("Couldn't get the structure for this AlphaFold database ID.")
    pdb_url = afdb_response.json()[0]['pdbUrl']
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as pdb_file_path:
        with requests.get(pdb_url, stream=True) as r:
            r.raise_for_status()
            for chunk in r.iter_content(chunk_size=8192):
                pdb_file_path.write(chunk)
        return pdb_file_path.name, None

def get_modelarchive_structure(database_id):
    ma_url = f"https://www.modelarchive.org/api/projects/{database_id}?type=basic__model_file_name"
    with tempfile.NamedTemporaryFile(delete=False, suffix=".cif") as cif_file_path:
        with requests.get(ma_url, stream=True) as r:
            r.raise_for_status()
            for chunk in r.iter_content(chunk_size=8192):
                cif_file_path.write(chunk)
        return cif_file_path.name, None

## Coloring Functions
def apply_alphafold_colors(object_name):
    from pymol import cmd
    cmd.set_color("n0", [0.051, 0.341, 0.827]); cmd.set_color("n1", [0.416, 0.796, 0.945]); cmd.set_color("n2", [0.996, 0.851, 0.212]); cmd.set_color("n3", [0.992, 0.490, 0.302])
    cmd.color("n0", f"{object_name} and b < 100"); cmd.color("n1", f"{object_name} and b < 90"); cmd.color("n2", f"{object_name} and b < 70"); cmd.color("n3", f"{object_name} and b < 50")

def apply_bfactor_colors(object_name):
    from pymol import cmd
    cmd.spectrum("b", palette="rainbow", selection=object_name)

## Main Dialog
def make_dialog():
    from pymol import cmd
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi

    dialog = QtWidgets.QDialog()
    uifile = os.path.join(os.path.dirname(__file__), 'widget.ui')
    form = loadUi(uifile, dialog)

    def get_uniprot_sequence():
        uniprot_id = form.input_uniprot_id.text()
        if not uniprot_id:
            QtWidgets.QMessageBox.warning(form, "Error", "Please enter a valid UniProt ID.")
            return
        try:
            uniprot_response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}", headers={"accept": "application/json"}, params={"fields": ["sequence"]})
            uniprot_response.raise_for_status()
            form.input_aa_seq.setPlainText(uniprot_response.json()['sequence']['value'])
        except Exception as e:
            QtWidgets.QMessageBox.warning(form, "Error", f"Couldn't get sequence for UniProt ID {uniprot_id}: {e}")

    def run():
        try:
            model_name = form.input_list_models.currentText()
            database_id = form.input_database_id.text()
            seed = int(form.input_seed.text())
            af_coloring = form.input_af_coloring.isChecked()
            bfactor_coloring = form.input_bfactor_coloring.isChecked()

            folded_path, additional_results = None, None

            if model_name == "Boltz":
                protein_seqs = form.input_boltz_protein_seqs.toPlainText()
                protein_ids = form.input_boltz_protein_ids.text()
                ligand_seqs = form.input_boltz_ligand_seqs.toPlainText()
                ligand_ids = form.input_boltz_ligand_ids.text()
                ligand_type = form.input_boltz_ligand_type.currentText()
                na_seqs = form.input_boltz_na_seqs.toPlainText()
                na_ids = form.input_boltz_na_ids.text()
                na_type = form.input_boltz_na_type.currentText()

                recycling_steps = int(form.input_boltz_recycling_steps.text())
                sampling_steps = int(form.input_boltz_sampling_steps.text())
                use_msa_server = form.input_boltz_use_msa_server.isChecked()
                is_cyclic = form.input_boltz_is_cyclic.isChecked()
                predict_affinity = form.input_boltz_predict_affinity.isChecked()
                binder_chain_id = form.input_boltz_binder_chain_id.text().strip()
                force_cpu = form.input_boltz_force_cpu.isChecked()
                no_kernels = form.input_boltz_no_kernels.isChecked()

                folded_path, additional_results = fold_boltz(
                    protein_seqs, protein_ids, ligand_seqs, ligand_ids, ligand_type,
                    na_seqs, na_ids, na_type, use_msa_server, recycling_steps,
                    sampling_steps, is_cyclic, predict_affinity, binder_chain_id, force_cpu, no_kernels
                )
            else:
                # Handle other models which use the original input fields
                system_definition = form.input_aa_seq.toPlainText()
                ligand_sequence = form.input_ligand.toPlainText().strip()
                ligand_type = form.input_ligand_type.currentText() if ligand_sequence else None

                if model_name not in ["AlphaFoldDB", "ModelArchive"] and not system_definition:
                    QtWidgets.QMessageBox.warning(form, "Error", "Please enter a sequence.")
                    return

                if model_name.startswith("esm3"):
                    esm_token = form.input_esm_token.text(); esm_temp = float(form.input_esm_temp.text()); esm_nsteps = int(form.input_esm_nsteps.text())
                    folded_path, additional_results = fold_esm(model_name, system_definition, temperature=esm_temp, num_steps=esm_nsteps, token=esm_token)
                elif model_name == "chai-1":
                    chai_recycling_steps = int(form.input_chai_recycling_steps.text()); chai_diffusion_steps = int(form.input_chai_diffusion_steps.text())
                    folded_path, additional_results = fold_chai(system_definition, ligand=ligand_sequence, ligand_type=ligand_type, num_trunk_recycles=chai_recycling_steps, num_diffn_timesteps=chai_diffusion_steps, seed=seed)
                elif model_name == "protenix":
                    protenix_use_msa = form.input_protenix_use_msa_server.isChecked()
                    folded_path, additional_results = fold_protenix(system_definition, ligand=ligand_sequence, ligand_type=ligand_type, use_msa_server=protenix_use_msa, seed=seed)
                elif model_name == "AlphaFoldDB":
                    folded_path, additional_results = get_afdb_structure(database_id)
                elif model_name == "ModelArchive":
                    folded_path, additional_results = get_modelarchive_structure(database_id)
                else:
                    QtWidgets.QMessageBox.critical(form, "Error", f"Not a supported model name: {str(model_name)}"); return

            if not folded_path:
                QtWidgets.QMessageBox.critical(form, "Error", "No folded structure was returned."); return
            
            object_name = f"folded_structure_{''.join(random.choices(string.ascii_lowercase + string.digits, k=3))}"
            cmd.load(folded_path, object_name)
            
            if af_coloring: apply_alphafold_colors(object_name)
            if bfactor_coloring: apply_bfactor_colors(object_name)
            
            success_msg = f"Structure loaded into PyMOL from {model_name}!"
            if additional_results:
                QtWidgets.QMessageBox.information(form, "Success & Affinity Results", f"{success_msg}\n\n{additional_results}")
            else:
                QtWidgets.QMessageBox.information(form, "Success", success_msg)
        
        except Exception as e:
            QtWidgets.QMessageBox.critical(form, "Error", f"An error occurred: {str(e)}")

    def update_ui():
        model_name = form.input_list_models.currentText()
        is_database_model = model_name in ["AlphaFoldDB", "ModelArchive"]
        is_boltz = model_name == "Boltz"
        
        # Switch between the new Boltz form and the standard input form
        form.input_stacked_widget.setCurrentIndex(1 if is_boltz else 0)

        # Hide/show UniProt fetcher (not used with new Boltz form)
        form.label_uniprot_id.setVisible(not is_database_model and not is_boltz)
        form.input_uniprot_id.setVisible(not is_database_model and not is_boltz)
        form.button_uniprot_id.setVisible(not is_database_model and not is_boltz)
        
        # Hide/show database ID input
        form.input_database_id.setVisible(is_database_model)
        form.label_database_id.setVisible(is_database_model)
        
        # Hide dedicated legacy ligand input for Boltz
        is_ligand_supported_legacy = model_name in ["chai-1", "protenix"]
        form.input_ligand.setVisible(is_ligand_supported_legacy)
        form.label_ligand.setVisible(is_ligand_supported_legacy)
        form.input_ligand_type.setVisible(is_ligand_supported_legacy)
        form.label_ligand_type.setVisible(is_ligand_supported_legacy)

        # Show/hide settings groups
        form.group_esm_settings.setVisible(model_name.startswith("esm3"))
        form.group_chai_settings.setVisible(model_name=="chai-1")
        form.group_boltz_settings.setVisible(is_boltz)
        form.group_protenix_settings.setVisible(model_name=="protenix")
        form.label_settings.setVisible(not is_database_model)

        # Enable the Boltz settings group so its children can be clicked
        form.group_boltz_settings.setEnabled(is_boltz)
        
        # Enable binder ID input only when affinity prediction is checked
        if is_boltz:
            form.input_boltz_binder_chain_id.setEnabled(form.input_boltz_predict_affinity.isChecked())

        # Configure main "Fold" button
        form.button_fold.setEnabled(not model_name.startswith("--"))
        form.button_fold.setText("Download" if is_database_model else "Fold")
        dialog.adjustSize()

    form.button_uniprot_id.clicked.connect(get_uniprot_sequence)
    form.button_fold.clicked.connect(run)
    form.button_close.clicked.connect(dialog.close)
    form.input_list_models.currentIndexChanged.connect(update_ui)
    form.input_boltz_predict_affinity.stateChanged.connect(update_ui)
    
    # Trigger initial UI update to set the correct state
    update_ui()
    return dialog
