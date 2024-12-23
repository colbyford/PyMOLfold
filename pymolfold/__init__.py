'''
PyMOL Protein Folding Plugin

By Colby T. Ford, Ph.D.
License: GPLv3
'''

from __future__ import absolute_import
from __future__ import print_function
import os, tempfile, random, string, sys, subprocess, json
from pathlib import Path

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


## Folding Functions

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

    ## Generate the protein structure
    structure_prediction_config = GenerationConfig(
        track="structure",
        num_steps=num_steps,
        temperature=temperature,
    )

    structure_prediction_prompt = ESMProtein(sequence=aa_sequence)

    structure_prediction = model.generate(
        structure_prediction_prompt,
        structure_prediction_config
    )

    structure_prediction_chain = structure_prediction.to_protein_chain()

    pdb_string = structure_prediction_chain.to_pdb_string()

    ## Save the output PDB file temporarily
    with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb:
        temp_pdb.write(pdb_string.encode())
        temp_pdb_path = temp_pdb.name

    return temp_pdb_path

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
    
    ## Start building FASTA content
    fasta_content = f">protein|name=chain_A\n{aa_sequence}\n"
    
    ## Add ligand if provided
    if ligand and ligand_type:
        fasta_content += f">ligand|name=chain_B\n{ligand}\n"

    ## Create temp fasta file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
        temp_fasta.write(fasta_content.encode())
        temp_fasta_path = temp_fasta.name

    ## Create temp output directory
    output_dir = tempfile.mkdtemp()

    ## Detect devices
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    ## Run inferencing
    candidates = run_inference(
        fasta_file=Path(temp_fasta_path),
        output_dir=Path(output_dir),
        num_trunk_recycles=num_trunk_recycles,
        num_diffn_timesteps=num_diffn_timesteps,
        seed=seed,
        device=device,
        use_esm_embeddings=True
    )

    if not candidates.cif_paths:
        raise ValueError("No structure files were generated")

    return candidates.cif_paths[0]

## Boltz Folding
def fold_boltz(aa_sequence:str, ligand:str=None, ligand_type:str=None, use_msa_server:bool=False, recycling_steps:int=3, sampling_steps:int=200):
    """
    Protein folding using Boltz-1 model
    """
    try:
        import boltz
        import torch
    except ModuleNotFoundError as e:
        raise Exception(f"Could not import required module: {str(e)}")

    ## Start building FASTA content
    fasta_content = f">A|protein|empty\n{aa_sequence}\n"
    
    ## Add ligand if provided
    if ligand and ligand_type:
        fasta_content += f">B|{ligand_type}|\n{ligand}\n"

    ## Create temp fasta file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
        temp_fasta.write(fasta_content.encode())
        temp_fasta_path = temp_fasta.name
        temp_fasta_filename = os.path.basename(temp_fasta_path).replace(".fasta", "")

    ## Create temp output directory
    output_dir = tempfile.mkdtemp()
    
    ## Set device
    device = "gpu" if torch.cuda.is_available() else "cpu"
    
    try:
        ## Run boltz command
        cmd = [
            "boltz",
            "predict",
            temp_fasta_path,
            "--out_dir", output_dir,
            "--accelerator", device,
            "--output_format", "pdb",
            "--recycling_steps", str(recycling_steps),
            "--sampling_steps", str(sampling_steps)
        ]

        print("Running Boltz with command:", " ".join(cmd))

        if use_msa_server:
            cmd.append("--use_msa_server")
        
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    except Exception as e:
        raise Exception(f"Error during structure prediction: {str(e)}")
    
    ## Get the path to the folded PDB file
    folded_pdb_path = os.path.join(output_dir,
                        f"boltz_results_{temp_fasta_filename}",
                        "predictions",
                        temp_fasta_filename,
                        f"{temp_fasta_filename}_model_0.pdb")
    
    if not os.path.exists(folded_pdb_path):
        raise Exception(f"Expected output file not found: {folded_pdb_path}")
    
    return folded_pdb_path

def fold_protenix(aa_sequence:str, ligand:str=None, ligand_type:str="smiles", use_msa_server:bool=False, seed:int=1337):
    """
    Protein folding using Protenix model
    """
    try:
        import protenix
        import torch
    except ModuleNotFoundError as e:
        raise Exception(f"Could not import required module: {str(e)}")
    
    ## Build the JSON body
    json_content = [{"sequences": [
            {"proteinChain": {
                "sequence": aa_sequence,
                "count": 1
                }
            }
        ],
        "name": "pymolfold"
        }]
    
    ## If ligand provided, add it to the JSON
    if ligand:
        if ligand_type.upper() == "CCD":
            ligand = f"CCD_{ligand[:4]}" if ligand.upper().startswith("CCD_") else f"CCD_{ligand}"
        ligand_dict = {
            "ligand": {
                "ligand": ligand,
                "count": 1
            }
        }
        json_content[0]["sequences"].append(ligand_dict)

    ## Create temp json file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".json") as temp_json:
        temp_json.write(json.dumps(json_content).encode())
        temp_json_path = temp_json.name

    ## Create temp output directory
    output_dir = tempfile.mkdtemp()
        
    try:
        ## Run protenix command
        cmd = [
            "protenix",
            "predict",
            "--input", temp_json_path,
            "--out_dir", output_dir,
            "--seeds", str(seed)
        ]
        if use_msa_server:
            cmd.append("--use_msa_server")

        result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    except Exception as e:
        raise Exception(f"Error during structure prediction: {str(e)}")

    cif_files = [str(file) for file in Path(os.path.join(output_dir, "pymolfold")).rglob("*")
                 if file.is_file() and str(file).endswith("_sample_0.cif")]
    if len(cif_files) > 0:
        folded_cif_path = cif_files[0]
    else:
        raise Exception(f"No valid result found in {os.path.join(output_dir, 'pymolfold')}")
    if not os.path.exists(folded_cif_path):
        raise Exception(f"Expected output file not found: {folded_cif_path}")
    
    return folded_cif_path

def apply_alphafold_colors(object_name):
    """
    Apply AlphaFold-style coloring to the structure
    Credit: Konstantin Korotkov
    """
    from pymol import cmd
    
    ## Define AlphaFold colors
    cmd.set_color("n0", [0.051, 0.341, 0.827])  # High confidence (dark blue)
    cmd.set_color("n1", [0.416, 0.796, 0.945])  # Good confidence (light blue)
    cmd.set_color("n2", [0.996, 0.851, 0.212])  # Medium confidence (yellow)
    cmd.set_color("n3", [0.992, 0.490, 0.302])  # Low confidence (orange)
    
    ## Apply coloring based on B-factor ranges
    cmd.color("n0", f"{object_name} and b < 100")
    cmd.color("n1", f"{object_name} and b < 90")
    cmd.color("n2", f"{object_name} and b < 70")
    cmd.color("n3", f"{object_name} and b < 50")


## Main Dialog
def make_dialog():
    ## Entrypoint to the PyMOL API
    from pymol import cmd

    ## Pymol.Qt provides the PyQt5 interface, but may support PyQt4
    ## and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi

    ## Create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'widget.ui')
    form = loadUi(uifile, dialog)

    ## Hide UI controls
    form.input_ligand.setVisible(False)
    form.label_ligand.setVisible(False)
    form.input_ligand_type.setVisible(False)
    form.label_ligand_type.setVisible(False)
    form.group_esm_settings.setVisible(False)
    form.group_chai_settings.setVisible(False)
    form.group_boltz_settings.setVisible(False)
    form.group_protenix_settings.setVisible(False)

    ## Callback for the "Fold" button
    def run():

        ## Get basic form data
        model_name = form.input_list_models.currentText()
        aa_sequence = form.input_aa_seq.toPlainText()
        ligand_sequence = form.input_ligand.toPlainText().strip()
        ligand_type = form.input_ligand_type.currentText() if ligand_sequence else None
        seed = int(form.input_seed.text())
        af_coloring = form.input_af_coloring.isChecked()

        if not aa_sequence:
            QtWidgets.QMessageBox.warning(form, "Error", "Please enter a valid amino acid sequence.")
            return

        try:
            if model_name.startswith("esm3"):
                ## ESM Parameters
                esm_token = form.input_esm_token.text()
                esm_temp = float(form.input_esm_temp.text())
                esm_nsteps = int(form.input_esm_nsteps.text())

                folded_pdb_path = fold_esm(model_name,
                                           aa_sequence,
                                           temperature=esm_temp,
                                           num_steps=esm_nsteps,
                                           token=esm_token)
                
            elif model_name == "boltz-1":
                ## Boltz Parameters
                boltz_recycling_steps = int(form.input_boltz_recycling_steps.text())
                boltz_sampling_steps = int(form.input_boltz_sampling_steps.text())
                boltz_use_msa_server = form.input_boltz_use_msa_server.isChecked()

                folded_pdb_path = fold_boltz(aa_sequence,
                                             ligand=ligand_sequence,
                                             ligand_type=ligand_type,
                                             use_msa_server=boltz_use_msa_server,
                                             recycling_steps=boltz_recycling_steps,
                                             sampling_steps=boltz_sampling_steps)
                
            elif model_name == "chai-1":
                ## Chai Parameters
                chai_recycling_steps = int(form.input_chai_recycling_steps.text())
                chai_diffusion_steps = int(form.input_chai_diffusion_steps.text())

                folded_pdb_path = fold_chai(aa_sequence,
                                            ligand=ligand_sequence,
                                            ligand_type=ligand_type,
                                            num_trunk_recycles=chai_recycling_steps,
                                            num_diffn_timesteps=chai_diffusion_steps,
                                            seed=seed)
            
            elif model_name == "protenix":
                ## Protenix Parameters
                protenix_use_msa = form.input_protenix_use_msa_server.isChecked()

                folded_pdb_path = fold_protenix(aa_sequence,
                                                ligand=ligand_sequence,
                                                ligand_type=ligand_type,
                                                use_msa_server=protenix_use_msa,
                                                seed=seed)

            else:
                QtWidgets.QMessageBox.critical(form, "Error", f"Not a supported model name: {str(model_name)}")
                return

            ## Load the folded structure into PyMOL
            if not folded_pdb_path:
                QtWidgets.QMessageBox.critical(form, "Error", "No folded structure was returned.")
                return
            
            ## Generate a unique object name
            object_name = f"folded_structure_{''.join(random.choices(string.ascii_lowercase + string.digits, k=3))}"
            cmd.load(folded_pdb_path, object_name)
            
            ## Apply AlphaFold-style coloring
            if af_coloring:
                apply_alphafold_colors(object_name)
            
            QtWidgets.QMessageBox.information(form, "Success", f"Sequence folded with {model_name} and structure loaded into PyMOL!")
        
        except Exception as e:
            QtWidgets.QMessageBox.critical(form, "Error", f"An error occurred: {str(e)}")

    def update_ui():
        ## Update the UI based on the selected model
        model_name = form.input_list_models.currentText()

        ## Ligand supported models
        ligand_supported_models = ["chai-1", "boltz-1", "protenix"]
        form.input_ligand.setVisible(model_name in ligand_supported_models)
        form.label_ligand.setVisible(model_name in ligand_supported_models)
        form.input_ligand_type.setVisible(model_name in ligand_supported_models)
        form.label_ligand_type.setVisible(model_name in ligand_supported_models)

        ## Group boxes for settings
        form.group_esm_settings.setVisible(model_name.startswith("esm3"))
        form.group_chai_settings.setVisible(model_name=="chai-1")
        form.group_boltz_settings.setVisible(model_name=="boltz-1")
        form.group_protenix_settings.setVisible(model_name=="protenix")

    ## Button callbacks
    form.button_fold.clicked.connect(run)
    form.button_close.clicked.connect(dialog.close)

    ## Model selection callback
    form.input_list_models.currentIndexChanged.connect(update_ui)

    return dialog