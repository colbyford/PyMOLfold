'''
PyMOL Protein Folding Plugin

(c) Tuple, LLC by Colby T. Ford, Ph.D.
License: GPLv3
'''

from __future__ import absolute_import
from __future__ import print_function
import os, tempfile, random, string

def __init_plugin__(app=None):
    '''
    Add an entry to the PyMOL "Plugin" menu
    '''
    from pymol.plugins import addmenuitemqt
    addmenuitemqt('PyMOLfold', run_plugin_gui)


## global reference to avoid garbage collection of our dialog
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
def fold_esm(model_name, aa_sequence, token):
    """
    Protein folding using ESM models
    """
    from esm.sdk import client
    from esm.sdk.api import ESMProtein, GenerationConfig

    try:
        model = client(model=model_name, url="https://forge.evolutionaryscale.ai", token=token)
    except Exception as e:
        raise Exception(f"Error getting ESM model with token: {str(e)}")

    ## Generate the protein structure
    structure_prediction_config = GenerationConfig(
        track="structure",
        num_steps=8,
        temperature=0.7,
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
def fold_chai(aa_sequence):
    """
    Protein folding using Chai models
    """
    from chai_lab.chai1 import run_inference
    import torch

    fasta_line = f">aa_sequence\n{aa_sequence}"

    ## Create temp fasta file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
        temp_fasta.write(fasta_line.encode())
        temp_fasta_path = temp_fasta.name

    ## Create temp output directory
    output_dir = tempfile.mkdtemp()

    ## Detect devices
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    ## Run inferencing
    candidates = run_inference(
        fasta_file=temp_fasta_path,
        output_dir=output_dir,
        num_trunk_recycles=3,
        num_diffn_timesteps=200,
        seed=1337,
        device=device,
        use_esm_embeddings=True
    )

    cif_paths = candidates.cif_paths

    return cif_paths[0]


## Boltz Folding
def fold_boltz(aa_sequence):
    """
    Protein folding using Boltz-1 model
    """
    fasta_line = f">A|protein|empty\n{aa_sequence}"

    ## Create temp fasta file
    with tempfile.NamedTemporaryFile(delete=False, suffix=".fasta") as temp_fasta:
        temp_fasta.write(fasta_line.encode())
        temp_fasta_path = temp_fasta.name
        temp_fasta_filename = os.path.basename(temp_fasta_path).replace(".fasta", "")

    ## Create temp output directory
    output_dir = tempfile.mkdtemp()

    ## Run Inferencing
    import torch
    accelerator = torch.device("gpu" if torch.cuda.is_available() else "cpu")
    boltz_run = os.system(f"boltz predict {temp_fasta_path} --out_dir {output_dir} --output_format pdb --use_msa_server --accelerator {accelerator}")

    if boltz_run != 0:
        raise Exception("Error running Boltz Model.")
    
    ## Get the path to the folded PDB file
    folded_pdb_path = os.path.join(output_dir, f"boltz_results_{temp_fasta_filename}", "predictions", temp_fasta_filename, f"{temp_fasta_filename}_model_0.pdb")

    return folded_pdb_path


## Main Dialog
def make_dialog():
    ## Entrypoint to the PyMOL API
    from pymol import cmd

    # pymol.Qt provides the PyQt5 interface, but may support PyQt4
    # and/or PySide as well
    from pymol.Qt import QtWidgets
    from pymol.Qt.utils import loadUi

    ## create a new Window
    dialog = QtWidgets.QDialog()

    # populate the Window from our *.ui file which was created with the Qt Designer
    uifile = os.path.join(os.path.dirname(__file__), 'widget.ui')
    form = loadUi(uifile, dialog)

    ## callback for the "Fold" button
    def run():

        ## Supported Models
        models = [
            "esm3-small-2024-08",
            "esm3-open-2024-03",
            "esm3-medium-2024-08",
            "esm3-large-2024-03",
            "esm3-medium-multimer-2024-09",
            "boltz-1",
            "chai-1"
        ]
        
        ## Get form data
        # model_name = form.input_model_name.text()
        model_name = models[form.input_list_models.currentRow()]
        aa_sequence = form.input_aa_seq.toPlainText()
        token = form.input_token.text()

        if not aa_sequence:
            QtWidgets.QMessageBox.warning(form, "Error", "Please enter a valid amino acid sequence.")
            return

        try:
            if model_name.startswith("esm3"):
                folded_pdb_path = fold_esm(model_name, aa_sequence, token)
            elif model_name == "chai-1":
                folded_pdb_path = fold_chai(aa_sequence)
            elif model_name == "boltz-1":
                folded_pdb_path = fold_boltz(aa_sequence)
            else:
                QtWidgets.QMessageBox.critical(form, "Error", f"Not a supported model name: {str(model_name)}")

            ## Load the folded structure into PyMOL
            if not folded_pdb_path:
                QtWidgets.QMessageBox.critical(form, "Error", "No folded structure was returned.")
                return
            cmd.load(folded_pdb_path, f"folded_structure_{''.join(random.choices(string.ascii_lowercase + string.digits, k=3))}")
            QtWidgets.QMessageBox.information(form, "Success", "Structure folded and loaded into PyMOL!")
        
        except Exception as e:
            QtWidgets.QMessageBox.critical(form, "Error", f"An error occurred: {str(e)}")

    ## Button callbacks
    form.button_fold.clicked.connect(run)
    form.button_close.clicked.connect(dialog.close)

    return dialog