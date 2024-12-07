'''
PyMOL Protein Folding Plugin

(c) Tuple, LLC by Colby T. Ford, Ph.D.
License: GPLv3
'''

from __future__ import absolute_import
from __future__ import print_function
import os


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
        
        ## Get form data
        aa_sequence = form.input_aa_seq.text()
        token = form.input_token.text()

        if not aa_sequence:
            QtWidgets.QMessageBox.warning(form, "Error", "Please enter a valid amino acid sequence.")
            return

        try:
            ## Import necessary modules
            import tempfile
            from esm.sdk import client
            from esm.sdk.api import ESMProtein, GenerationConfig

            model = client(model="esm3-small-2024-08", url="https://forge.evolutionaryscale.ai", token=token)

            structure_prediction_config = GenerationConfig(
                track="structure",
                num_steps=8,
                temperature=0.7,
            )

            structure_prediction_prompt = ESMProtein(sequence=aa_sequence)

            structure_prediction = model.generate(
                structure_prediction_prompt, structure_prediction_config
            )

            structure_prediction_chain = structure_prediction.to_protein_chain()

            pdb_string = structure_prediction_chain.to_pdb_string()

            ## Save the output PDB file temporarily
            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdb") as temp_pdb:
                temp_pdb.write(pdb_string.encode())
                temp_pdb_path = temp_pdb.name

            # Load the folded structure into PyMOL
            cmd.load(temp_pdb_path, "folded_structure")
            QtWidgets.QMessageBox.information(form, "Success", "Structure folded and loaded into PyMOL!")
        
        except Exception as e:
            QtWidgets.QMessageBox.critical(form, "Error", f"An error occurred: {str(e)}")

    ## Button callbacks
    form.button_fold.clicked.connect(run)
    form.button_close.clicked.connect(dialog.close)

    return dialog