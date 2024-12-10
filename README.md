# PyMOLfold
Plugin for folding sequences directly in PyMOL using EvolutionaryScale's ESM3 model.

<h3 align="right">Colby T. Ford, Ph.D. | Tuple, LLC</h3>

![Plugin Screenshot](screenshot.png)

## Plugin Installation

1. Download the [pymolfold.zip](pymolfold.zip) file from this repository.
2. Open PyMOL and in the menu bar, go to `Plugin` > `Plugin Manager`.
3. On the `Install New Plugin` tab, click the `Choose file...` button under "Install from local file".
4. Locate the .zip file and open it.

Once installed, navigate to `Plugin` > `PyMOLfold` in the menu bar.
Then, in the dialog box, simply paste in your API key and amino acid sequence you wish to fold.

> [!NOTE]
> Note: This plugin does not include any of the model package(s). You'll need to install them into the PyMOL conda environment.


## Package Installation

You can install PyPI libraries from the PyMOL command line. However, you can run the pip install commands from your normal terminal, if you know the path to PyMOL's Python executable. To find it, use the `sys` library as follows:

```python
PyMOL>import sys
PyMOL>print(sys.executable)
/Applications/PyMOL.app/Contents/bin/python
```
Then, to run a command in your normal terminal, simply point to the PyMOL Python executable. For example:

```bash
/Applications/PyMOL.app/Contents/bin/python -m pip install <package_name>
```

Here are the example install commands for the supported models:

- ESM3: `pip install esm`

> [!NOTE]
> Note: To use an ESM3 model, you'll need an API key from [Evolutionary Scale Forge](https://forge.evolutionaryscale.ai/). Otherwise, the plugin can't download the folding model.


## Feature Roadmap

- [X] Basic folding with `esm3-small-2024-08`
- [ ] UI controls for temperature and steps.
- [ ] Support for multiple chains (or FASTA input format).
- [ ] UI Dropdown to select different models.
