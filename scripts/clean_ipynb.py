# This script seves to clean the output of all .ipynb
# as they might contain sensitive information.
# the .gitattributes file is configured to run this script
# using: *.ipynb filter=clean_ipynb
# where clean_ipynb is defined in .git/config as:
# git config filter.clean_ipynb.clean 'python clean_ipynb.py %f'

import nbformat


def clean_notebook(notebook_path):
    with open(notebook_path, "r") as f:
        nb = nbformat.read(f, as_version=4)

    # Remove output cells
    for cell in nb.cells:
        if "outputs" in cell:
            cell["outputs"] = []
        if "execution_count" in cell:
            cell["execution_count"] = None

    # Write the modified notebook back to file
    with open(notebook_path, "w") as f:
        nbformat.write(nb, f)


if __name__ == "__main__":
    import sys

    notebook_path = sys.argv[1]
    clean_notebook(notebook_path)
