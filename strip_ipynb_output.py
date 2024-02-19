# This file deletes all the output cells from a Jupyter notebook.
# whenever a commit is made, this script is called to clean the notebook
# and remove all the output cells. This is done to protect sensitive data.


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
