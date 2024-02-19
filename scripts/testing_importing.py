import os

print("path:", os.getcwd())
# from src.initialisation import load_surfplan

import sys
import os

# Add the parent directory of `src` to sys.path
script_dir = os.path.dirname(__file__)  # Directory of the script
parent_dir = os.path.dirname(script_dir)  # Parent directory (project_root)
sys.path.append(parent_dir)

from src.initialisation import load_surfplan

print("path:", os.getcwd())
