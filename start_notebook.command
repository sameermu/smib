#!/usr/bin/env bash
# Double-click this file in Finder to launch the smib Jupyter notebook.
#
# What it does:
#   1. cd into the smib repo (wherever this file lives)
#   2. Activate the .venv created by setup.sh
#   3. Launch Jupyter, landing on the notebooks/ file browser so you
#      can pick any phase notebook.
#
# Press Ctrl+C twice in this Terminal window to stop Jupyter when done.

set -e

cd "$(dirname "$0")"

if [ ! -d ".venv" ]; then
    echo ">> No .venv found. Run 'bash setup.sh' first."
    read -n 1 -s -r -p "Press any key to close..."
    exit 1
fi

# shellcheck disable=SC1091
source .venv/bin/activate

# Land on the notebooks/ file browser so the user can choose between
# phase1_gencls.ipynb, phase2_0_genrou.ipynb, etc.
jupyter notebook --notebook-dir=notebooks
