#!/usr/bin/env bash
# Double-click this file in Finder to launch the smib Jupyter notebook.
#
# What it does:
#   1. cd into the smib repo (wherever this file lives)
#   2. Activate the .venv created by setup.sh
#   3. Launch Jupyter on the Phase 1 notebook (will open your browser)
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

# Launch on the Phase 1 notebook by default. Edit this line if you
# want to land on a different notebook in the future, or remove the
# trailing path to land on the file browser.
jupyter notebook notebooks/phase1_gencls.ipynb
