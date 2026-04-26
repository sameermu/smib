#!/usr/bin/env bash
# One-shot setup for the smib project on macOS / Linux.
#
# Run from the smib repo root:
#     bash setup.sh
#
# What it does:
#   1. Creates a venv at  ./.venv   (using whichever python3 is on PATH)
#   2. Activates the venv and installs requirements.txt
#   3. Runs the test suite to confirm everything works
#   4. Prints what to do next in PyCharm

set -e

cd "$(dirname "$0")"

if [ ! -d ".venv" ]; then
    echo ">> Creating venv at ./.venv"
    python3 -m venv .venv
else
    echo ">> Reusing existing ./.venv"
fi

# shellcheck disable=SC1091
source .venv/bin/activate

echo ">> Upgrading pip"
python -m pip install --upgrade pip --quiet

echo ">> Installing requirements"
pip install -r requirements.txt --quiet

echo ">> Installing smib in editable mode (so notebooks can import it)"
pip install -e . --quiet

echo ">> Running test battery"
pytest tests/ -q || {
    echo "!! Tests failed. Stopping."
    exit 1
}

echo ""
echo "================================================================"
echo "Setup complete."
echo ""
echo "Next steps in PyCharm:"
echo "  1. File -> Open    ->  $(pwd)"
echo "       (open the SMIB ROOT, not the notebooks subfolder)"
echo "  2. PyCharm will auto-detect .venv and prompt you to use it."
echo "       Click 'Use this interpreter' / 'OK'."
echo "  3. Open notebooks/phase1_gencls.ipynb and run all cells."
echo "================================================================"
