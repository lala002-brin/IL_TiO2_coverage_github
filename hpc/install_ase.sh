#!/bin/bash
set -e
python3 -m pip install ase
python3 -c "import ase; print('ASE OK:',ase.__version__)"
