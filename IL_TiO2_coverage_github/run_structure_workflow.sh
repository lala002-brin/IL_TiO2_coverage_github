#!/usr/bin/env bash
set -e

python3 scripts/01_check_IL_TiO2.py
python3 scripts/03_test_distance.py
python3 scripts/04_make_1IL_TiO2.py
python3 scripts/05_test_2IL_TiO2.py
python3 scripts/06_check_periodic_IL.py
python3 scripts/07_make_2IL_TiO2.py
python3 scripts/08_test_3IL_monolayer.py
python3 scripts/09_search_3IL_staggered.py
python3 scripts/10_test_3IL_multilayer.py
python3 scripts/11_make_3IL_multilayer.py
