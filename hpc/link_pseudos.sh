#!/bin/bash
set -e
PSEUDO_SOURCE='/mgpfs/home/lala002/pseudo'
mkdir -p pseudo
ln -sfn "$PSEUDO_SOURCE/Ti.pbe-spn-kjpaw_psl.1.0.0.UPF" pseudo/Ti.UPF
ln -sfn "$PSEUDO_SOURCE/O.pbe-n-kjpaw_psl.1.0.0.UPF" pseudo/O.UPF
ln -sfn "$PSEUDO_SOURCE/C.pbe-n-kjpaw_psl.1.0.0.UPF" pseudo/C.UPF
ln -sfn "$PSEUDO_SOURCE/H.pbe-kjpaw_psl.1.0.0.UPF" pseudo/H.UPF
ln -sfn "$PSEUDO_SOURCE/N.pbe-n-kjpaw_psl.1.0.0.UPF" pseudo/N.UPF
ls -lh pseudo/
