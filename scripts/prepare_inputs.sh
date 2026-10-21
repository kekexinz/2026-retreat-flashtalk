#!/usr/bin/env bash
# Inputs for the 2DTM movie, built from what is in data/lysate_12781 into work/:
#   1. unpack the micrograph (.mrc.xz) and check its SHA256 against manifest.json
#   2. simulate the 60S template from the PDB with cisTEM simulate (448^3, 1.66 A,
#      per-atom PDB B-factors, no added B), then crop it to the 320^3 Accuracy box
#   3. prep.py: micrograph PNG, template mesh, detections -> work/lysate_12781
# Steps 1 and 2 are skipped when their outputs already exist.
set -euo pipefail
cd "$(dirname "$0")/.."

PY=${PY:-python}
D=data/lysate_12781
I=work/inputs
W=work/lysate_12781
MIC=$I/micrograph_12781_native_0p83A.mrc
TMPL=$I/template_5lzv_60S_accuracy_1p66A_320.mrc
mkdir -p $I

if [ ! -f $MIC ]; then
  xz -dkc $D/micrograph_12781_native_0p83A.mrc.xz > $MIC.part && mv $MIC.part $MIC
fi
echo "728b51fc837bfe47e8eec2cf77dc5c600abdaf886d511be26eb950408eb3fef5  $MIC" | sha256sum -c -

if [ ! -f $TMPL ]; then
  # template_simulator.py writes next to the PDB, so give it a copy in work/
  cp $D/5lzv_60S_minCtermL24.pdb $I/
  "$PY" scripts/template_simulator.py --pdb $I/5lzv_60S_minCtermL24.pdb \
    --output_size 448 --pixel_size 1.66 --bfactor 0 --pdb_bfactor_scaling 1
  "$PY" scripts/crop_template.py $I/5lzv_60S_minCtermL24_pdbB1_bf0_448_1.66Apix.mrc $TMPL
fi

"$PY" scripts/prep.py \
  --micrograph $MIC --template $TMPL \
  --star $D/accuracy_peaks.star --match "" \
  --defocus 15589 --mass-kda 2200 --mesh-lowpass 0 \
  --out $W
