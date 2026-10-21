import argparse
import os
from subprocess import Popen, PIPE
from collections import namedtuple

SimulateParameters2 = namedtuple("SimulateParameters2", [
    "output_filename",
    "make_3d",
    "output_size",
    "n_threads",
    "input_pdb_file",
    "para1",
    "pixel_size",
    "linear_scaling_of_PDB_bfactors",
    "per_atom_bfactor",
    "exposure_per_frame",
    "n_frames",
    "expert_options",
    "para4",
    "para5",
    "KV",
    "CS",
    "OA",
    "para6",
    "pre_exposure",
    "para7",
    "para8",
    "para9",
    "para10",
    "para11",
    "para12",
    "para13",
    "para14",
    "para15",
    "para16",
    "para17",
    "para18"
], defaults=[
    "out.mrc",
    "yes",
    512,
    1,
    "in.pdb",
    "no",
    0.8,
    0.0,
    55.0,
    1.0,
    30.0,
    "yes",
    "1",
    0.0,
    300.0,
    2.7,
    100.0,
    0.0,
    0.0,
    32,
    10,
    2048,
    0.0,
    2.0,
    0.1,
    0.0001,
    0.0,
    0.0,
    0.0,
    0.0,
    0.0
])

parser = argparse.ArgumentParser(description="Run cisTEM simulator on a PDB file.")
parser.add_argument("--pdb", required=True, help="Path to input PDB file")
parser.add_argument("--output_size", type=int, required=True, help="Output volume size in pixels")
parser.add_argument("--pixel_size", type=float, required=True, help="Pixel size in Angstroms")
parser.add_argument("--bfactor", type=float, required=True, help="B-factor added to every atom")
parser.add_argument("--pdb_bfactor_scaling", type=float, default=0.0,
                    help="Linear scaling of the per-atom B-factors in the PDB (0 = ignore them; "
                         "1 = use them as is, the Zottig et al. eLife 110114 5LZV recipe with --bfactor 0)")
args = parser.parse_args()

# Build output path: {pdb_stem}[_pdbB{scaling}]_bf{bfactor}_{output_size}_{pixel_size}Apix.mrc
pdb_dir = os.path.dirname(args.pdb)
pdb_stem = os.path.splitext(os.path.basename(args.pdb))[0]
scaling_tag = f"_pdbB{args.pdb_bfactor_scaling:g}" if args.pdb_bfactor_scaling else ""
output_filename = os.path.join(
    pdb_dir,
    f"{pdb_stem}{scaling_tag}_bf{int(args.bfactor)}_{args.output_size}_{args.pixel_size}Apix.mrc"
)

para = SimulateParameters2(
    output_filename=output_filename,
    output_size=args.output_size,
    input_pdb_file=args.pdb,
    pixel_size=args.pixel_size,
    linear_scaling_of_PDB_bfactors=args.pdb_bfactor_scaling,
    per_atom_bfactor=args.bfactor,
)

# UMass used a je_add_metal_ions_to_simulator build; this is the cisTEM 2.0.0-alpha-249 static build installed on cedar.
SIMULATE = os.environ.get("CISTEM_SIMULATE", "/data/software/cistem/part07-2d-template-matching/cistem/simulate")

proc = Popen(
    f"{SIMULATE} --only-modify-signal-3d=2 --wgt=0.225 --water-shell-only",
    shell=True,
    stdin=PIPE,
)
proc.communicate(input=str.encode('\n'.join(map(str, list(para)))))
