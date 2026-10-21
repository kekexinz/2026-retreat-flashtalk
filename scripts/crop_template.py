"""Central crop of the simulated 448^3 template to the 320^3 Accuracy box (voxels 64:384, 1.66 A/voxel).

  python scripts/crop_template.py in_448.mrc out_320.mrc
"""
import sys

import mrcfile
import numpy as np

src, dst = sys.argv[1:3]
with mrcfile.open(src, permissive=True) as f:
    v = np.asarray(f.data, dtype=np.float32)
    apix = float(f.voxel_size.x)
if v.shape != (448, 448, 448):
    raise ValueError(f'expected a 448^3 map, got {v.shape}')
with mrcfile.new(dst, overwrite=True) as f:
    f.set_data(np.ascontiguousarray(v[64:384, 64:384, 64:384]))
    f.voxel_size = apix
print('wrote', dst)
