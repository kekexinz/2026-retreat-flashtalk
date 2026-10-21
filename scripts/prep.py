"""Preprocess inputs for the 2DTM Blender movie.

Inputs : micrograph (.mrc), template map (.mrc), detections (.star, cisTEM columns)
Outputs (in --out):
  micrograph.png   display image (binned, low-pass filtered, contrast stretched)
  template.ply     iso-surface mesh of the template, centred, units = Angstrom
  scene.json       micrograph geometry + detections (pixel xy, rotation matrix)
  check.png        sanity check: CTF-filtered projections vs micrograph patches

Orientation convention (verified by cross-correlation on EMPIAR-10998):
  pixel x = OriginalXPosition / apix, pixel row y = OriginalYPosition / apix
  (row 0 = first row of the MRC), and the template, in (x=col, y=row, z) coords,
  is rotated by R = (Rz(phi) Ry(theta) Rz(psi))^T and projected along z.
"""
import argparse
import json
import os
import sys

import mrcfile
import numpy as np
from scipy.ndimage import gaussian_filter, zoom
from PIL import Image

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from common import read_cistem_star, euler_zyz, project, ctf2d  # noqa: E402


def detection_matrix(phi, theta, psi):
    return euler_zyz(psi, theta, phi).T


def bin_nd(a, b):
    if b == 1:
        return a
    sl = tuple(slice(0, (s // b) * b) for s in a.shape)
    a = a[sl]
    shp = []
    for s in a.shape:
        shp += [s // b, b]
    return a.reshape(shp).mean(axis=tuple(range(1, 2 * a.ndim, 2)))


def make_micrograph_png(mic, apix, out, width, lowpass_A):
    b = max(1, int(round(mic.shape[1] / width)))
    m = bin_nd(mic, b)
    m = gaussian_filter(m, lowpass_A / (apix * b) / 2.355)
    # flatten slow background (ice gradient) so contrast stays even
    m = m - gaussian_filter(m, 600 / (apix * b))
    lo, hi = np.percentile(m, [0.5, 99.5])
    img = np.clip((m - lo) / (hi - lo), 0, 1)
    # MRC row 0 is drawn at the bottom of the plane in Blender (y up), so flip
    # rows for the PNG whose first row is the top of the image.
    Image.fromarray((img[::-1] * 255).astype(np.uint8)).save(out)
    return b


def make_mesh(vol, apix, out, bin_factor, lowpass_A, level, mass_kda=None):
    from skimage.measure import marching_cubes
    n0 = vol.shape[0]
    v = bin_nd(vol, bin_factor)
    a = apix * bin_factor
    v = gaussian_filter(v, lowpass_A / a / 2.355)
    if level is None and mass_kda:
        # threshold enclosing the molecular volume (1.21 A^3 per Da)
        nvox = mass_kda * 1e3 * 1.21 / a ** 3
        level = float(np.sort(v, axis=None)[-int(nvox)])
    elif level is None:
        level = float(np.percentile(v, 96))
    verts, faces, normals, _ = marching_cubes(v, level=level)
    # verts are (z, y, x) in binned voxels -> (x, y, z) Angstrom about centre
    c = (n0 - 1) / 2.0
    xyz = verts[:, ::-1] * bin_factor + (bin_factor - 1) / 2.0 - c
    xyz *= apix
    # (z,y,x) -> (x,y,z) is a mirror swap, which already turns the winding of
    # marching_cubes faces outward for Blender; do not flip again.
    with open(out, 'w') as f:
        f.write('ply\nformat ascii 1.0\n')
        f.write(f'element vertex {len(xyz)}\nproperty float x\nproperty float y\nproperty float z\n')
        f.write(f'element face {len(faces)}\nproperty list uchar int vertex_indices\nend_header\n')
        np.savetxt(f, xyz, fmt='%.3f')
        np.savetxt(f, np.c_[np.full(len(faces), 3), faces], fmt='%d')
    return level, len(xyz), len(faces)


def resample(a, f):
    """Scale array sampling by factor f (<1 shrinks); integer bins when possible."""
    if abs(f - 1) < 1e-3:
        return a
    inv = 1 / f
    if abs(inv - round(inv)) < 1e-3:
        return bin_nd(a, int(round(inv)))
    return zoom(a, f, order=1)


def check_detections(vol, tapix, mic, apix, dets, out, k=8):
    """Correlate CTF-filtered template projections with the micrograph at the
    first k detections, for the assumed convention and for the alternatives.
    A clearly higher SNR for 'used' confirms positions and angles line up."""
    w = max(tapix, 2 * apix)
    v = resample(vol, tapix / w)
    m = resample(mic, apix / w)
    n = v.shape[0] - v.shape[0] % 2
    v = v[:n, :n, :n]
    norm = lambda z: np.clip((z - z.mean()) / (3 * z.std()) * 0.5 + 0.5, 0, 1)
    alts = {'used': lambda p, t, s: euler_zyz(s, t, p).T,
            'inverse': lambda p, t, s: euler_zyz(s, t, p),
            'swapped': lambda p, t, s: euler_zyz(p, t, s).T,
            'y-flipped': None}
    snr = {key: [] for key in alts}
    tiles = []
    for d in dets[:k]:
        for key, fn in alts.items():
            y = (mic.shape[0] - 1 - d['y']) if key == 'y-flipped' else d['y']
            cx, cy = int(round(d['x'] * apix / w)), int(round(y * apix / w))
            if not (n // 2 <= cx < m.shape[1] - n // 2 and n // 2 <= cy < m.shape[0] - n // 2):
                continue
            patch = m[cy - n // 2:cy + n // 2, cx - n // 2:cx + n // 2]
            patch = patch - patch.mean()
            R = np.array(d['R']) if fn is None else fn(d['phi'], d['theta'], d['psi'])
            proj = project(v, R)
            t = np.real(np.fft.ifft2(np.fft.fft2(proj) * ctf2d(n, w, d['df1'], d['df2'], d['ast'])))
            t -= t.mean()
            t /= np.linalg.norm(t)
            cc = np.real(np.fft.ifft2(np.fft.fft2(patch) * np.conj(np.fft.fft2(t))))
            cc = np.fft.fftshift(cc)[n // 2 - 4:n // 2 + 5, n // 2 - 4:n // 2 + 5]
            snr[key].append(cc.max() / patch.std())
            if key == 'used':
                tiles.append(np.hstack([norm(gaussian_filter(patch, 1.0)),
                                        norm(gaussian_filter(t, 1.0)), norm(proj)]))
    for key, vals in snr.items():
        print(f'  check {key:10s} mean SNR {np.mean(vals):6.2f}  ', np.round(vals, 1))
    # columns: micrograph patch | CTF-filtered projection | raw projection
    if tiles:
        Image.fromarray((np.vstack(tiles) * 255).astype(np.uint8)).save(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--micrograph', required=True)
    ap.add_argument('--template', required=True)
    ap.add_argument('--star', required=True)
    ap.add_argument('--match', default=None,
                    help='substring selecting star rows for this micrograph '
                         '(default: micrograph file stem; use "" for all rows)')
    ap.add_argument('--apix', type=float, default=None, help='default: MRC header')
    ap.add_argument('--min-score', type=float, default=-1e9)
    ap.add_argument('--max-detections', type=int, default=0)
    ap.add_argument('--png-width', type=int, default=2880)
    ap.add_argument('--png-lowpass', type=float, default=12.0, help='Angstrom FWHM')
    ap.add_argument('--template-apix', type=float, default=None, help='default: MRC header')
    ap.add_argument('--defocus', type=float, default=10000.0,
                    help='defocus (A) for the check when the star has no defocus columns')
    ap.add_argument('--mesh-bin', type=int, default=None, help='default: 1 (full sampling)')
    ap.add_argument('--mesh-lowpass', type=float, default=4.0, help='Angstrom FWHM')
    ap.add_argument('--iso', type=float, default=None, help='iso level (binned, filtered map)')
    ap.add_argument('--mass-kda', type=float, default=None,
                    help='set iso level to enclose this mass (e.g. 2200 for a 60S)')
    ap.add_argument('--no-check', action='store_true')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)

    with mrcfile.open(a.micrograph, permissive=True) as f:
        mic = np.asarray(f.data, dtype=np.float32)
        apix = a.apix or float(f.voxel_size.x)
    if mic.ndim == 3:
        mic = mic[0]
    with mrcfile.open(a.template, permissive=True) as f:
        vol = np.asarray(f.data, dtype=np.float32)
        tapix = a.template_apix or float(f.voxel_size.x) or apix

    match = os.path.splitext(os.path.basename(a.micrograph))[0] if a.match is None else a.match
    rows = read_cistem_star(a.star)
    for k, r in enumerate(rows):
        r['_row'] = k   # file order; autocorrMatchedReferenceIndex refers to it
    rows = [r for r in rows if match in str(r.get('cisTEMOriginalImageFilename', ''))]
    rows = [r for r in rows if r.get('cisTEMScore', 0) >= a.min_score]
    rows.sort(key=lambda r: -r.get('cisTEMScore', 0))
    if a.max_detections:
        rows = rows[:a.max_detections]
    dets = []
    for r in rows:
        phi, th, psi = r['cisTEMAnglePhi'], r['cisTEMAngleTheta'], r['cisTEMAnglePsi']
        dets.append(dict(
            x=r['cisTEMOriginalXPosition'] / apix, y=r['cisTEMOriginalYPosition'] / apix,
            phi=phi, theta=th, psi=psi, score=r.get('cisTEMScore', 0.0),
            df1=r.get('cisTEMDefocus1', a.defocus), df2=r.get('cisTEMDefocus2', a.defocus),
            ast=r.get('cisTEMDefocusAngle', 0.0),
            matched=int(r.get('autocorrMatchedReferenceIndex', 0)) >= 0,
            ref_index=int(r.get('autocorrMatchedReferenceIndex', -1)),
            star_row=r['_row'],
            R=detection_matrix(phi, th, psi).tolist()))
    print(f'{len(dets)} detections selected (match="{match}")')

    b = make_micrograph_png(mic, apix, os.path.join(a.out, 'micrograph.png'), a.png_width, a.png_lowpass)
    mesh_bin = a.mesh_bin or 1
    level, nv, nf = make_mesh(vol, tapix, os.path.join(a.out, 'template.ply'), mesh_bin, a.mesh_lowpass, a.iso, a.mass_kda)
    print(f'mesh: iso={level:.4g} verts={nv} faces={nf}')
    scene = dict(apix=apix, mic_w=mic.shape[1], mic_h=mic.shape[0], template_box=vol.shape[0] * tapix / apix,
                 png_bin=b, template_path=os.path.abspath(a.template), template_apix=tapix,
                 iso=level, mesh_bin=mesh_bin, mesh_lowpass=a.mesh_lowpass, detections=dets)
    with open(os.path.join(a.out, 'scene.json'), 'w') as f:
        json.dump(scene, f, indent=1)
    if not a.no_check:
        check_detections(vol, tapix, mic, apix, dets, os.path.join(a.out, 'check.png'))
    print('wrote', a.out)


if __name__ == '__main__':
    main()
