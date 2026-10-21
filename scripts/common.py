"""Shared helpers: cisTEM star parsing, Euler rotations, projections, CTF."""
import numpy as np
import shlex


def read_cistem_star(path):
    """Return list of dicts from a cisTEM/RELION-style single loop star file."""
    cols, rows = [], []
    with open(path) as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith('#') or s.startswith('data_') or s == 'loop_':
                continue
            if s.startswith('_'):
                cols.append(s.split()[0][1:])
                continue
            vals = shlex.split(s)
            if len(vals) != len(cols):
                continue
            row = {}
            for c, v in zip(cols, vals):
                try:
                    row[c] = float(v)
                except ValueError:
                    row[c] = v
            rows.append(row)
    return rows


def rot_z(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def rot_y(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def euler_zyz(phi, theta, psi):
    """Active rotation R = Rz(psi) Ry(theta) Rz(phi), angles in degrees."""
    p, t, s = np.radians([phi, theta, psi])
    return rot_z(s) @ rot_y(t) @ rot_z(p)


def project(vol, R, order=1):
    """Rotate vol (z,y,x) actively by R acting on (x,y,z) coords, sum along z.

    Returns image indexed [y, x].
    """
    from scipy.ndimage import affine_transform
    n = vol.shape[0]
    c = (n - 1) / 2.0
    # output coord o (xyz) samples input at R^T o. Arrays are (z,y,x) -> permute.
    P = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0]])  # xyz <-> zyx
    M = P @ R.T @ P
    off = c - M @ np.full(3, c)
    rv = affine_transform(vol, M, offset=off, order=order)
    return rv.sum(axis=0)


def ctf2d(n, apix, df1, df2, ast_deg, kv=300.0, cs_mm=2.7, ac=0.07):
    lam = 12.2643247 / np.sqrt(kv * 1e3 * (1 + kv * 1e3 * 0.978466e-6))
    f = np.fft.fftfreq(n, d=apix)
    fx, fy = np.meshgrid(f, f)
    k2 = fx ** 2 + fy ** 2
    ang = np.arctan2(fy, fx)
    df = 0.5 * (df1 + df2 + (df1 - df2) * np.cos(2 * (ang - np.radians(ast_deg))))
    cs = cs_mm * 1e7
    chi = np.pi * lam * df * k2 - 0.5 * np.pi * cs * lam ** 3 * k2 ** 2
    return -(np.sqrt(1 - ac ** 2) * np.sin(chi) + ac * np.cos(chi))
