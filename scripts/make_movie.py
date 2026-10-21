"""Build (and optionally render) the 2DTM sliding-window movie in Blender.

  blender -b -P scripts/make_movie.py -- --work work/placeholder [--render] [--save scene.blend]

Reads micrograph.png, template.ply and scene.json written by prep.py.
Layout (orthographic top view, 1 Blender unit = 100 A):
  left : the 3D template, tumbling through orientations
  right: the micrograph, flat, seen from the top
A window the size of the template box slides over the micrograph in a
serpentine raster. When it reaches a detection from the star file, a copy of
the template is placed there in the detected orientation (red overlay).
The first few detections are shown slowly: the window locks on and the
template turns into the detected orientation before the overlay drops in.
"""
import argparse
import json
import math
import os
import random
import sys

import bpy
from mathutils import Matrix, Quaternion, Vector

# ----------------------------------------------------------------------------
# settings
argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
ap = argparse.ArgumentParser()
ap.add_argument('--work', required=True)
ap.add_argument('--out', default=None, help='frame output dir (default: <work>/frames)')
ap.add_argument('--render', action='store_true')
ap.add_argument('--save', default=None, help='save .blend to this path')
ap.add_argument('--fps', type=int, default=24)
ap.add_argument('--res', default='1920x1080')
ap.add_argument('--samples', type=int, default=32)
ap.add_argument('--featured', type=int, default=2, help='detections shown slowly')
ap.add_argument('--intro', type=float, default=1.0, help='seconds before the scan starts')
ap.add_argument('--outro', type=float, default=2.0, help='seconds held after the scan')
ap.add_argument('--row-gap', type=float, default=0.2, help='seconds to move to the next row')
ap.add_argument('--lock-seconds', default='0.25,0.6,0.2,0.5,0.25',
                help='featured lock-on: move in, turn, pop, hold, move out')
ap.add_argument('--row-seconds', default='2,1.2,1,0.8',
                help='sweep duration per row; last value repeats')
ap.add_argument('--window-px', type=float, default=None,
                help='sliding window size in micrograph pixels (default: template box)')
ap.add_argument('--frames', default=None, help='render only a:b (for previews)')
ap.add_argument('--step', type=int, default=1, help='render every n-th frame')
ap.add_argument('--mn', action='store_true',
                help='render the template with Molecular Nodes (density surface + AO material)')
ap.add_argument('--mn-ao-distance', type=float, default=25.0, help='AO sampling distance (A)')
ap.add_argument('--mn-ao-exponent', type=float, default=3.0)
ap.add_argument('--mn-threshold', type=float, default=None, help='default: iso from prep.py')
ap.add_argument('--row0', choices=['top', 'bottom'], default='bottom',
                help='where MRC row 0 is drawn: top (numpy/imshow, PNG previews) or '
                     'bottom (IMOD/cisTEM). "top" is the same scene turned 180 deg about x, '
                     'so molecules keep their handedness.')
ap.add_argument('--no-static-text', action='store_true',
                help='leave out the fixed panel titles (drawn by a sharp overlay instead)')
ap.add_argument('--static-text-png', default=None,
                help='render only the fixed panel titles to this transparent PNG and exit')
ap.add_argument('--text-scale', type=float, default=2.0, help='resolution multiplier for --static-text-png')
ap.add_argument('--no-text', action='store_true',
                help='no text at all (titles, status, counter); e.g. when a web page adds its own')
ap.add_argument('--transparent', action='store_true',
                help='transparent background (RGBA frames, e.g. for a website)')
ap.add_argument('--border', default=None,
                help='render only this region, cropped: xmin,xmax,ymin,ymax as frame fractions (y from bottom)')
ap.add_argument('--seed', type=int, default=1)
A = ap.parse_args(argv)

WORK = os.path.abspath(A.work)
S = json.load(open(os.path.join(WORK, 'scene.json')))
FPS = A.fps
RES_X, RES_Y = map(int, A.res.split('x'))
random.seed(A.seed)

U = 0.01                     # Blender units per Angstrom
APIX = S['apix']
MIC_W = S['mic_w'] * APIX * U
MIC_H = S['mic_h'] * APIX * U
WIN_PX = A.window_px or S['template_box']
WIN = WIN_PX * APIX * U

# colours (linear-ish RGB)
BG = (0.012, 0.014, 0.02)
C_TEMPLATE = (0.40, 0.45, 0.95)   # blue/purple
C_HIT = (0.90, 0.12, 0.08)        # red
C_WINDOW = (1.0, 0.78, 0.18)
C_TEXT = (0.92, 0.92, 0.92)

# frame geometry: micrograph fills ~82% of the height, on the right
FRAME_H = MIC_H / 0.80
FRAME_W = FRAME_H * RES_X / RES_Y
MARGIN = 0.035 * FRAME_W
MIC_X1 = FRAME_W / 2 - MARGIN
MIC_X0 = MIC_X1 - MIC_W
MIC_Y0 = -MIC_H / 2 - 0.02 * FRAME_H
LEFT_CX = (-FRAME_W / 2 + MIC_X0) / 2
LEFT_W = MIC_X0 - (-FRAME_W / 2)


ROW0_TOP = A.row0 == 'top'
# turning the scene 180 deg about x (row 0 at the top) = viewing the specimen
# from the other side: positions get y -> -y, orientations get F @ R
FLIP = Matrix.Rotation(math.pi, 3, 'X') if ROW0_TOP else Matrix.Identity(3)


def mic_to_world(x, y, z=0.0):
    if ROW0_TOP:
        y = S['mic_h'] - y
    return Vector((MIC_X0 + x * APIX * U, MIC_Y0 + y * APIX * U, z))


def det_quat(d):
    return (FLIP @ Matrix(d['R'])).to_quaternion()


# ----------------------------------------------------------------------------
# scene reset
bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.render.fps = FPS
scene.render.resolution_x, scene.render.resolution_y = RES_X, RES_Y
scene.render.resolution_percentage = 100
for eng in ('BLENDER_EEVEE', 'BLENDER_EEVEE_NEXT'):
    try:
        scene.render.engine = eng
        break
    except TypeError:
        pass
try:
    scene.eevee.taa_render_samples = A.samples
except AttributeError:
    pass
for attr, val in (('use_raytracing', True), ('use_shadows', True), ('use_gtao', True)):
    if hasattr(scene.eevee, attr):
        setattr(scene.eevee, attr, val)
scene.view_settings.view_transform = 'Standard'
scene.view_settings.look = 'None'

world = bpy.data.worlds.new('World')
scene.world = world
world.use_nodes = True
# camera rays see the dark background; objects are lit by a bright, even
# ambient (with AO this gives the soft ChimeraX 'full lighting' look)
wnt = world.node_tree
bg_cam = wnt.nodes['Background']
bg_cam.inputs['Color'].default_value = (*BG, 1)
bg_amb = wnt.nodes.new('ShaderNodeBackground')
bg_amb.inputs['Color'].default_value = (1, 1, 1, 1)
bg_amb.inputs['Strength'].default_value = 0.25
lp = wnt.nodes.new('ShaderNodeLightPath')
wmix = wnt.nodes.new('ShaderNodeMixShader')
wnt.links.new(lp.outputs['Is Camera Ray'], wmix.inputs['Fac'])
wnt.links.new(bg_amb.outputs[0], wmix.inputs[1])
wnt.links.new(bg_cam.outputs[0], wmix.inputs[2])
wnt.links.new(wmix.outputs[0], wnt.nodes['World Output'].inputs['Surface'])

cam_data = bpy.data.cameras.new('Camera')
cam_data.type = 'ORTHO'
cam_data.ortho_scale = FRAME_W
cam_data.clip_start, cam_data.clip_end = 0.1, 1000
cam = bpy.data.objects.new('Camera', cam_data)
cam.location = (0, 0, 200)
scene.collection.objects.link(cam)
scene.camera = cam

# The camera looks straight down, so the key light comes in low from the
# upper left to give the surface relief; energies stay low enough that the
# 'Standard' view transform does not clip the colours flat.
def sun_light(name, energy, rx, ry, angle=3.0):
    d = bpy.data.lights.new(name, 'SUN')
    d.energy = energy
    d.angle = math.radians(angle)
    ob = bpy.data.objects.new(name, d)
    ob.rotation_euler = (math.radians(rx), math.radians(ry), 0)
    scene.collection.objects.link(ob)
    return ob


sun_light('Key', 3.2, -45, -40, angle=6)
sun_light('Fill', 0.5, 30, 35, angle=30)


# ----------------------------------------------------------------------------
# materials
def emission_mat(name, color, strength=1.0, alpha=1.0):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    em = nt.nodes.new('ShaderNodeEmission')
    em.inputs['Color'].default_value = (*color, 1)
    em.inputs['Strength'].default_value = strength
    if alpha < 1.0:
        tr = nt.nodes.new('ShaderNodeBsdfTransparent')
        mix = nt.nodes.new('ShaderNodeMixShader')
        mix.inputs['Fac'].default_value = alpha
        nt.links.new(tr.outputs[0], mix.inputs[1])
        nt.links.new(em.outputs[0], mix.inputs[2])
        nt.links.new(mix.outputs[0], out.inputs['Surface'])
        set_blended(m)
    else:
        nt.links.new(em.outputs[0], out.inputs['Surface'])
    return m


def set_blended(m):
    for attr, val in (('surface_render_method', 'BLENDED'), ('blend_method', 'BLEND')):
        try:
            setattr(m, attr, val)
        except (AttributeError, TypeError):
            pass
    for attr in ('use_transparency_overlap', 'show_transparent_back'):
        if hasattr(m, attr):
            setattr(m, attr, False)


def surface_mat(name, color, ao_dist, alpha=1.0, emit=0.0):
    m = bpy.data.materials.new(name)
    m.use_nodes = True
    nt = m.node_tree
    b = nt.nodes['Principled BSDF']
    # ambient occlusion darkens grooves and crevices
    ao = nt.nodes.new('ShaderNodeAmbientOcclusion')
    ao.inputs['Distance'].default_value = ao_dist
    ao.samples = 16
    ao.inputs['Color'].default_value = (*color, 1)
    nt.links.new(ao.outputs['Color'], b.inputs['Base Color'])
    b.inputs['Roughness'].default_value = 0.7
    for nm in ('Specular IOR Level', 'Specular'):
        if nm in b.inputs:
            b.inputs[nm].default_value = 0.2
    b.inputs['Emission Color'].default_value = (*color, 1)
    b.inputs['Emission Strength'].default_value = emit
    b.inputs['Alpha'].default_value = alpha
    if alpha < 1.0:
        set_blended(m)
    m.use_backface_culling = True
    return m


# ----------------------------------------------------------------------------
# micrograph plane
img = bpy.data.images.load(os.path.join(WORK, 'micrograph.png'))
img.colorspace_settings.name = 'sRGB'
mic_mat = bpy.data.materials.new('Micrograph')
mic_mat.use_nodes = True
nt = mic_mat.node_tree
nt.nodes.clear()
tex = nt.nodes.new('ShaderNodeTexImage')
tex.image = img
tex.interpolation = 'Cubic'
em = nt.nodes.new('ShaderNodeEmission')
out = nt.nodes.new('ShaderNodeOutputMaterial')
nt.links.new(tex.outputs['Color'], em.inputs['Color'])
nt.links.new(em.outputs[0], out.inputs['Surface'])

me = bpy.data.meshes.new('MicrographPlane')
me.from_pydata([(0, 0, 0), (MIC_W, 0, 0), (MIC_W, MIC_H, 0), (0, MIC_H, 0)], [], [(0, 1, 2, 3)])
uv = me.uv_layers.new()
# micrograph.png has MRC row 0 as its last (bottom) row
uvs = [(0, 1), (1, 1), (1, 0), (0, 0)] if ROW0_TOP else [(0, 0), (1, 0), (1, 1), (0, 1)]
for li, co in zip(range(4), uvs):
    uv.data[li].uv = co
me.materials.append(mic_mat)
mic = bpy.data.objects.new('Micrograph', me)
mic.location = (MIC_X0, MIC_Y0, 0)
scene.collection.objects.link(mic)

# ----------------------------------------------------------------------------
# template: either the prep.py mesh (units Angstrom), or Molecular Nodes
# density objects (units nm, i.e. 0.1 Blender unit per Angstrom). Objects are
# scaled so that `scale` below is always Blender units per Angstrom.
KINDS = {'template': C_TEMPLATE, 'window': C_TEMPLATE, 'hit': C_HIT}
if A.mn:
    import importlib
    bpy.ops.preferences.addon_enable(module='bl_ext.blender_org.molecularnodes')
    mn = importlib.import_module('bl_ext.blender_org.molecularnodes')
    Grids = importlib.import_module('bl_ext.blender_org.molecularnodes.entities.density.grids').Grids
    MN_UNIT = 0.1    # MN Blender units per Angstrom
    dens = Grids.load(S['template_path'], style='density_surface', center=True)
    mn_src = dens.object
    mn_tree = mn_src.modifiers[0].node_group
    mn_tree.nodes['Density Style Surface'].inputs['Threshold'].default_value = \
        A.mn_threshold if A.mn_threshold is not None else S['iso']
    bpy.context.view_layer.update()
    mn_eval = mn_src.evaluated_get(bpy.context.evaluated_depsgraph_get())
    mn_geo = mn_eval.evaluated_geometry()
    geo = mn_geo.mesh
    n_faces = len(geo.polygons)
    radius = max(v.co.length for v in geo.vertices) / MN_UNIT
    print(f'MN surface: {n_faces} faces, radius {radius:.0f} A')
    # bake the MN surface once; every copy shares a mesh instead of
    # re-running the density style on the full volume
    import numpy as np
    co = np.empty(len(geo.vertices) * 3, dtype=np.float32)
    geo.vertices.foreach_get('co', co)
    loop_vert = np.empty(len(geo.loops), dtype=np.int32)
    geo.loops.foreach_get('vertex_index', loop_vert)
    loop_start = np.empty(len(geo.polygons), dtype=np.int32)
    geo.polygons.foreach_get('loop_start', loop_start)
    mn_mesh = bpy.data.meshes.new('MN_Surface')
    mn_mesh.vertices.add(len(geo.vertices))
    mn_mesh.vertices.foreach_set('co', co)
    mn_mesh.loops.add(len(loop_vert))
    mn_mesh.loops.foreach_set('vertex_index', loop_vert)
    mn_mesh.polygons.add(len(loop_start))
    mn_mesh.polygons.foreach_set('loop_start', loop_start)
    mn_mesh.update(calc_edges=True)
    mn_mesh.validate()
    mn_mesh.shade_smooth()
    del mn_eval, mn_geo, geo
    mn_src.hide_render = True
    mn_src.hide_set(True)
else:
    MN_UNIT = 1.0
    bpy.ops.wm.ply_import(filepath=os.path.join(WORK, 'template.ply'))
    tmpl_src = bpy.context.selected_objects[0]
    tmpl_mesh = tmpl_src.data
    tmpl_mesh.name = 'TemplateMesh'
    for p in tmpl_mesh.polygons:
        p.use_smooth = True
    bpy.data.objects.remove(tmpl_src)
    radius = max(v.co.length for v in tmpl_mesh.vertices)

LEFT_SCALE = 0.80 * min(LEFT_W, FRAME_H * 0.6) / (2 * radius)
KIND_SCALE = {'template': LEFT_SCALE, 'window': U, 'hit': U}


MN_HIT_FACES = 250000
MN_KIND_MESH = {}


def mn_kind_mesh(kind):
    """Baked MN surface for this kind: decimated for the small copies, with the
    kind's colour as the 'Color' attribute and an MN ambient-occlusion
    material whose sampling distance matches the object's size in the scene."""
    if kind in MN_KIND_MESH:
        return MN_KIND_MESH[kind]
    me = mn_mesh.copy()
    me.name = 'MN_' + kind
    if kind != 'template' and len(me.polygons) > MN_HIT_FACES:
        tmp = bpy.data.objects.new('tmp', me)
        scene.collection.objects.link(tmp)
        dec = tmp.modifiers.new('dec', 'DECIMATE')
        dec.ratio = MN_HIT_FACES / len(me.polygons)
        me = bpy.data.meshes.new_from_object(tmp.evaluated_get(bpy.context.evaluated_depsgraph_get()))
        bpy.data.objects.remove(tmp)
        me.name = 'MN_' + kind
        me.shade_smooth()
    col = me.color_attributes.new('Color', 'FLOAT_COLOR', 'POINT')
    col.data.foreach_set('color', [*KINDS[kind], 1.0] * len(me.vertices))
    mat = mn.material.AmbientOcclusion(name='MN_AO_' + kind)
    mat.distance = A.mn_ao_distance / 10 * KIND_SCALE[kind] / MN_UNIT   # A -> nm, scaled
    mat.exponent = A.mn_ao_exponent
    me.materials.clear()
    me.materials.append(mat.material)
    MN_KIND_MESH[kind] = me
    return me


def template_object(name, kind, mesh=None):
    if A.mn:
        ob = bpy.data.objects.new(name, mn_kind_mesh(kind))
    else:
        mesh = mesh or tmpl_mesh
        if len(mesh.materials) == 0:
            mesh.materials.append(MATS[kind])
        ob = bpy.data.objects.new(name, mesh)
        ob.material_slots[0].link = 'OBJECT'
        ob.material_slots[0].material = MATS[kind]
    scene.collection.objects.link(ob)
    ob.rotation_mode = 'QUATERNION'
    ob.scale = (KIND_SCALE[kind] / MN_UNIT,) * 3
    return ob


if not A.mn:
    hit_mesh = tmpl_mesh
    HIT_FACES = 200000
    if len(tmpl_mesh.polygons) > HIT_FACES:
        tmp = bpy.data.objects.new('tmp', tmpl_mesh.copy())
        scene.collection.objects.link(tmp)
        mod = tmp.modifiers.new('dec', 'DECIMATE')
        mod.ratio = HIT_FACES / len(tmpl_mesh.polygons)
        dg = bpy.context.evaluated_depsgraph_get()
        hit_mesh = bpy.data.meshes.new_from_object(tmp.evaluated_get(dg))
        hit_mesh.name = 'TemplateMeshLow'
        bpy.data.objects.remove(tmp)
    MATS = {'template': surface_mat('Template', C_TEMPLATE, 0.3 * radius * LEFT_SCALE),
            'window': surface_mat('WindowCopy', C_TEMPLATE, 0.3 * radius * U, alpha=0.9),
            'hit': surface_mat('Hit', C_HIT, 0.3 * radius * U, alpha=0.9)}
else:
    hit_mesh = None


big = template_object('Template', 'template')
big.location = (LEFT_CX, -0.02 * FRAME_H, 20)


# ----------------------------------------------------------------------------
# text
def text(name, body, loc, size, color=C_TEXT, align='CENTER'):
    cu = bpy.data.curves.new(name, 'FONT')
    cu.body = body
    cu.size = size
    cu.align_x = align
    cu.align_y = 'CENTER'
    ob = bpy.data.objects.new(name, cu)
    ob.location = loc
    cu.materials.append(emission_mat(name + 'Mat', color, 1.0))
    scene.collection.objects.link(ob)
    return ob


TS = 0.045 * FRAME_H
title_y = MIC_Y0 + MIC_H + 0.055 * FRAME_H
if not A.no_static_text:
    text('LabelTemplate', 'Template (3D map)', (LEFT_CX, title_y, 50), TS)
    text('LabelMicrograph', 'Micrograph', (MIC_X0 + MIC_W / 2, title_y, 50), TS)
status_y = big.location.y - 0.36 * FRAME_H


# ----------------------------------------------------------------------------
# sliding window: square outline + translucent fill + template copy
def square_ring(name, side, width):
    h, g = side / 2, side / 2 - width
    verts = [(-h, -h, 0), (h, -h, 0), (h, h, 0), (-h, h, 0),
             (-g, -g, 0), (g, -g, 0), (g, g, 0), (-g, g, 0)]
    faces = [(0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    m = bpy.data.meshes.new(name)
    m.from_pydata(verts, [], faces)
    return m


win_mat = emission_mat('WindowEdge', C_WINDOW, 2.0)
win = bpy.data.objects.new('Window', square_ring('WindowRing', WIN, WIN * 0.025))
win.data.materials.append(win_mat)
scene.collection.objects.link(win)
fill_me = bpy.data.meshes.new('WindowFill')
h = WIN / 2
fill_me.from_pydata([(-h, -h, 0), (h, -h, 0), (h, h, 0), (-h, h, 0)], [], [(0, 1, 2, 3)])
fill_me.materials.append(emission_mat('WindowFillMat', C_WINDOW, 1.0, alpha=0.10))
win_fill = bpy.data.objects.new('WindowFill', fill_me)
win_fill.parent = win
win_fill.location = (0, 0, -0.5)
scene.collection.objects.link(win_fill)
win_copy = template_object('WindowTemplate', 'window')
win_copy.parent = win
win_copy.location = (0, 0, 2)
con = win_copy.constraints.new('COPY_ROTATION')
con.target = big
win.location = (0, 0, 30)


# ----------------------------------------------------------------------------
# keyframe helpers
def key(ob, path, frame, value=None, interp='BEZIER'):
    if value is not None:
        setattr(ob, path, value)
    bpy.context.preferences.edit.keyframe_new_interpolation_type = interp
    ob.keyframe_insert(path, frame=frame)


def random_quat():
    u1, u2, u3 = random.random(), random.random(), random.random()
    return Quaternion((math.sqrt(1 - u1) * math.sin(2 * math.pi * u2),
                       math.sqrt(1 - u1) * math.cos(2 * math.pi * u2),
                       math.sqrt(u1) * math.sin(2 * math.pi * u3),
                       math.sqrt(u1) * math.cos(2 * math.pi * u3)))


def sec(s):
    return s * FPS


# ----------------------------------------------------------------------------
# plan the scan
dets = S['detections']
n_rows = max(1, math.ceil(S['mic_h'] / WIN_PX))
y_lo, y_hi = WIN_PX / 2, S['mic_h'] - WIN_PX / 2
row_y = [y_hi - i * (y_hi - y_lo) / max(1, n_rows - 1) for i in range(n_rows)]  # top first
x_lo, x_hi = WIN_PX / 2, S['mic_w'] - WIN_PX / 2
row_secs = [float(v) for v in A.row_seconds.split(',')]
for d in dets:
    d['row'] = min(range(n_rows), key=lambda i: abs(row_y[i] - d['y']))

T_INTRO = sec(A.intro)
LOCK_IN, TURN, POP, HOLD, LOCK_OUT = (sec(float(v)) for v in A.lock_seconds.split(','))
win_keys = []          # (frame, x_px, y_px)
pops = []              # (frame, det)
flashes = []           # frames
locks = []             # (start, turned, end, det)
t = T_INTRO
for r in range(n_rows):
    forward = r % 2 == 0
    xa, xb = (x_lo, x_hi) if forward else (x_hi, x_lo)
    speed = (x_hi - x_lo) / sec(row_secs[min(r, len(row_secs) - 1)])
    y = row_y[r]
    if r > 0:
        t += sec(A.row_gap)   # move down to the next row
    win_keys.append((t, xa, y))
    cur = xa
    row_dets = sorted([d for d in dets if d['row'] == r],
                      key=lambda d: d['x'] if forward else -d['x'])
    for d in row_dets:
        xd = min(max(d['x'], x_lo), x_hi)
        t += abs(xd - cur) / speed
        cur = xd
        win_keys.append((t, xd, y))
        if len(locks) < A.featured:
            t0 = t
            win_keys.append((t0 + LOCK_IN, d['x'], d['y']))
            turned = t0 + LOCK_IN + TURN
            pops.append((turned, d))
            flashes.append(turned)
            end = turned + POP + HOLD
            win_keys.append((end, d['x'], d['y']))
            t = end + LOCK_OUT
            win_keys.append((t, xd, y))
            locks.append((t0, turned, end, d))
        else:
            pops.append((t, d))
            flashes.append(t)
    t += abs(xb - cur) / speed
    win_keys.append((t, xb, y))
T_SCAN_END = t
T_END = T_SCAN_END + sec(A.outro)
scene.frame_start = 1
scene.frame_end = int(math.ceil(T_END))
print(f'rows={n_rows} detections={len(dets)} featured={len(locks)} '
      f'frames={scene.frame_end} ({scene.frame_end / FPS:.1f}s)')

# window path
win.scale = (0.001,) * 3
key(win, 'scale', T_INTRO - sec(0.4))
key(win, 'scale', T_INTRO, (1, 1, 1))
for f, x, y in win_keys:
    key(win, 'location', f, mic_to_world(x, y, 30), 'LINEAR')
key(win, 'scale', T_SCAN_END, (1, 1, 1))
key(win, 'scale', T_SCAN_END + sec(0.4), (0.001,) * 3)

# window edge flashes
emis = win_mat.node_tree.nodes['Emission'].inputs['Strength']
emis.default_value = 2.0
emis.keyframe_insert('default_value', frame=1)
for f in flashes:
    for df, v in ((-1, 2.0), (0, 8.0), (sec(0.2), 2.0)):
        emis.default_value = v
        emis.keyframe_insert('default_value', frame=f + df)

# template tumbling, interrupted by the lock-ons
q_prev = Quaternion()


def key_q(f, q):
    global q_prev
    if q.dot(q_prev) < 0:
        q = -q
    q_prev = q
    key(big, 'rotation_quaternion', f, q)


big.scale = (0.01,) * 3
key(big, 'scale', 1)
key(big, 'scale', sec(0.8), (LEFT_SCALE / MN_UNIT,) * 3)
f = 1
step_slow, step_fast = sec(1.0), sec(0.45)
lock_iter = iter(locks)
nxt = next(lock_iter, None)
while f < T_END:
    if nxt and f + (step_fast if f > T_INTRO else step_slow) > nxt[0]:
        t0, turned, end, d = nxt
        key_q(t0, random_quat())
        qd = det_quat(d)
        key_q(turned, qd)
        key_q(end, qd)
        f = end + sec(0.4)
        nxt = next(lock_iter, None)
        continue
    if f >= T_SCAN_END:
        break
    key_q(f, random_quat())
    f += step_fast if f > T_INTRO else step_slow
# finish slowly turning
key_q(T_SCAN_END + sec(1.5), random_quat())
key_q(T_END, random_quat())

# hit overlays
for i, (f, d) in enumerate(pops):
    ob = template_object(f'Hit{i:03d}', 'hit', hit_mesh)
    ob.rotation_quaternion = det_quat(d)
    ob.location = mic_to_world(d['x'], d['y'], 5)
    ob.scale = (1e-4 * U / MN_UNIT,) * 3
    key(ob, 'scale', f)
    key(ob, 'scale', f + sec(0.12), (1.35 * U / MN_UNIT,) * 3)
    key(ob, 'scale', f + sec(0.3), (U / MN_UNIT,) * 3)

# status text and detection counter, updated by a frame handler that is also
# stored in the .blend (text block, registered on load when scripts may run)
status = text('Status', ' ', (LEFT_CX, status_y, 50), TS * 0.8, color=(0.7, 0.7, 0.7))
status_match = text('StatusMatch', ' ', (LEFT_CX, status_y, 50), TS * 0.8, color=C_HIT)
count_y = MIC_Y0 - 0.045 * FRAME_H
counter = text('Counter', ' ', (MIC_X1, count_y, 50), TS * 0.85, align='RIGHT')
timeline = dict(
    intro=T_INTRO, scan_end=T_SCAN_END,
    pops=sorted(f for f, _ in pops),
    locks=[(turned, end + LOCK_OUT, f'match!  SNR {d["score"]:.1f}') for _, turned, end, d in locks])
HANDLER = """import bpy, bisect, json
TL = json.loads(%r)


def movie_text_update(scene, *args):
    f = scene.frame_current
    st, sm, ct = (bpy.data.objects.get(n) for n in ('Status', 'StatusMatch', 'Counter'))
    if not (st and sm and ct):
        return
    match = next((txt for a, b, txt in TL['locks'] if a <= f < b), None)
    sm.data.body = match or ' '
    st.data.body = ' ' if (match or not TL['intro'] <= f < TL['scan_end']) else 'searching orientations\u2026'
    ct.data.body = ('detections: %%d' %% bisect.bisect_right(TL['pops'], f)) if f >= TL['intro'] else ' '


for h in list(bpy.app.handlers.frame_change_pre):
    if getattr(h, '__name__', '') == 'movie_text_update':
        bpy.app.handlers.frame_change_pre.remove(h)
bpy.app.handlers.frame_change_pre.append(movie_text_update)
""" % json.dumps(timeline)
tb = bpy.data.texts.new('movie_handlers.py')
tb.write(HANDLER)
tb.use_module = True
exec(compile(HANDLER, 'movie_handlers.py', 'exec'), {})
scene.render.use_lock_interface = True

# ----------------------------------------------------------------------------
out_dir = os.path.abspath(A.out or os.path.join(WORK, 'frames'))
scene.render.filepath = os.path.join(out_dir, 'f_')
if A.static_text_png:
    # same scene and camera, so the titles land exactly where the movie expects them
    for ob in scene.objects:
        if ob.type != 'CAMERA' and ob.name not in ('LabelTemplate', 'LabelMicrograph'):
            ob.hide_render = True
    scene.render.film_transparent = True
    scene.render.image_settings.color_mode = 'RGBA'
    scene.render.resolution_x = int(RES_X * A.text_scale)
    scene.render.resolution_y = int(RES_Y * A.text_scale)
    scene.frame_set(1)
    scene.render.filepath = os.path.abspath(A.static_text_png)
    bpy.ops.render.render(write_still=True)
    raise SystemExit
scene.render.image_settings.file_format = 'PNG'
if A.no_text:
    for ob in scene.objects:
        if ob.type == 'FONT':
            ob.hide_render = True
if A.transparent:
    scene.render.film_transparent = True
    scene.render.image_settings.color_mode = 'RGBA'
if A.border:
    b = [float(v) for v in A.border.split(',')]
    scene.render.use_border = scene.render.use_crop_to_border = True
    (scene.render.border_min_x, scene.render.border_max_x,
     scene.render.border_min_y, scene.render.border_max_y) = b
if A.save:
    bpy.ops.wm.save_as_mainfile(filepath=os.path.abspath(A.save))
if A.render:
    if A.frames:
        a, b = A.frames.split(':')
        scene.frame_start, scene.frame_end = int(a), min(int(b), scene.frame_end)
    scene.frame_step = A.step
    bpy.ops.render.render(animation=True)
