"""Render text lines to a transparent PNG with Blender's default font (same as the movie).

  blender -b -P scripts/text_overlay.py -- --spec overlay.json --out overlay.png

spec: {"size": [1920, 1080], "lines": [{"text": "...", "x": 960, "y": 100, "size": 30,
        "color": [0.9, 0.9, 0.9], "align": "CENTER"}, ...]}   (x, y in pixels, y from top)
"""
import argparse
import json
import sys

import bpy

argv = sys.argv[sys.argv.index('--') + 1:] if '--' in sys.argv else []
ap = argparse.ArgumentParser()
ap.add_argument('--spec', required=True)
ap.add_argument('--out', required=True)
ap.add_argument('--scale', type=float, default=1.0, help='output resolution multiplier')
A = ap.parse_args(argv)
spec = json.load(open(A.spec))
W, H = spec.get('size', [1920, 1080])

bpy.ops.wm.read_factory_settings(use_empty=True)
scene = bpy.context.scene
scene.render.engine = 'BLENDER_EEVEE'
scene.eevee.taa_render_samples = 16
scene.render.resolution_x, scene.render.resolution_y = int(W * A.scale), int(H * A.scale)
scene.render.film_transparent = True
scene.render.image_settings.file_format = 'PNG'
scene.render.image_settings.color_mode = 'RGBA'
scene.view_settings.view_transform = 'Standard'

cam_data = bpy.data.cameras.new('Camera')
cam_data.type = 'ORTHO'
cam_data.ortho_scale = W
cam = bpy.data.objects.new('Camera', cam_data)
cam.location = (W / 2, H / 2, 100)
scene.collection.objects.link(cam)
scene.camera = cam

def make_text(name, body, size, color, align, shear=0.0):
    cu = bpy.data.curves.new(name, 'FONT')
    cu.body = body
    cu.size = size
    cu.align_x = align
    cu.align_y = 'CENTER'
    cu.shear = shear
    m = bpy.data.materials.new(name + 'm')
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    em = nt.nodes.new('ShaderNodeEmission')
    em.inputs['Color'].default_value = (*color, 1)
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    nt.links.new(em.outputs[0], out.inputs['Surface'])
    cu.materials.append(m)
    ob = bpy.data.objects.new(name, cu)
    scene.collection.objects.link(ob)
    return ob


# lines with "parts": [{"text", "shear"}] are laid out left to right and centred on x
# (shear ~0.2 gives an italic look; the default font has no italic face)
for i, ln in enumerate([l for l in spec['lines'] if 'parts' in l]):
    obs = [make_text(f'p{i}_{k}', p['text'], ln.get('size', 30), ln.get('color', [0.92] * 3), 'LEFT',
                     p.get('shear', 0.0)) for k, p in enumerate(ln['parts'])]
    bpy.context.view_layer.update()
    gap = ln.get('size', 30) * 0.28
    widths = [o.dimensions.x for o in obs]
    x = ln['x'] - (sum(widths) + gap * (len(obs) - 1)) / 2
    for o, w in zip(obs, widths):
        o.location = (x, H - ln['y'], 0)
        x += w + gap
spec['lines'] = [l for l in spec['lines'] if 'parts' not in l]

for i, ln in enumerate(spec['lines']):
    cu = bpy.data.curves.new(f't{i}', 'FONT')
    cu.body = ln['text']
    cu.size = ln.get('size', 30)
    cu.align_x = ln.get('align', 'CENTER')
    cu.align_y = 'CENTER'
    m = bpy.data.materials.new(f'm{i}')
    m.use_nodes = True
    nt = m.node_tree
    nt.nodes.clear()
    em = nt.nodes.new('ShaderNodeEmission')
    em.inputs['Color'].default_value = (*ln.get('color', [0.92, 0.92, 0.92]), 1)
    out = nt.nodes.new('ShaderNodeOutputMaterial')
    nt.links.new(em.outputs[0], out.inputs['Surface'])
    cu.materials.append(m)
    ob = bpy.data.objects.new(f't{i}', cu)
    ob.location = (ln['x'], H - ln['y'], 0)
    scene.collection.objects.link(ob)

scene.render.filepath = A.out
bpy.ops.render.render(write_still=True)
