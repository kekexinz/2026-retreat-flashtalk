#!/usr/bin/env bash
# Website version of the 2DTM movie: no text (the page adds its own labels),
# transparent background, cropped to the content, rendered at 2x and
# downscaled so the micrograph and maps stay sharp. --gpu-backend vulkan puts
# EEVEE on the RTX 5080; the OpenGL default lands on the AMD iGPU (~15x slower).
#   web/2dtm_lysate_light.mp4 / _dark.mp4 : H.264 on the site's light (#ffffff)
#       and dark (#1c1c1d) page colours; plays everywhere incl. iPhone/Safari.
#       Dark uses #1d1d1e so it decodes to within 1 level of #1c1c1d.
#   web/2dtm_lysate.webm : VP9 with alpha (Chrome, Firefox, Edge only)
#   web/2dtm_lysate_poster.png : first frame with content, transparent
set -euo pipefail
cd "$(dirname "$0")"

BLENDER=${BLENDER:-/data/software/blender-5.2.2-linux-x64/blender}
W=work/lysate_12781
RES=${RES:-3840x2160}
WIDTH=${WIDTH:-1800}
OUT=web
mkdir -p $OUT

[ -f $W/scene.json ] || scripts/prepare_inputs.sh

if [ "${SKIP_RENDER:-0}" != 1 ]; then
  rm -rf $W/frames_web
  "$BLENDER" -b --gpu-backend vulkan -P scripts/make_movie.py -- --work $W --mn --no-text --transparent \
    --res $RES --border 0.026,0.974,0.065,0.894 --out $W/frames_web --render
fi

if [ -n "$WIDTH" ]; then SCALE="scale=${WIDTH}:-2:flags=lanczos"; else SCALE="crop=trunc(iw/2)*2:trunc(ih/2)*2"; fi
ffmpeg -v error -y -framerate 24 -i $W/frames_web/f_%04d.png -vf "$SCALE,format=yuva420p" \
  -c:v libvpx-vp9 -pix_fmt yuva420p -crf 30 -b:v 0 -row-mt 1 -deadline good -cpu-used 1 \
  -auto-alt-ref 0 -an $OUT/2dtm_lysate.webm
# overlay in packed RGB: the planar (gbrp) colour source leaves a black
# 8 px stripe at the right edge at this width
for v in light:0xffffff dark:0x1d1d1e; do
  ffmpeg -v error -y -framerate 24 -i $W/frames_web/f_%04d.png -f lavfi -i color=${v#*:}:s=1800x886:r=24 \
    -filter_complex "[0]$SCALE,format=rgba[v];[1]format=rgba[bg];[bg][v]overlay=shortest=1:format=rgb,format=rgb24,scale=out_color_matrix=bt709:out_range=tv:flags=accurate_rnd+full_chroma_int,format=yuv420p" \
    -c:v libx264 -profile:v high -preset slow -crf 22 -color_range tv \
    -bsf:v h264_metadata=colour_primaries=1:transfer_characteristics=1:matrix_coefficients=1 \
    -movflags +faststart -an $OUT/2dtm_lysate_${v%%:*}.mp4
done
ffmpeg -v error -y -i $W/frames_web/f_0200.png -vf "$SCALE" -frames:v 1 $OUT/2dtm_lysate_poster.png
ls -la $OUT
