#!/usr/bin/env bash
# Talk version of the 2DTM movie (slide 2 of index.html):
#   media/02_2dtm_movie.mp4       1920x1080, 24 fps, without the fixed panel titles
#   media/02_2dtm_movie_text.png  text at 2x (3840x2160), laid over the video by index.html
#                                 so it stays sharp: the panel titles from make_movie.py
#                                 plus the captions in scripts/talk_captions.json
set -euo pipefail
cd "$(dirname "$0")"

BLENDER=${BLENDER:-/data/software/blender-5.2.2-linux-x64/blender}
W=work/lysate_12781

[ -f $W/scene.json ] || scripts/prepare_inputs.sh

rm -rf $W/frames_notext
"$BLENDER" -b -P scripts/make_movie.py -- --work $W --mn --no-static-text --out $W/frames_notext --render
"$BLENDER" -b -P scripts/make_movie.py -- --work $W --mn --static-text-png $W/titles_4k.png --text-scale 2
"$BLENDER" -b -P scripts/text_overlay.py -- --spec scripts/talk_captions.json --out $W/captions_4k.png --scale 2
ffmpeg -v error -y -i $W/titles_4k.png -i $W/captions_4k.png \
  -filter_complex "[0]format=rgba[a];[1]format=rgba[b];[a][b]overlay=format=rgb" -frames:v 1 media/02_2dtm_movie_text.png

# two-pass at 20 Mb/s keeps a ~15 s movie well under 50 MB
for pass in 1 2; do
  ffmpeg -v error -y -framerate 24 -i $W/frames_notext/f_%04d.png -c:v libx264 -pix_fmt yuv420p \
    -preset slow -b:v 20M -maxrate 25M -bufsize 40M -pass $pass -passlogfile $W/x264n \
    -movflags +faststart $( [ $pass = 1 ] && echo "-f mp4 /dev/null" || echo media/02_2dtm_movie.mp4 )
done
echo "wrote media/02_2dtm_movie.mp4 media/02_2dtm_movie_text.png"
