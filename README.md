# 2026 retreat flash talk

Open `index.html` in Chrome (keep the `media/` folder next to it).
F = fullscreen · arrows / space / clicker = next/previous · R = replay movie.
`flashtalk_backup.pdf` is a still-image backup of the same slides.

`web/` has the 2DTM movie for websites: no text, 1800×886, 15 s loop, H.264
(plays everywhere, including Safari and iPhone). `2dtm_lysate_light.mp4` is on
white (#ffffff), `2dtm_lysate_dark.mp4` on #1c1c1d (al-folio dark theme); show
the one matching the theme in `<video autoplay loop muted playsinline>`.
`2dtm_lysate_poster.png` is a transparent still for `poster`.

## Rebuilding the 2DTM movie

The movie shows a 2DTM search window sliding over one BSC-1 lysate micrograph
(EMPIAR-13577, image 12781) and placing the 60S template (PDB 5LZV) at the
19 Accuracy-mode detections.

### Software

| Tool | Version used | Used for |
|---|---|---|
| [Blender](https://www.blender.org/download/) | 5.2.2 LTS | scene and rendering (`scripts/make_movie.py`) |
| [Molecular Nodes](https://extensions.blender.org/add-ons/molecularnodes/) Blender extension | 520.2.0 | template density surface (`--mn`); install from Blender's Get Extensions |
| [cisTEM](https://cistem.org) `simulate` | 2.0.0-alpha-249 (static build) | simulating the template map from the PDB |
| Python | 3.12 | `prep.py`, `template_simulator.py`, `crop_template.py` |
| numpy, scipy, mrcfile, Pillow, scikit-image | numpy 2.5, scipy 1.18, mrcfile 1.5, Pillow 12.3 | micrograph PNG, template mesh, detections |
| ffmpeg | 8.0 (libx264, libvpx-vp9) | frames → MP4 / WebM |
| xz | any | unpacking the micrograph |

```
pip install numpy scipy mrcfile pillow scikit-image
```

Paths default to the lab workstation; override them with environment variables:

| Variable | Default |
|---|---|
| `BLENDER` | `/data/software/blender-5.2.2-linux-x64/blender` |
| `CISTEM_SIMULATE` | `/data/software/cistem/part07-2d-template-matching/cistem/simulate` |
| `PY` | `python` |

`run_web.sh` renders with `--gpu-backend vulkan` so EEVEE runs on the discrete GPU.

### Data (`data/lysate_12781/`)

| File | What it is |
|---|---|
| `micrograph_12781_native_0p83A.mrc.xz` | the micrograph, 5760×4092 at 0.83 Å/px, xz-compressed (lossless; the unpacked `.mrc` matches the SHA256 in `manifest.json`) |
| `5lzv_60S_minCtermL24.pdb` | template model: 5LZV 60S with the L24 C-terminal tail removed (Zottig et al., eLife 110114) |
| `accuracy_peaks.star` | the 19 Accuracy-mode detections shown in the movie (cisTEM columns, post-Adam poses) |
| `speed_peaks.star`, `reference_peaks.star` | Speed-mode detections and the 18 deposited reference particles (not used by the movie) |
| `manifest.json`, `README.md` | provenance, matching statistics and checksums of the original bundle |

The template map is not committed; it is simulated from the PDB (below).

### Simulating the template map

`scripts/template_simulator.py` (from `umass_scripts_archive`) runs cisTEM
`simulate` on the PDB. The Zottig et al. recipe is 448³ at 1.66 Å/px, the
per-atom B-factors from the PDB as is, and no added B-factor:

```
python scripts/template_simulator.py --pdb data/lysate_12781/5lzv_60S_minCtermL24.pdb \
    --output_size 448 --pixel_size 1.66 --bfactor 0 --pdb_bfactor_scaling 1
# -> data/lysate_12781/5lzv_60S_minCtermL24_pdbB1_bf0_448_1.66Apix.mrc

python scripts/crop_template.py data/lysate_12781/5lzv_60S_minCtermL24_pdbB1_bf0_448_1.66Apix.mrc \
    template_5lzv_60S_accuracy_1p66A_320.mrc
```

The output is written next to the PDB. `crop_template.py` takes the central
320³ (the box of the Accuracy search). The result correlates at 0.9999 with
the map used for the talk (Ximena Zottig's original simulation, cropped the
same way), so it is not byte-identical and its SHA256 differs from
`manifest.json`.

`scripts/prepare_inputs.sh` does all of this inside `work/inputs/`, then
runs `scripts/prep.py` to write the micrograph PNG, template mesh and
`scene.json` to `work/lysate_12781/`.

### Rendering

```
./run_talk.sh   # media/02_2dtm_movie.mp4 + media/02_2dtm_movie_text.png (talk)
./run_web.sh    # web/2dtm_lysate_{light,dark}.mp4, web/2dtm_lysate.webm, poster (website)
```

Both call `scripts/prepare_inputs.sh` first if `work/lysate_12781/scene.json`
does not exist. Everything intermediate goes to `work/` (not in git, several GB
once frames are rendered).
