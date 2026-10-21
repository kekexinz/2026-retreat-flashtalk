# Lysate12781 visualization bundle

Selected by metadata from the frozen200,not detection performance:18 reference particles equals the cohort median;among equally close counts,choose mean defocus closest to the cohort median. This frame mean defocus15589.3 Å.

Files:
- micrograph_12781_native_0p83A.mrc:original native image,5760×4092 pixels,0.83 Å/pixel.
- reference_peaks.star:18 deposited refined reference particles.
- accuracy_peaks.star:all19 Accuracy detections,original frozen cutoff7.180006.
- speed_peaks.star:all15 Speed detections,retained cutoff7.144050598144531.
- template_5lzv_60S_accuracy_1p66A_320.mrc:actual320³ full-bandwidth Accuracy template,1.66 Å/voxel. Speed uses the same molecular source with a different cropped/resampled156³ template at2.979487 Å/voxel;the provided map is for the Accuracy template/3D molecular illustration.

Both detection STARs preserve the final Adam poses/positions and all detections,including false positives. Positions are Å in cisTEMOriginalXPosition/YPosition. autocorrOriginalXPixel/YPixel provide native0.83 Å image coordinates,so all3 STARs overlay the same native image. cisTEMPixelSize is0.83;autocorrSearchPixelSize distinguishes the actual search sampling. Image paths point to the image in this bundle. No image rotation/flip or coordinate transformation was applied.

Angles cisTEMAnglePhi/Theta/Psi are degrees,in the pipeline's ZYZ/cisTEM convention. Convert explicitly for Blender's axes/view convention;do not interpret these as BlenderXYZ Euler angles. STAR has no estimated particleZ position from this2D search.

Detection scores are the frozen pre-Adam selection scores;angles/centers are post-Adam. autocorrMatchedReferenceIndex uses the original strict<5 Å/<2° one-to-one assignment,-1 means unmatched under that criterion. Match indices are zero-based rows in reference_peaks.star. The manifest separately provides strict<3 Å/<2° evaluation;matching labels have not been silently changed.

Per-image pipeline timing:Accuracy30.18s,Speed10.40s onRTX5080,excluding prior binning/setup/STAR export. This single frame is an illustration,not a replacement for the200-image pooled statistics. manifest.json records units,source paths,matching statistics and SHA256 checksums.
