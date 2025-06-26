import os
import sys
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import subprocess

# Step 1: Define paths and input parameters
mean_image_path = sys.argv[1]
processed_image_path = sys.argv[2]
slice_numbers = list(map(int, sys.argv[3:]))  # List of slice numbers as input
output_dir = "overlay_screenshots"
movie_output = f"Signal_Change_Map_Slice{slice_numbers}.mp4"

# Step 2: Load the mean image and processed image
mean_img = nib.load(mean_image_path)
processed_img = nib.load(processed_image_path)

# Checking orientations of the "to be processed" NIFTI files
affine_mean = mean_img.affine
orientation_mean = nib.orientations.aff2axcodes(affine_mean)
print('Orientations of the mean image:', orientation_mean)

affine_processed_img = processed_img.affine
orientation = nib.orientations.aff2axcodes(affine_processed_img)
print('Orientations of the processed image:', orientation)

mean_data = mean_img.get_fdata()  # 3D mean image
processed_data = processed_img.get_fdata()  # 4D processed image

# Ensure the dimensions match (except for the time/volume dimension)
if mean_data.shape != processed_data.shape[:3]:
    raise ValueError("The dimensions of the mean image and processed image do not match!")

# Get the number of volumes in the processed image
num_volumes = processed_data.shape[3]

# Step 3: Process each slice number individually
for slice_number in slice_numbers:
    print(f"Processing slice {slice_number}...")

    # Create a separate output directory for each slice
    slice_output_dir = os.path.join(output_dir, f"slice_{slice_number}")
    os.makedirs(slice_output_dir, exist_ok=True)

    # Define the movie output filename for the current slice
   

    # Step 4: Scroll through each volume and overlay on the mean image
    for vol_idx in range(num_volumes):
        # Extract the current volume of the processed image
        processed_volume = processed_data[..., vol_idx]

        # Overlay: Plot the mean image and processed volume in the same plot
        fig, ax = plt.subplots(figsize=(6, 6))

        # Flip the slices left-right (flip along the horizontal axis, axis 1)
        mean_slice = np.flip(mean_data[..., slice_number], axis=1)
        processed_slice = np.flip(processed_volume[..., slice_number], axis=1)

        # Plot the mean image (underlay)
        ax.imshow(mean_slice, cmap="gray", vmin=np.min(mean_data), vmax=np.max(mean_data))

        # Plot the processed slice overlay
        im_combined = ax.imshow(processed_slice, cmap="hot", alpha=0.4, vmin=4, vmax=20)

        # Add a colorbar for signal change
        cbar_combined = fig.colorbar(im_combined, ax=ax, orientation='vertical', fraction=0.046, pad=0.04)
        cbar_combined.set_label('Percent Signal Change', rotation=270, labelpad=15)

        # Save the screenshot
        screenshot_path = os.path.join(slice_output_dir, f"frame_{vol_idx}.png")
        plt.savefig(screenshot_path, bbox_inches='tight', pad_inches=0)
        plt.close()

        print(f"Captured overlay for volume {vol_idx} in slice {slice_number}")

    # Step 5: Use ffmpeg to combine the frames into an MP4 movie for each slice
    subprocess.run([
        "ffmpeg", "-framerate", "4", "-i", f"{slice_output_dir}/frame_%d.png",
        "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2, transpose=0",  # Rotate 90 degrees anti-clockwise
        "-c:v", "libx264", "-pix_fmt", "yuv420p", movie_output
    ])

    print(f"Movie created and saved as {movie_output} for slice {slice_number}")

    # Optionally, clean up frame images
    # for filename in os.listdir(slice_output_dir):
    #     os.remove(os.path.join(slice_output_dir, filename))
