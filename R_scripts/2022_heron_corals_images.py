library(reticulate)

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import os
from natsort import natsorted

# Define the directory containing the images
image_dir = 'D:\\1_UQ\\Images\\2022_12 Heron\\2022_12_11_russID_site\\all_ID'

# Get the list of image files and sort them naturally
img_paths = natsorted([os.path.join(image_dir, f) for f in os.listdir(image_dir) if f.endswith(('png', 'JPG', 'jpeg'))])

# Extract the image labels (filenames without the directory path)
img_labels = [os.path.basename(f) for f in img_paths]

# Number of images
num_images = len(img_paths)

# Calculate the grid size (assuming we want to display them in a square grid)
grid_size = int(num_images**0.5) + 1

# Create a figure and axis objects
fig, axes = plt.subplots(7, 4, figsize=(15, 15))

# Flatten axes for easy iteration if they are not already flat
axes = axes.flatten()

# Iterate over the images and display them
for idx, img_path in enumerate(img_paths):
    # Load the image
    img = mpimg.imread(img_path)
    
    # Display the image on the grid
    axes[idx].imshow(img)
    axes[idx].axis('off')  # Hide the axis
    # Add the image label
    axes[idx].annotate(img_labels[idx], (10, 10), color='white', fontsize=12, ha='left', va='top',
                       bbox=dict(facecolor='black', alpha=0.5))

# Hide any remaining axes if the number of images is not a perfect square
for idx in range(num_images, len(axes)):
    axes[idx].axis('off')

# Adjust layout
plt.tight_layout()
plt.show()



############
# Define the path to your Python script
python_script <- "path/to/display_images.py"

# Use reticulate to source and run the Python script
source_python(python_script)

