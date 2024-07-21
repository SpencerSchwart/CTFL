import os
from moviepy.editor import ImageSequenceClip

# Directory containing your images
image_folder = '/home/spencer/basilisk/CTFL/data/blur/blur3D/mode4/lvl9/img/f'

# List all files in the directory
image_files = [f for f in os.listdir(image_folder) if f.endswith('.png')]

# Sort files based on the numerical value in the filename
sorted_files = sorted(image_files, key=lambda x: float(x.split('-')[1].replace('.png', '')))

# Create full paths to the image files
image_paths = [os.path.join(image_folder, f) for f in sorted_files]

# Define frames per second
fps = 10

# Create the video clip
clip = ImageSequenceClip(image_paths, fps=fps)

# Write the video file
clip.write_videofile('output.mp4', codec='libx264')
