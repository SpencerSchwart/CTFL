# Set the terminal and output
set terminal pngcairo size 800,600 enhanced font "Arial,10"
set output 'main_with_inset_fixed.png'

# Set the main plot properties
set xrange [0:1]
set yrange [0:1.2]
set xlabel "X-axis"
set ylabel "Y-axis"
set style line 1 lc rgb "blue" lt 1 lw 2

# Begin multiplot mode
set multiplot

# Main plot
plot x**2 with lines ls 1 title 'Main plot'

# Define the inset
set origin 0.6,0.6  # Position of the inset (x, y)
set size 0.3,0.3    # Size of the inset (width, height)
set xrange [0.1:0.3]   # Adjusted range for x-axis in the inset
set yrange [0.01:0.09] # Adjusted range for y-axis in the inset
unset xlabel
unset ylabel
set border lw 1  # Adjust border line width for inset

# Plot inside the inset
plot x**2 with lines ls 1 notitle

# End multiplot mode
unset multiplot

# Close the output
set output
