"""
Program: plot_heatmap.py
Pau Amaro Seoane
Berlin 17 August 2024

Description:
------------
This Python script generates a heatmap that visualizes the core density of stellar models 
as a function of various simulation parameters, such as relative velocity and periastron separation. 
The script reads the simulation results from a file named 'results.txt', which is generated 
by running the accompanying simulation script (e.g., run_simulations.sh). 

The heatmap is rendered using LaTeX for text formatting, ensuring high-quality labels and 
annotations. The script can be customized to adjust font sizes and other visual aspects 
of the plot.

The heatmap is displayed with core density on the color bar and the axes labeled with 
velocity and periastron separation.

Requirements:
-------------
1. Python 3.x
2. Required Python libraries:
   - pandas: For handling and processing the data.
   - seaborn: For generating the heatmap.
   - matplotlib: For plotting and LaTeX rendering.

   These can be installed via pip:

   $ pip install pandas seaborn matplotlib


3. A 'results.txt' file in the current directory. This file should be structured as follows:
- Model1: Stellar model for the first star.
- Model2: Stellar model for the second star.
- Velocity: Relative velocity at infinity (in km/s).
- Periastron: Periastron separation (in units of R_1 + R_2).
- Initial Separation: Initial separation (normalized to the sum of the parent star radii).
- Core Radius: Radius at which the density starts to decrease significantly, indicating the core boundary (in R_sun).
- Core Density: Density at the core radius (in g/cm^3).

Each row of this file corresponds to a simulation, with values separated by spaces.

Usage:
------
1. Ensure that the 'results.txt' file is present in the current directory.
2. Run the script:

    $ python plot_heatmap.py


3. The heatmap will be displayed, showing core density as a function of the chosen parameters. 
The plot can be customized by editing the script, such as changing font sizes or plot dimensions.

Customization:
--------------
- Font Size: You can adjust the font size for labels and annotations by modifying the `fontsize` 
and `annot_kws` parameters in the script.
- LaTeX Rendering: Ensure that LaTeX is installed on your system for proper rendering of text 
within the plot. If LaTeX is not installed, remove or comment out the LaTeX-related lines 
in the script.
- Axes and Labels: The script can be further customized to plot different parameters by adjusting 
the pivot table creation or the heatmap plotting sections.
"""

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

# Load the results file into a DataFrame
df = pd.read_csv('./results.txt', delim_whitespace=True, comment='#',
                 names=['Model1', 'Model2', 'Velocity', 'Periastron', 'Initial_Separation', 'Core_Radius', 'Core_Density'])

# Convert Core_Density to float (remove scientific notation)
df['Core_Density'] = df['Core_Density'].astype(float)

# Pivot table to create a 2D grid for heatmap
# We will plot Core_Density as a function of Velocity and Periastron
heatmap_data = df.pivot_table(index='Periastron', columns='Velocity', values='Core_Density', aggfunc=np.mean)

# Set up LaTeX for rendering
plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)  # Change 'size' to adjust the global font size

# Plot the heatmap without annotations in the cells
plt.figure(figsize=(12, 10))
sns.heatmap(heatmap_data, cmap="YlGnBu", annot=False, linewidths=.5, 
            cbar_kws={"label": r'$\mathrm{Core\ Density\ (g/cm^3)}$'})  # Simplified LaTeX syntax

# Customize the axes labels with LaTeX and adjustable font size
plt.title(r'$\mathrm{Core\ density\ as\ a\ function\ of\ velocity\ and\ periapsis distance}$', fontsize=18)
plt.xlabel(r'$\mathrm{Velocity\ (km/s)}$', fontsize=16)
plt.ylabel(r'$\mathrm{Periapsis distance\ (R_1\ +\ R_2)}$', fontsize=16)

# Display the plot
plt.show()

