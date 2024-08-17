"""
Program: plot_density.py
Author: Pau Amaro Seoane
Date: 17 August 2024

Description:
------------
This Python script generates a heatmap that visualizes the maximum core density of stellar models 
as a function of various simulation parameters, such as relative velocity and periastron separation. 
The script reads the simulation results from a file named 'results.txt', which is generated by running 
the accompanying simulation script (e.g., run_simulations.sh). 

The heatmap is rendered using LaTeX for text formatting, ensuring high-quality labels and 
annotations. The script can be customized to adjust font sizes and other visual aspects 
of the plot.

The heatmap displays core density on the color bar with axes labeled by velocity and periastron separation.

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
- Max Density: Maximum density in the simulation (in g/cm^3).

Each row of this file corresponds to a simulation, with values separated by spaces.

Usage:
------
1. Ensure that the 'results.txt' file is present in the current directory.
2. Run the script:

    $ python plot_density.py

3. The heatmap will be displayed, showing maximum core density as a function of the chosen parameters. 
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
                 names=['Model1', 'Model2', 'Velocity', 'Periastron', 'Initial_Separation', 'Max_Density'])

# Convert Max_Density to float
df['Max_Density'] = df['Max_Density'].astype(float)

# Pivot table to create a 2D grid for heatmap
# We will plot Max_Density as a function of Velocity and Periastron
heatmap_data = df.pivot_table(index='Periastron', columns='Velocity', values='Max_Density', aggfunc="mean")

# Check if the pivot table is empty
if heatmap_data.empty:
    print("Error: No data available for plotting. The pivot table is empty.")
else:
    # Set up LaTeX for rendering
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=20)  # Change 'size' to adjust the global font size

    # Plot the heatmap with a homogeneous palette
    plt.figure(figsize=(12, 10))
    sns.heatmap(heatmap_data, cmap="Reds", annot=False, linewidths=0.0,
                cbar_kws={"label": r'$\mathrm{Max\ density\ (g/cm^3)}$'})  # Simplified LaTeX syntax

    # Customize the axes labels with LaTeX and adjustable font size
    plt.title(r'$\mathrm{Max\ density\ as\ a\ function\ of\ velocity\ and\ periapsis\ distance}$', fontsize=24, pad=13)
    plt.xlabel(r'$\mathrm{Velocity\ (km/s)}$', fontsize=24)
    plt.ylabel(r'$\mathrm{Periapsis\ distance\ (R_1\ +\ R_2)}$', fontsize=24)

    # Display the plot
    plt.show()
