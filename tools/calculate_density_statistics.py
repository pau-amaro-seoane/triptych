import numpy as np
import sys

# Script: calculate_densities.py
# Author: Pau Amaro Seoane
# Date: 21 August 2024
#
# Description:
# ------------
# This Python script calculates the average, maximum, and minimum densities from a specified results file.
# The results file is expected to contain simulation data in the following format:
#
# Model 1 | Model 2 | Velocity (km/s) | Periastron | Initial Separation | Max Density (g/cm^3)
#
# Each row corresponds to a specific simulation, and the script extracts the "Max Density" values 
# (the last column) to compute the statistics.
#
# The script reads the file, processes the density values, and then calculates:
# 1. The average density
# 2. The maximum density
# 3. The minimum density
#
# The script outputs these statistics in scientific notation.
#
# Usage:
# ------
# To run the script, you can use the following command:
#
#   $ python calculate_densities.py results.txt
#
# Where `results.txt` is the name of your input file containing the simulation results.
#
# If no file is provided, the script will print an error message and exit.

def read_densities(filename):
    """
    Reads the densities from the specified results file.

    Parameters:
    filename (str): The name of the results file to read.

    Returns:
    list: A list of density values extracted from the file.
    """
    densities = []

    with open(filename, 'r') as file:
        for line in file:
            # Skip comments, empty lines, and lines with "Simulation failed"
            if line.startswith('#') or line.strip() == '' or "Simulation failed" in line:
                continue

            # Extract the density (last column in the file)
            try:
                density = float(line.split()[-1])
                densities.append(density)
            except ValueError:
                print(f"Warning: Skipping invalid density value in line: {line.strip()}")
    
    return densities

def calculate_density_stats(densities):
    """
    Calculates and prints the average, maximum, and minimum densities.

    Parameters:
    densities (list): A list of density values.
    """
    if not densities:
        print("Error: No valid density values found.")
        return

    # Convert the list to a NumPy array for easy calculation
    densities_array = np.array(densities)

    # Calculate the statistics
    avg_density = np.mean(densities_array)
    max_density = np.max(densities_array)
    min_density = np.min(densities_array)

    # Print the results
    print(f"Average density: {avg_density:.8e} g/cm^3")
    print(f"Maximum density: {max_density:.8e} g/cm^3")
    print(f"Minimum density: {min_density:.8e} g/cm^3")

def main():
    """
    Main function that runs the script to calculate and display density statistics.
    """
    # Check if the filename is provided as a command-line argument
    if len(sys.argv) != 2:
        print("Usage: python calculate_densities.py <results_file>")
        sys.exit(1)

    # Get the filename from the command-line argument
    filename = sys.argv[1]

    # Read the densities from the file
    densities = read_densities(filename)

    # Calculate and print the statistics
    calculate_density_stats(densities)

# Run the script
if __name__ == "__main__":
    main()

