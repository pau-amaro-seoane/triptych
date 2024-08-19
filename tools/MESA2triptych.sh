#!/bin/bash

# Script: MESA2triptych.sh
# Pau Amaro Seoane
# 19 August 2024, Berlin
#
# Description:
# ------------
# This script converts an ASCII output file from the MESA stellar evolution code
# (specifically a `profile*.data` file) into a format compatible with the triptych
# hydrodynamics code. The triptych format is designed to describe the internal
# structure of a star, where each row corresponds to a specific mass shell within
# the star, starting from the center and moving outward.
#
# Input:
# ------
# The input file should be a MESA `profile*.data` file, which contains data on
# the internal structure of a star at a given time in its evolution. The file is
# typically organized in columns, with each column representing a different physical
# quantity (e.g., mass, radius, pressure, density, etc.) and each row representing
# a different zone or mass shell within the star. The first six lines of the file
# are usually headers or metadata and should be skipped.
#
# Output:
# -------
# The output file is a triptych-compatible file with the following structure:
#
#   # 1:Enclosed mass m
#   # 2:Radius@m
#   # 3:Pressure@m
#   # 4:Density@m
#   # 5:H fraction chem. abundance by mass
#   # 6: Same metals
#   # 7: Same He3
#   # 8: Same C12
#   # 9: Same C13
#   # 10:Same N14
#   # 11:Same N15
#   # 12:Same O16
#   # 13:Same O17
#   # 14:Same O18
#
# The rows in the triptych file represent mass shells within the star, starting
# from the center and moving outward. The script removes the first row in the
# output file, as it typically contains a mass of zero, which could cause issues
# in subsequent analyses.
#
# Procedure:
# ----------
# 1. **Remove the first six lines**: The script uses `tail -n +7` to skip the first
#    six lines of the MESA input file. These lines typically contain headers or
#    metadata that are not needed for the triptych format.
#
# 2. **Column Mapping and Conversion**: The script extracts specific columns from the MESA file
#    and maps them to the corresponding columns in the triptych format. Specifically:
#
#    triptych Column   | MESA Column   | Description
#    ================= | ============= | ====================================
#    1: Enclosed mass  | Column 2      | Mass enclosed within the shell, converted from solar masses to grams
#    2: Radius@m       | Column 3      | Radius at the shell, converted from logR (10^(logR) in Rsun) to cm
#    3: Pressure@m     | Column 6      | Pressure at the shell (converted from logP)
#    4: Density@m      | Column 5      | Density at the shell (converted from logRho)
#    5: H fraction     | Column 7      | Hydrogen fraction by mass
#    6: Metals fraction| Column 9      | Metals fraction by mass
#    7: He3 fraction   | Column 11     | Helium-3 fraction by mass
#    8: C12 fraction   | Column 13     | Carbon-12 fraction by mass
#    9: C13 fraction   | Column 14     | Carbon-13 fraction by mass
#    10: N14 fraction  | Column 15     | Nitrogen-14 fraction by mass
#    11: N15 fraction  | Column 16     | Nitrogen-15 fraction by mass
#    12: O16 fraction  | Column 17     | Oxygen-16 fraction by mass
#    13: O17 fraction  | Column 18     | Oxygen-17 fraction by mass
#    14: O18 fraction  | Column 19     | Oxygen-18 fraction by mass
#
# 3. **Conversion of Logarithmic Values and Mass**:
#    - The MESA file typically stores some quantities, such as radius, pressure, and density, in logarithmic form
#      (logR, logP, logRho). The script converts these values back to their linear form using the expression 10^(log_value).
#    - The enclosed mass is converted from solar masses to grams using the conversion factor:
#      1 solar mass = 1.98847e33 grams.
#    - The radius is converted from Rsun to cm using the conversion factor:
#      1 solar radius = 6.957e10 cm.
#
# 4. **Order Reversal**: The triptych format expects the data to start at the
#    center of the star and move outward. MESA files, however, might list data 
#    from the surface inward. The script reverses the order of the lines using 
#    the `tac` command to ensure that the triptych file starts at the star's center.
#
# 5. **Remove First Line**: The first line in the processed data, which usually 
#    contains a mass of zero, is removed to avoid potential issues in further 
#    analyses.
#
# 6. **Output**: The script writes the processed data to the specified output file,
#    including a header that describes the columns. The output is formatted in
#    scientific notation with eight decimal places for precision.
#
# Example:
# --------
# To convert a MESA output file named `profile8.data` to a triptych format file
# named `triptych_input.txt`, you would use the following command:
#
#   ./MESA2triptych.sh profile8.data triptych_input.txt
#
# This will create the file `triptych_input.txt` with the necessary data in the
# correct format for triptych.
#
# Notes:
# ------
# - Ensure that the input file is properly formatted according to MESA's standards,
#   with the relevant quantities located in the expected columns.
# - This script assumes that the MESA file starts from the surface and moves inward.
#   If your file is ordered differently, additional adjustments might be necessary.
#
# - Note that columns 7-14 required by triptych are not in the default MESA
#   profile*.data files.  To get the information you need, you have to add this
#   line to the `inlist` file of your simulation directory:
#  
#    &star_job
# 
#     read_extra_star_job_inlist(1) = .true.
#     extra_star_job_inlist_name(1) = 'inlist_project'
#     profile_columns_file = 'profile_columns.list' ! <--- This one 
# 
# Then create a `profile_columns.list` with the following contents:
# 
# ```! profile_columns.list -- determines the contents of star model profiles
# ! you can use a non-standard version by setting profile_columns_file in your inlist
# 
# ! units are cgs unless otherwise noted.
# 
#    zone       ! numbers start with 1 at the surface
#    mass       ! m/Msun. mass coordinate of outer boundary of cell.
#    logR       ! log10(radius/Rsun) at outer boundary of zone
#    logT       ! log10(temperature) at center of zone 
#    logRho     ! log10(density) at center of zone
#    logP       ! log10(pressure) at center of zone
#    x_mass_fraction_H
#    y_mass_fraction_He
#    z_mass_fraction_metals
#    h1
#    he3
#    he4
#    c12
#    c13
#    n14
#    n15
#    o16
#    o17
#    o18
#    ne20
#    mg24
#    q
# 
# ! everything below this line is deactivated``` 
################################################################################

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: ./MESA2triptych.sh input_file output_file"
    exit 1
fi

input_file=$1
output_file=$2

# Define conversion factors
msun_to_g=1.98847e33
rsun_to_cm=6.957e10

# Add header to the output file
echo "# 1:Enclosed mass m 2:Radius@m 3:Pressure@m 4:Density@m 5:H fraction chem. abundance by mass 6: Same metals 7: Same He3 8: Same C12 9: Same C13 10:Same N14 11:Same N15 12:Same O16 13:Same O17 14:Same O18" > "$output_file"

# Process the MESA file, skip the first 6 lines, reverse the order of lines, and convert the necessary values
tail -n +7 "$input_file" | awk -v msun_to_g="$msun_to_g" -v rsun_to_cm="$rsun_to_cm" '{
    mass = $2 * msun_to_g;  # Convert from solar masses to grams
    radius = 10^($3) * rsun_to_cm;  # Convert from logR in Rsun to cm
    pressure = 10^($6);  # Convert from logP to pressure in cgs units (dyne/cm^2)
    density = 10^($5);  # Convert from logRho to density in cgs units (g/cm^3)
    h_fraction = $7;  # Hydrogen fraction by mass
    metals = $9;  # Metals fraction by mass
    he3 = $11;  # Helium-3 fraction by mass
    c12 = $13;  # Carbon-12 fraction by mass
    c13 = $14;  # Carbon-13 fraction by mass
    n14 = $15;  # Nitrogen-14 fraction by mass
    n15 = $16;  # Nitrogen-15 fraction by mass
    o16 = $17;  # Oxygen-16 fraction by mass
    o17 = $18;  # Oxygen-17 fraction by mass
    o18 = $19;  # Oxygen-18 fraction by mass
    
    # Print the converted values in scientific notation
    printf("%.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e %.8e\n", 
        mass, radius, pressure, density, h_fraction, metals, he3, c12, c13, n14, n15, o16, o17, o18);
}' | tac | tail -n +2 >> "$output_file"  # Reverse the order and remove the first line

echo "Conversion complete. Output saved to $output_file."

