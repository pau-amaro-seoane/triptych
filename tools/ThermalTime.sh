#!/bin/bash
#
# Script: calculate_timescale.sh
# Author: Pau Amaro Seoane
# Date: 17 August 2024
#
# Description:
# ------------
# This script calculates the core luminosity and the thermal (Kelvin-Helmholtz) timescale
# for a given stellar core using data from a file named `product.dat`. The script takes 
# as input the radius of the stellar core (R_c) in units of solar radii (R_sun) and outputs 
# the core luminosity in solar luminosities (L_sun) and the thermal timescale in years.
#
# The data file (`product.dat`) should be structured as follows:
#   - Column 1: Mass (m) in units of solar masses (M_sun)
#   - Column 2: Pressure (P) in units of dyne/cm^2
#   - Column 3: Density (rho) in units of g/cm^3
#   - Column 4: Radius (r) in units of solar radii (R_sun)
#   - Column 5: Hydrogen fraction (X)
#   - Column 6: Metallicity (Z)
#   - Column 7: Specific angular momentum (j) in cm^2/s
#
# The script performs the following steps:
# 1. Prompts the user to enter the core radius R_c in units of R_sun.
# 2. Uses `awk` to find the row in `product.dat` that corresponds to the core radius R_c.
#    This row provides the core properties: mass (m_c), pressure (P_c), density (rho_c),
#    and radius (R_c). It also retrieves the next row, corresponding to the properties at
#    radius R_c+1, for use in the calculation.
# 3. Converts the mass and radius values from solar units to CGS units:
#    - Mass (M_c) is converted from M_sun to grams.
#    - Radius (R_c) is converted from R_sun to centimeters.
# 4. Calculates the core temperature (T_c) and the temperature at the next layer (T_c+1)
#    using the provided pressure-temperature relationship:
#    T_c = (3 * P_c / a)^(1/4), where 'a' is the radiation constant.
# 5. Computes the core luminosity (L_c) using the formula:
#    L_c = - (16πac)/(3κ) * (R_c^2)/(rho_c) * T_c^3 * (T_c - T_c+1)/(R_c - R_c+1)
#    The result is in erg/s, and it is then converted to solar luminosities (L_sun).
# 6. Calculates the thermal (Kelvin-Helmholtz) timescale (t_KH) using the formula:
#    t_KH = G(M_c)^2/(2R_cL_c), where 'G' is the gravitational constant.
#    The timescale is initially calculated in seconds and then converted to years.
# 7. Outputs the core luminosity in solar luminosities (L_sun) and the thermal timescale
#    in years.
#
# Requirements:
# -------------
# - The `product.dat` file must be present in the same directory as the script.
# - The script requires a Unix-like environment with `awk` available.
#
# Usage:
# ------
# 1. Make the script executable: chmod +x calculate_timescale.sh
# 2. Run the script: ./calculate_timescale.sh
# 3. When prompted, input the core radius (R_c) in R_sun.
#
# Example:
# --------
# $ ./calculate_timescale.sh
# Enter the core radius R_c (in R_sun): 0.05
# Core Luminosity L_c: 2.35e+02 L_sun
# Thermal (Kelvin-Helmholtz) Timescale t_KH: 3.15e+05 years

# Constants
G=6.67430e-8    # Gravitational constant in cm^3 g^-1 s^-2
a=7.564e-15     # Radiation constant in erg cm^-3 K^-4
c=2.998e10      # Speed of light in cm/s
kappa=0.268     # Opacity in cm^2/g
M_sun=1.989e33  # Solar mass in grams
R_sun=6.957e10  # Solar radius in cm
L_sun=3.828e33  # Solar luminosity in erg/s
seconds_in_year=3.154e7  # Number of seconds in a year

# Prompt for core radius input in R_sun
read -p "Enter the core radius R_c (in R_sun): " R_c_input

# Extract the row corresponding to the core radius R_c
line=$(awk -v Rc="$R_c_input" 'BEGIN{min_diff=1e9} {diff=$4-Rc; if(diff<0) diff=-diff; if(diff<min_diff) {min_diff=diff; closest=$0}} END{print closest}' product.dat)

# Extract the next row for R_(c+1)
next_line=$(awk -v Rc="$R_c_input" 'BEGIN{found=0} {if(found==1) {print $0; exit} if($4==Rc) found=1}' product.dat)

# Extract variables from the lines
m_c=$(echo $line | awk '{print $1}')
P_c=$(echo $line | awk '{print $2}')
rho_c=$(echo $line | awk '{print $3}')
R_c=$(echo $line | awk '{print $4}')
m_c_next=$(echo $next_line | awk '{print $1}')
P_c1=$(echo $next_line | awk '{print $2}')
R_c1=$(echo $next_line | awk '{print $4}')

# Convert units
M_c=$(awk "BEGIN {print $m_c * $M_sun}")
R_c=$(awk "BEGIN {print $R_c * $R_sun}")
R_c1=$(awk "BEGIN {print $R_c1 * $R_sun}")

# Calculate the temperatures
T_c=$(awk "BEGIN {print ((3 * $P_c / $a) ^ 0.25)}")
T_c1=$(awk "BEGIN {print ((3 * $P_c1 / $a) ^ 0.25)}")

# Calculate the luminosity L_c in erg/s
L_c=$(awk "BEGIN {print - (16 * 3.141592653589793 * $a * $c) / (3 * $kappa) * ($R_c ^ 2) / $rho_c * $T_c ^ 3 * ($T_c - $T_c1) / ($R_c - $R_c1)}")

# Convert the luminosity to solar units (L_sun)
L_c_Lsun=$(awk "BEGIN {print $L_c / $L_sun}")

# Calculate the thermal (Kelvin-Helmholtz) timescale t_KH in seconds
t_KH=$(awk "BEGIN {print $G * $M_c ^ 2 / (2 * $R_c * $L_c)}")

# Convert the timescale to years
t_KH_years=$(awk "BEGIN {print $t_KH / $seconds_in_year}")

# Output the results
echo "Core Luminosity L_c: $L_c_Lsun L_sun"
echo "Thermal (Kelvin-Helmholtz) Timescale t_KH: $t_KH_years years"

