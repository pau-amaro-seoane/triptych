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
# It also checks that R_c does not exceed 500 R_sun and calculates the thermal timescale
# using the Eddington luminosity. Furthermore, the script calculates the opacity using either:
# (1) A simple Thomson scattering formula (electron scattering), where 
#     kappa = 0.2 * (1 + X) cm^2/g, or 
# (2) An alternative opacity formula provided by Buchler and Yueh (1976), 
#     which is more detailed and takes into account various physical effects.
# 
# The user is prompted to select which opacity calculation to use.
# The hydrogen fraction (X) is extracted from column 5 in the data file, from the row
# corresponding to the specified core radius (R_c).
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
# 3. When prompted, input the core radius (R_c) in R_sun and select the opacity calculation.
#
# Example:
# --------
# $ ./calculate_timescale.sh
# Enter the core radius R_c (in R_sun): 0.05
# What opacity do you want to use?
# (1) kappa as in Thomson (electron scattering), k = 0.2(1+X) cm^2/g
# (2) kappa as estimated by Buchler and Yueh (1976)
# Enter 1 or 2: 1
# Core luminosity L_c: 3027.3 L_sun
# Core Eddington Luminosity L_Edd_c: xxxx L_sun
# Thermal timescale calculated with L_Edd_c: xxxx years
# Real thermal timescale t_KH: 501.015 years

# Constants
G=6.67430e-8    # Gravitational constant in cm^3 g^-1 s^-2
a=7.564e-15     # Radiation constant in erg cm^-3 K^-4
c=2.998e10      # Speed of light in cm/s
M_sun=1.989e33  # Solar mass in grams
R_sun=6.957e10  # Solar radius in cm
L_sun=3.828e33  # Solar luminosity in erg/s
seconds_in_year=3.154e7  # Number of seconds in a year

# Prompt for core radius input in R_sun
read -p "Enter the core radius R_c (in R_sun): " R_c_input

# Check if R_c exceeds 500 R_sun
if [ "$(awk 'BEGIN{print ('$R_c_input' > 500)}')" -eq 1 ]; then
  echo "Error: The core radius R_c exceeds 500 R_sun. Please enter a smaller value."
  exit 1
fi

# Extract the row corresponding to the core radius R_c
line=$(awk -v Rc="$R_c_input" 'BEGIN{min_diff=1e9} {diff=$4-Rc; if(diff<0) diff=-diff; if(diff<min_diff) {min_diff=diff; closest=$0}} END{print closest}' product.dat)

# Extract the next row for R_(c+1)
next_line=$(awk -v Rc="$R_c_input" 'BEGIN{found=0} {if(found==1) {print $0; exit} if($4==Rc) found=1}' product.dat)

# Extract variables from the lines
m_c=$(echo $line | awk '{print $1}')
P_c=$(echo $line | awk '{print $2}')
rho_c=$(echo $line | awk '{print $3}')
R_c=$(echo $line | awk '{print $4}')
X_c=$(echo $line | awk '{print $5}')
P_c1=$(echo $next_line | awk '{print $2}')
R_c1=$(echo $next_line | awk '{print $4}')

# Convert units
M_c=$(awk -v m_c="$m_c" -v M_sun="$M_sun" 'BEGIN {print m_c * M_sun}')
R_c=$(awk -v R_c="$R_c" -v R_sun="$R_sun" 'BEGIN {print R_c * R_sun}')
R_c1=$(awk -v R_c1="$R_c1" -v R_sun="$R_sun" 'BEGIN {print R_c1 * R_sun}')

# Calculate the temperatures
T_c=$(awk -v P_c="$P_c" -v a="$a" 'BEGIN {print (3 * P_c / a) ^ 0.25}')
T_c1=$(awk -v P_c1="$P_c1" -v a="$a" 'BEGIN {print (3 * P_c1 / a) ^ 0.25}')

# Prompt for the opacity calculation method
echo "What opacity do you want to use?"
echo "(1) kappa as in Thomson (electron scattering), k = 0.2(1+X) cm^2/g"
echo "(2) kappa as estimated by Buchler and Yueh (1976)"
read -p "Enter 1 or 2: " opacity_choice

# Calculate the opacity based on the user's choice
if [ "$opacity_choice" -eq 1 ]; then
  kappa=$(awk -v X="$X_c" 'BEGIN {print 0.2 * (1 + X)}')
  echo "To calculate the opacity we have used kappa with X=$X_c."
elif [ "$opacity_choice" -eq 2 ]; then
  kappa=$(awk -v X="$X_c" -v rho="$rho_c" -v T="$T_c" 'BEGIN {
      term1 = 0.2 * (1 + X);
      term2 = 1 + 2.7e11 * rho / (T * T);
      term3 = 1 + (T / 4.5e8) ^ 0.86;
      kappa_BY = term1 * (1 / term2) * (1 / term3);
      print kappa_BY
  }')
  echo "To calculate the opacity we have used kappa_BY with X=$X_c."
else
  echo "Invalid choice. Exiting."
  exit 1
fi

# Calculate the luminosity L_c in erg/s using the selected kappa
L_c=$(awk -v a="$a" -v c="$c" -v kappa="$kappa" -v R_c="$R_c" -v rho_c="$rho_c" -v T_c="$T_c" -v T_c1="$T_c1" -v R_c1="$R_c1" 'BEGIN {
    L = - (16 * 3.141592653589793 * a * c) / (3 * kappa) * (R_c ^ 2) / rho_c * (T_c ^ 3) * (T_c - T_c1) / (R_c - R_c1);
    print L
}')

# Convert the luminosity to solar units (L_sun)
L_c_Lsun=$(awk -v L_c="$L_c" -v L_sun="$L_sun" 'BEGIN {print L_c / L_sun}')

# Calculate the Eddington luminosity L_Edd in L_sun
L_Edd_Lsun=$(awk -v m_c="$m_c" 'BEGIN {print 3.8e4 * m_c}')

# Calculate the thermal (Kelvin-Helmholtz) timescale t_KH using the Eddington luminosity in seconds
t_KH_Edd=$(awk -v G="$G" -v M_c="$M_c" -v R_c="$R_c" -v L_Edd_Lsun="$L_Edd_Lsun" -v L_sun="$L_sun" 'BEGIN {
    t = G * (M_c ^ 2) / (2 * R_c * L_Edd_Lsun * L_sun);
    print t
}')

# Convert the timescale to years for Eddington luminosity
t_KH_Edd_years=$(awk -v t_KH_Edd="$t_KH_Edd" -v seconds_in_year="$seconds_in_year" 'BEGIN {print t_KH_Edd / seconds_in_year}')

# Calculate the real thermal (Kelvin-Helmholtz) timescale t_KH in seconds
t_KH=$(awk -v G="$G" -v M_c="$M_c" -v R_c="$R_c" -v L_c="$L_c" 'BEGIN {
    t = G * (M_c ^ 2) / (2 * R_c * L_c);
    print t
}')

# Convert the timescale to years for the real luminosity
t_KH_years=$(awk -v t_KH="$t_KH" -v seconds_in_year="$seconds_in_year" 'BEGIN {print t_KH / seconds_in_year}')

# Output the results
echo "Core luminosity L_c: $L_c_Lsun L_sun"
echo "Core Eddington Luminosity L_Edd_c: $L_Edd_Lsun L_sun"
echo "Thermal timescale calculated with L_Edd_c: $t_KH_Edd_years years"
echo "Real thermal timescale t_KH: $t_KH_years years"
