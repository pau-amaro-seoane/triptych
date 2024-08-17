#!/bin/bash
#
# Description:
# ------------
# This script automates the process of running stellar collision simulations using various combinations
# of stellar models, relative velocities, periastron separations, and initial separations.
#
# The stellar models are located in the `../stellar_models/` directory. The script creates symbolic links
# to these models in the current working directory using the `ln -sf` command, so that the simulations
# can access them directly. After the simulations are complete, the temporary files are removed, but the 
# symbolic links to the stellar models remain in place.
#
# The user specifies the ranges for velocities, periastron separations, and initial separations. The script 
# then calculates the total number of simulations based on these ranges and asks for user confirmation 
# before proceeding.
#
# The simulation is launched using the command:
#     timeout 10 ../bin/triptych < ./input > ./LOG.dat
#
# If the simulation runs for more than 4 seconds, it is terminated, and the script moves on to the next one.
#
# The results are saved in the `./results.txt` file with the following format:
#     Model1 | Model2 | Velocity (km/s) | Periastron | Initial Separation | Core Radius (R_sun) | Core Density (g/cm^3)
#
# Usage:
# ------
# 1. Make the script executable: `chmod +x run_simulations.sh`
# 2. Run the script and specify the ranges for the parameters:
#    `./run_simulations.sh`
#
# Example:
# --------
# $ ./run_simulations.sh
#
# The script will calculate the number of simulations and ask for confirmation before proceeding.
#

# Directory for stellar models
model_dir="../stellar_models/"

# Check if the directory exists
if [ ! -d "$model_dir" ]; then
    echo "Error: Stellar models directory '$model_dir' not found."
    exit 1
fi

# Create symbolic links to stellar models in the current directory
for model in $model_dir/*.mdl; do
    ln -sf "$model" .
done

# Find all .mdl files in the current directory (linked models)
models=($(ls *.mdl))

# Get ranges from the user
read -p "Enter the minimum velocity (km/s): " vel_min
read -p "Enter the maximum velocity (km/s): " vel_max
read -p "Enter the velocity step (km/s): " vel_step
read -p "Enter the minimum periastron separation (e.g., 0.1): " peri_min
read -p "Enter the maximum periastron separation (e.g., 0.7): " peri_max
read -p "Enter the periastron step (e.g., 0.1): " peri_step
read -p "Enter the minimum initial separation (e.g., 2): " sep_min
read -p "Enter the maximum initial separation (e.g., 18): " sep_max
read -p "Enter the initial separation step (e.g., 4): " sep_step

# Calculate the number of values for each parameter
num_models=${#models[@]}
num_velocities=$(seq $vel_min $vel_step $vel_max | wc -l)
num_periastron_seps=$(seq $peri_min $peri_step $peri_max | wc -l)
num_initial_seps=$(seq $sep_min $sep_step $sep_max | wc -l)

# Calculate total number of simulations
total_simulations=$((num_models * num_models * num_velocities * num_periastron_seps * num_initial_seps))

# Display the total number of simulations and ask for confirmation
echo "Total number of simulations to be run: $total_simulations"
read -p "Do you want to proceed? (yes/no): " proceed

if [ "$proceed" != "yes" ]; then
    echo "Aborting simulations."
    exit 0
fi

# Generate parameter lists
velocities=($(seq $vel_min $vel_step $vel_max))
periastron_seps=($(seq $peri_min $peri_step $peri_max))
initial_seps=($(seq $sep_min $sep_step $sep_max))

# Output file for results
output_file="./results.txt"

# Write the header to the results file
cat > $output_file <<- EOM
# Simulation Results Summary
# 
# Each row corresponds to a specific simulation with the following parameters:
# 
# 1. Model 1: Stellar model for the first star.
# 2. Model 2: Stellar model for the second star.
# 3. Velocity: Relative velocity at infinity (in km/s).
# 4. Periastron: Periastron separation (in units of R_1 + R_2).
# 5. Initial Separation: Initial separation (normalized to the sum of the parent star radii).
# 6. Core Radius: Radius at which the density starts to decrease significantly, indicating the core boundary (in R_sun).
# 7. Core Density: Density at the core radius (in g/cm^3).
#
# Format:
# Model 1 | Model 2 | Velocity (km/s) | Periastron | Initial Separation | Core Radius (R_sun) | Core Density (g/cm^3)
#
EOM

# Function to extract the core radius and density
extract_core_properties() {
    local product_file=$1

    # Extract the core radius and density
    core_radius=$(awk 'BEGIN {prev_rho = 0} $3 < prev_rho {print $4; exit} {prev_rho = $3}' $product_file)
    core_density=$(awk -v rc=$core_radius '$4 == rc {print $3}' $product_file)

    echo "$core_radius $core_density"
}

# Loop over all combinations of models, velocities, periastron separations, and initial separations
for model1 in "${models[@]}"; do
    for model2 in "${models[@]}"; do
        for velocity in "${velocities[@]}"; do
            for periastron_sep in "${periastron_seps[@]}"; do
                for initial_sep in "${initial_seps[@]}"; do

                    # Create the input file for the simulation
                    input_file="./input"
                    cat > $input_file <<- EOM
$model1  ! input file for parent star 1
$model2  ! input file for parent star 2
$velocity          ! relative velocity at infinity, in km/s
$periastron_sep    ! periastron separation (in units of R_1 + R_2)
$initial_sep       ! initial separation, normalized to the sum of the parent star radii
'dynamics.dat'     ! output data file name for dynamics
'product.dat'      ! output data file name for hydrodynamics
'evolution.dat'    ! output data file name for evolution
EOM

                    # Run the simulation with a timeout of 4 seconds
                    timeout 4 ../bin/triptych < $input_file > ./LOG.dat

                    # Check if the simulation produced the product.dat file
                    if [ -f "product.dat" ]; then
                        core_properties=$(extract_core_properties "product.dat")
                        echo "$model1 $model2 $velocity $periastron_sep $initial_sep $core_properties" >> "$output_file"
                    else
                        echo "Simulation for $model1 $model2 with velocity $velocity km/s and periastron $periastron_sep failed (timeout or other issue)." >> "$output_file"
                    fi

                    # Clean up temporary files
                    rm -f input dynamics.dat product.dat evolution.dat LOG.dat

                done
            done
        done
    done
done

echo "Simulations completed. Results are saved in $output_file."

