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
# The user specifies the ranges for velocities, periastron separations, initial separations, and the 
# stellar models to use. If the user does not input any values, the script takes default values.
#
# The simulation is launched using the command:
#     timeout 3 ../bin/triptych < ./input > ./LOG.dat
#
# If the simulation runs for more than 3 seconds, it is terminated, and the script moves on to the next one.
#
# The results are saved in the `./results.txt` file with the following format:
#     Model1 | Model2 | Velocity (km/s) | Periastron | Initial Separation | Max Density (g/cm^3)
#
# Usage:
# ------
# 1. Make the script executable: `chmod +x run_simulations.sh`
# 2. Run the script and specify the ranges for the parameters, or use the default values:
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
available_models=($(ls *.mdl))

# List available models
echo "Available models:"
for model in "${available_models[@]}"; do
    echo "  - $(basename $model)"
done

# Function to prompt for user input with default values
prompt_with_default() {
    local prompt="$1"
    local default="$2"
    read -p "$prompt [$default]: " input
    # If user provides input, use it; otherwise, use the default value
    echo "${input:-$default}"
}

# Prompt the user to select models to use
selected_models=$(prompt_with_default "Enter the models you want to use (comma-separated list or 'all')" "all")

# If the user selects "all", use all available models
if [ "$selected_models" == "all" ]; then
    models=("${available_models[@]}")
else
    # Otherwise, use only the specified models
    IFS=',' read -r -a models <<< "$selected_models"
fi

# Get ranges from the user or use default values
vel_min=$(prompt_with_default "Enter the minimum velocity (km/s)" "1")
vel_max=$(prompt_with_default "Enter the maximum velocity (km/s)" "7")
vel_step=$(prompt_with_default "Enter the velocity step (km/s)" "0.5")
peri_min=$(prompt_with_default "Enter the minimum periastron separation" "0.1")
peri_max=$(prompt_with_default "Enter the maximum periastron separation" "0.7")
peri_step=$(prompt_with_default "Enter the periastron step" "0.2")
sep_min=$(prompt_with_default "Enter the minimum initial separation" "2")
sep_max=$(prompt_with_default "Enter the maximum initial separation" "5")
sep_step=$(prompt_with_default "Enter the initial separation step" "2")

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
# 6. Max Density: Maximum density in the simulation (in g/cm^3).
#
# Format:
# Model 1 | Model 2 | Velocity (km/s) | Periastron | Initial Separation | Max Density (g/cm^3)
#
EOM

# Function to extract the maximum density
extract_max_density() {
    local product_file=$1

    # Extract the maximum density from the second line, 3rd column of product file (skipping comment lines)
    max_density=$(awk 'NR==2 {print $3}' $product_file)

    echo "$max_density"
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

                    # Run the simulation with a timeout of 3 seconds
                    timeout 3 ../bin/triptych < $input_file > ./LOG.dat

                    # Check if the simulation produced the product.dat file
                    if [ -f "product.dat" ]; then
                        max_density=$(extract_max_density "product.dat")
                        echo "$model1 $model2 $velocity $periastron_sep $initial_sep $max_density" >> "$output_file"
                    else
                        echo "$model1 $model2 $velocity $periastron_sep $initial_sep Simulation failed (timeout or other issue)." >> "$output_file"
                    fi

                    # Clean up temporary files
                    rm -f input dynamics.dat product.dat evolution.dat LOG.dat

                done
            done
        done
    done
done

echo "Simulations completed. Results are saved in $output_file."

