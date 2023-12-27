#!/bin/bash

# List of directories
directories=("40_2_1000" "80_2_500" "80_2_1000" "80_2_1000Mean" "80_4_500" "80_4_1000")

# Function to run scripts in the background
run_scripts() {
    cd "$1" || exit 1
    python 2D_directivity.py &
    python 3D_directivity.py &
    python Spectra_SPL.py &
    cd ..
}

# Loop through each directory and run the scripts in the background
for dir in "${directories[@]}"; do
    echo "Running scripts in $dir..."
    run_scripts "$dir"
done

# Wait for all background processes to finish
wait

echo "Script execution completed in all directories."
