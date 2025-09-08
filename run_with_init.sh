#!/bin/bash

# Create directories for results
mkdir -p results
mkdir -p assignments

# Path to the executable
executable="./src/lima_vns_64"

# Path to initial solutions directory
init_solutions_dir="initial_solutions"

# Variables to track missing initial solutions
total_dynamic_solutions=0
declare -A missing_solutions_by_dataset

# Print the current directory to debug
echo "Current directory: $(pwd)"
echo "Checking if executable exists: $executable"
if [ -f "$executable" ]; then
	echo "Executable found!"
else
	echo "Executable not found! Searching for it..."
	find . -name "lima_vns_64" 2>/dev/null
	echo "If no executable found, please compile the code first."
	exit 1
fi

# Check if initial solutions directory exists
if [ -d "$init_solutions_dir" ]; then
	echo "Using initial solutions from: $init_solutions_dir"
	echo "Files in initial_solutions directory:"
	ls -la "$init_solutions_dir/"
else
	echo "Warning: Initial solutions directory not found. Will generate solutions dynamically."
	init_solutions_dir=""
	total_dynamic_solutions=$((total_dynamic_solutions + 10 * ${#datasets[@]}))  # All solutions will be dynamic
fi

# Dataset parameters (format: filename,clusters,time_limit)
datasets=(
  "iris.csv,3,0.42"
  "wine.csv,3,0.60"
  "glass.csv,7,2.54"
  "thyroid.csv,3,1.45"
  "ionosphere.csv,2,3.65"
  "libra.csv,15,8.36"
  "user_knowledge.csv,4,12.35"
  "body.csv,2,10.51"          	# Corrected from body_measurements.csv
  "water.csv,13,24.04"        	# Corrected from water_treatment.csv
  "breast_cancer.csv,2,22.12"
  "synthetic_control.csv,6,16.98"
  "vehicle.csv,6,45.66"
  "vowel.csv,11,148.36"
  "yeast.csv,10,394.09"
  "multiple_features_reduced.csv,10,2233.19"  # Corrected from multiple_features.csv
  "image_segmentation.csv,7,1827.98"
)

# Initialize counters for each dataset
for dataset_info in "${datasets[@]}"; do
  IFS=',' read -r filename clusters time_limit <<< "$dataset_info"
  dataset_name=$(basename "$filename" .csv)
  missing_solutions_by_dataset["$dataset_name"]=0
done

# Total datasets to run
total_datasets=${#datasets[@]}
current_dataset=1

# Run 10 times for each dataset
for dataset_info in "${datasets[@]}"; do
  IFS=',' read -r filename clusters time_limit <<< "$dataset_info"
 
  # Extract dataset name without extension for result files
  dataset_name=$(basename "$filename" .csv)
 
  echo "========================================================"
  echo "[$current_dataset/$total_datasets] Running experiments for $dataset_name dataset ($clusters clusters, $time_limit seconds)"
  echo "========================================================"
 
  for run in {1..10}; do
	echo "  Run $run/10..."
	seed=$((1000 + run))
    
    # More verbose check for initial solution files
    solutions_file="$init_solutions_dir/${dataset_name}-init${run}.bin"
    if [ -n "$init_solutions_dir" ] && [ -f "$solutions_file" ]; then
        echo "  Using pregenerated initial solution: $solutions_file"
        solutions_param="$init_solutions_dir"
    else
        echo "  Could not find initial solution: $solutions_file"
        echo "  Files in $init_solutions_dir matching this dataset:"
        ls -la "$init_solutions_dir/${dataset_name}*" 2>/dev/null || echo "  No matching files found!"
        echo "  Using dynamically generated initial solution"
        solutions_param=""
        
        # Track missing solutions
        total_dynamic_solutions=$((total_dynamic_solutions + 1))
        missing_solutions_by_dataset["$dataset_name"]=$((missing_solutions_by_dataset["$dataset_name"] + 1))
    fi
    
	$executable "datasets/$filename" $clusters $time_limit 1 $seed "results/${dataset_name}-run$run" "assignments/${dataset_name}-run$run" $solutions_param
    
	# If the command fails, print an error but continue to the next run
	if [ $? -ne 0 ]; then
  	    echo "Error: Command failed for $dataset_name (run $run)! Check if the dataset file exists."
  	    ls -la "datasets/$filename" 2>/dev/null || echo "File not found!"
  	    # No break here - continue to the next run
	fi
    
	echo "----------------------------------------"
  done
 
  echo "Completed all runs for $dataset_name"
  echo "========================================================"
 
  # Increment dataset counter
  ((current_dataset++))
done

# Print summary of missing initial solutions
echo "========================================================"
echo "SUMMARY OF DYNAMICALLY GENERATED SOLUTIONS"
echo "========================================================"
echo "Total dynamically generated solutions: $total_dynamic_solutions"
echo ""
echo "Breakdown by dataset:"
for dataset_name in "${!missing_solutions_by_dataset[@]}"; do
  count=${missing_solutions_by_dataset["$dataset_name"]}
  if [ $count -gt 0 ]; then
    echo "  $dataset_name: $count out of 10 runs"
  fi
done
echo "========================================================"

echo "All experiments completed!"
