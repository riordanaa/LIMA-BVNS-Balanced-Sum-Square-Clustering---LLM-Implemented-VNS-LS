#!/bin/bash

# Create directories for results
mkdir -p results
mkdir -p assignments

# Find the executable (checking multiple possible locations)
if [ -f "./lima_vns_64" ]; then
    executable="./lima_vns_64"
elif [ -f "./executables/lima_vns_64" ]; then
    executable="./executables/lima_vns_64"
elif [ -f "../executables/lima_vns_64" ]; then
    executable="../executables/lima_vns_64"
elif [ -f "./src/lima_vns_64" ]; then
    executable="./src/lima_vns_64"
else
    echo "Error: Could not find lima_vns_64 executable!"
    echo "Searching for executable in common locations..."
    find . -name "lima_vns_64" -type f 2>/dev/null
    echo "Please specify the correct path to the executable."
    exit 1
fi

echo "Using executable: $executable"

# Path to initial solutions
init_solutions_dir="initial_solutions"

# Check if initial solutions directory exists
if [ ! -d "$init_solutions_dir" ]; then
    echo "Warning: Initial solutions directory not found: $init_solutions_dir"
    echo "Continuing without initial solutions. Solutions will be generated randomly."
    init_solutions_dir=""
fi

# Dataset parameters (format: filename,clusters,time_limit)
datasets=(
  #"iris.csv,3,0.42"
  #"wine.csv,3,0.60"
  #"glass.csv,7,2.54"
  #"thyroid.csv,3,1.45"
  #"ionosphere.csv,2,3.65"
  #"libra.csv,15,8.36"
  #"user_knowledge.csv,4,12.35"
  #"body.csv,2,10.51"
  #"water.csv,13,24.04"
  #"breast_cancer.csv,2,22.12"
  #"synthetic_control.csv,6,16.98"
  #"vehicle.csv,6,45.66"
  #"vowel.csv,11,148.36"
  #"yeast.csv,10,394.09"
  "multiple_features_reduced.csv,10,2233.19"
  "image_segmentation.csv,7,1827.98"
)

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
  
  # Check if initial solutions exist for this dataset
  use_init_solutions=""
  if [ -n "$init_solutions_dir" ]; then
    if [ -f "$init_solutions_dir/${dataset_name}-init1.bin" ]; then
      use_init_solutions="$init_solutions_dir"
      echo "Using pre-generated initial solutions from $use_init_solutions"
    else
      echo "Warning: No initial solutions found for $dataset_name, will generate dynamically"
    fi
  fi
 
  for run in {1..10}; do
    echo "  Run $run/10..."
    seed=$((1000 + run))
    
    # Determine which initial solution file to use for this run
    init_solution_param=""
    if [ -n "$use_init_solutions" ]; then
      init_solution_param="$use_init_solutions"
    fi
    
    # Run the algorithm with the executable path we found, passing initial solution directory if available
    $executable "datasets/$filename" $clusters $time_limit 1 $seed "results/${dataset_name}-run$run" "assignments/${dataset_name}-run$run" $init_solution_param
    
    # If the command fails, print an error but continue to the next run
    if [ $? -ne 0 ]; then
      echo "Error: Command failed for $dataset_name (run $run)! Check if the dataset file exists."
      ls -la "datasets/$filename" 2>/dev/null || echo "File not found!"
    fi
    
    echo "----------------------------------------"
  done
 
  echo "Completed all runs for $dataset_name"
  echo "========================================================"
 
  # Increment dataset counter
  ((current_dataset++))
done

echo "All experiments completed!"
