#!/bin/bash

###############################################################################
# Script Name: run_pipeline.sh
# Description: This script executes a series of Python scripts as part of a
#              data processing pipeline. The script requires a configuration 
#              file as input and optionally takes a flag to run alternative 
#              steps in the pipeline.
#
# Usage:
#   ./run_pipeline.sh <config_file_path> [-nnk True]
#
# Parameters:
#   - <config_file_path>: (Required) Path to the configuration file to be used.
#   - -nnk True: (Optional) If provided, executes S2_NNK and S3_NNK instead 
#                of S2 and S3.
#
# Example:
#   ./run_pipeline.sh Python_Scripts/configs/example_config.py
#
# Author: Jaro Steindorff
###############################################################################

# Function to display usage information
usage() {
  echo "Usage: $0 <config_file_path> [-nnk True]"
  exit 1
}

# Check if the required parameter (config file path) is provided
if [ $# -lt 1 ]; then
  echo "Error: Missing required argument <config_file_path>"
  usage
fi

# Initialize variables
CONFIG_FILE_PATH=$1
NNK_FLAG=false

# Parse optional parameters
shift # Shift past the first argument (config file path)
while [[ "$#" -gt 0 ]]; do
  case $1 in
    -nnk)
      if [ "$2" == "True" ]; then
        NNK_FLAG=true
        shift 2
      else
        echo "Error: Invalid value for -nnk. Use 'True' to enable."
        usage
      fi
      ;;
    *)
      echo "Error: Unknown parameter $1"
      usage
      ;;
  esac
done

# Check if the configuration file exists
if [ ! -f "$CONFIG_FILE_PATH" ]; then
  echo "Error: Configuration file '$CONFIG_FILE_PATH' not found."
  exit 1
fi

# Update the main configuration file
echo "Setting up configuration from $CONFIG_FILE_PATH..."
cp "$CONFIG_FILE_PATH" Python_Scripts/config.py

# Execute the pipeline steps
if [ "$NNK_FLAG" = false ]; then
echo "Running S1.py..."
./Python_Scripts/S1.py

if [ "$NNK_FLAG" = true ]; then
  echo "Running S2_NNK.py..."
  ./Python_Scripts/S2_NNK.py
  echo "Running S3_NNK.py..."
  ./Python_Scripts/S3_NNK.py
else
  echo "Running S2.py..."
  ./Python_Scripts/S2.py
  echo "Running S3.py..."
  ./Python_Scripts/S3.py
fi

echo "Running S4.py..."
./Python_Scripts/S4.py

echo "Running S5.py..."
./Python_Scripts/S5.py

echo "Pipeline execution completed."
