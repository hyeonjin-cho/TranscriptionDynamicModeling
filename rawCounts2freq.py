'''
Use this code to make pdf file from raw count matrix

run in command line:
python-scanpy rawCounts2freq.py "path_to_input_csv_file" --output_dir <path_to_output_dir>
'''

import argparse
import os
import pandas as pd
import numpy as np

# Function to compute frequency distribution and include all count levels
def compute_frequency(data, gene_name, max_count_level):
    # Count occurrences of each level
    counts = data.value_counts().sort_index()

    # Create a full range of count levels (from 0 to max_count_level)
    full_range = pd.Series(0, index=np.arange(max_count_level + 1))

    # Update the full range with actual counts
    full_range.update(counts)

    # Calculate frequencies
    frequencies = full_range / full_range.sum()

    # Prepare result DataFrame
    result = pd.DataFrame({
        "Frequency": frequencies.values,
        "Count Level": frequencies.index
    })
    return result

# Set up argparse for handling input arguments
parser = argparse.ArgumentParser(description="Process a gene x cell count matrix.")
parser.add_argument("input_csv", help="Path to the input CSV file containing the gene x cell count matrix.")
parser.add_argument("--output_dir", default=".", help="Directory where output files will be saved (default: current directory).")
args = parser.parse_args()

# Load the gene x cell count matrix
input_file = args.input_csv
output_dir = args.output_dir
data = pd.read_csv(input_file, index_col=0)

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Process each gene
for gene in data.index:
    # Extract counts for the current gene
    counts = data.loc[gene, :].astype(int)  # Ensure counts are integers

    # Find the maximum count level in the dataset
    max_count_level = int(counts.max())

    # Compute frequency distribution
    frequency_data = compute_frequency(counts, gene, max_count_level)

    # Save to a text file
    output_file = os.path.join(output_dir, f"{gene}.txt")
    frequency_data.to_csv(output_file, sep="\t", index=False, header=False)
