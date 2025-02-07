#!/bin/bash

# Check if the genelist.txt file exists
# genelist.txt contains list of genes to be simulated
if [ ! -f path/to/genelist.txt ]; then
    echo "file not found!"
    exit 1
fi

swarm_file="fit_smFISH.swarm"

while IFS= read -r gene; do
    # Submit the job for each gene
   echo "julia -t 2 -p 11 fitscript_smFISH.jl      $gene > output_dir/$gene.o" >> "$swarm_file"

done < path/to/genelist.txt

echo "Swarm file created: $swarm_file"
