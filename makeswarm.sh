#!/bin/bash


# Check if the genelist.txt file exists
if [ ! -f path/to/genelist.txt ]; then
    echo "file not found!"
    exit 1
fi

swarm_file="fit_scRNA.swarm"

while IFS= read -r gene; do
    # Submit the job for each gene
   echo "julia -t 2 -p 11 fitscript_scRNA.jl      $gene > swarm_outputs/$gene.o" >> "$swarm_file"

done < path/to/genelist.txt

echo "Swarm file created: $swarm_file"
