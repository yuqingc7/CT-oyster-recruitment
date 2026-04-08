#!/bin/bash

input="geno_n1601_in_silico_wild_AQ_selected5_905.recode.strct_in"
wild_aq="centroid_refs.txt"

mkdir -p structure_inputs

# Save the first two lines (header)
head -n 2 "$input" > header.txt

# Extract WILD and AQ individuals
grep -E "^WILD-|^AQ-" "$input" > "$wild_aq"

# Extract remaining individuals
tail -n +3 "$input" | grep -v -E "^WILD-|^AQ-" > query.txt

# Loop through remaining individuals
while read -r individual; do
    # Get individual ID (first field)
    ind_id=$(echo "$individual" | awk '{print $1}')
    
    # Create STRUCTURE input for this run
    output="structure_run_${ind_id}.strct_in"

    # Write header
    cat header.txt > "$output"
    # Append WILD, AQ, and the current individual
    cat "$wild_aq" >> "$output"
    echo "$individual" >> "$output"

done < query.txt

mv structure_run* structure_inputs/
