#!/bin/bash

> output/all_query_Q_values.tsv
> output/all_query_Q_values_with_pop.tsv

for qfile in output/rep_*/structure_run_*_q; do
    # Get last line, trim leading/trailing whitespace
    last_line=$(tail -n 1 "$qfile" | xargs)
    
    # Extract replicate from folder name
    rep_id=$(basename $(dirname "$qfile"))   # rep_0, rep_1, etc.

    # Get individual ID from filename
    ind_id=$(basename "$qfile" | sed -E 's/structure_run_(.*)_q/\1/')

    # # Extract first Q value (corresponding to first cluster=AQ)
    # first_q=$(echo "$last_line" | awk '{print $3}')

    # Extract first Q value: column 3 if bigger than column 4, else column 4
    first_row=$(head -n 1 "$qfile" | xargs)
    # Conditional: if column 3 in first row > column 4, take column 3 from last line
    first_q=$(awk -v first="$first_row" -v last="$last_line" '
        BEGIN {
            split(first, f)
            split(last, l)
            if (f[3] > f[4]) {
                print l[3]
            } else {
                print l[4]
            }
        }')
    # Print or save to summary file
    echo -e "$ind_id\t$rep_id\t$first_q" >> output/all_query_Q_values.tsv

done

awk '{
    id=$1
    sub(/-[^-]+$/, "", id)       # remove the last "-something" part
    OFS="\t"
    print id, $1, $2, $3
}' output/all_query_Q_values.tsv > output/all_query_Q_values_with_pop.tsv
