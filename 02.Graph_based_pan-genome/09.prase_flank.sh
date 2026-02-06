#!/bin/bash

input_file="panda_annota.gff.gene"
output_file="3.5kbflk.bed"

awk -F'\t' '
$3 == "gene" {
    chrom = $1;
    start = $4 - 1;  # Convert to 0-base
    end = $5;
    strand = $7;
    
    # Calculate upstream and downstream areas
    if (strand == "+") {
        upstream_start = (start - 5000 > 0) ? start - 5000 : 0;
        upstream_end = start;
        downstream_start = end;
        downstream_end = end + 5000;
    } else {
        upstream_start = end;
        upstream_end = end + 5000;
        downstream_start = (start - 5000 > 0) ? start - 5000 : 0;
        downstream_end = start;
    }
    
    print chrom, upstream_start, upstream_end;
    print chrom, downstream_start, downstream_end;
}' "$input_file" > "$output_file"

echo "$output_file"
echo "Format per line: Chromosome Start Point End Point"