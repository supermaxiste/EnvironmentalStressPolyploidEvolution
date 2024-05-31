#!/bin/bash

# This script takes a bed file and peptide fasta and gives them the same name. 
## It also removes splicing variants and keeps only the primary gene model (require for GENESPACE). 

# Check if correct number of arguments are provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 key peptide outpep"
    exit 1
fi

# Assign filenames to variables
key="$1"
peptide="$2"
renamed_peptide="$3"

# Check if the output files already exist
if [ -e "$renamed_peptide" ]; then
        echo "ERROR: file $renamed_peptide exists."
        exit 1 
fi 


# Some variables to keep track of parsing progression (script is very inefficient):
total=$(cat $key | wc -l) 
line_c=0 # current line counter
percentage=0 # percentage of genes processed


# Loop through each line in the bed file
while IFS= read -r line; do
    # Record and print the progress
    line_c=$((line_c + 1))
    frac=$(echo "scale=2; $line_c / $total" | bc)
    if [ $frac != $percentage ]; then
        percy=$(echo "$frac*100" | bc)
        echo $percy"%"
        percentage=$frac
    fi
    # Extract the transcript name & add a ">" at the start to match name in peptide file 
    transcript=$(echo "$line" | awk '{print $1}') 
    gene=$(echo "$line" | awk '{print $2}') 
    
    # awk script that prints the name and peptide sequence of the current bed entry. 
    # Since the bed only has the main transcript name the command ignores splicing variants.
    awk -v name=$transcript -v gene=$gene '$0 == ">"name {p=1; print ">"gene; next} p && /^>/ {p=0; next} p' $peptide >> $renamed_peptide

done < "$key"
