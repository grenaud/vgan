#!/bin/bash

# Check for input
if [ "$#" -ne 0 ] && [ "$#" -ne 1 ]; then
    echo "Usage: $0 [fastq_file]"
    echo "If no file is provided, input is taken from stdin"
    exit 1
fi

# Function to shuffle FASTQ
shuffle_fastq() {
    paste - - - - | \
    shuf | \
    sed 's/\t/\n/g'
}

# Use input file if provided, otherwise read from stdin
if [ "$#" -eq 1 ]; then
    input_file=$1
    shuffle_fastq < "$input_file"
else
    shuffle_fastq
fi

