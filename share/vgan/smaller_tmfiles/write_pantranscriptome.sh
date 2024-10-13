#!/bin/bash

# Check if input file is provided
if [ -z "$1" ]; then
    echo "Usage: $0 input_tsv_file"
    exit 1
fi

printf "Name\tLength\tTranscript\tHaplotypes\n" > pantranscriptome.txt
while read -r path_name path_length
do
    # Skip the header line
    if [ "$path_name" == "path_name" ]; then
        continue
    fi

    # Replace Node<number> with #<number>#
    modified_path_name=$(echo "$path_name" | sed 's/Node\([0-9]\+\)/#\1#/g')

    # Add _gbwt_ref_ prefix and _0_0 suffix to the original path name
    final_path_name="_gbwt_ref_${path_name}_0_0"

    # Generate a random string to be used as the transcript identifier
    transcript_id="ENST$RANDOM"

    # Print the final path name, path length from the input file, transcript identifier, and final path name again
    echo -e "$final_path_name\t$path_length\t$transcript_id\t$final_path_name"
done < "$1" >> pantranscriptome.txt

