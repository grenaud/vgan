#!/bin/bash

# Usage instructions
if [ "$#" -lt 2 ] || [ "$#" -gt 3 ]; then
    echo "Usage: $0 <FASTA_FILE> <prefix> [substring (optional)]"
    exit 1
fi

FASTA_FILE=$1
prefix=$2
substring=${3:-}  # Optional third argument

# Check if the FASTA file exists
if [ ! -f "$FASTA_FILE" ]; then
    echo "Error: FASTA file '$FASTA_FILE' not found."
    exit 1
fi

# Check if necessary Python scripts are available
if [ ! -f "rotate_by_substring.py" ]; then
    echo "Error: Python script 'rotate_by_substring.py' not found."
    exit 1
fi

if [ ! -f "extractAncSeq.py" ]; then
    echo "Error: Python script 'extractAncSeq.py' not found."
    exit 1
fi

# Ensure required tools are installed
required_tools=("mafft" "raxmlHPC-PTHREADS" "iqtree" "vg")
for tool in "${required_tools[@]}"; do
    if ! command -v $tool &> /dev/null; then
        echo "Error: Required tool '$tool' not found in PATH."
        exit 1
    fi
done

# Check for soibean_db.baseFreq
if [ ! -f "soibean_db.baseFreq" ]; then
    echo "Warning: 'soibean_db.baseFreq' file not found. Creating a new file."
    wget https://raw.githubusercontent.com/grenaud/vgan/refs/heads/main/share/vgan/soibean_dir/soibean_db.baseFreq
fi

echo "All checks passed! Starting the bioinformatics analysis now..."

# Clean and normalize sequences
awk '/^>/{print;next}{gsub(/[^ATCGNatcgn-]/,"N"); print}' "$FASTA_FILE" > "${prefix}_clean.fa"
awk '/^>/ { gsub(/ /, "_", $0); gsub(/[^a-zA-Z0-9_.>]/, "", $0); print; next } { print }' "${prefix}_clean.fa" > "${prefix}_clean2.fa"

# Call Python script for rotation (with or without optional substring)
if [ -n "$substring" ]; then
    python rotate_by_substring.py "${prefix}_clean2.fa" "$prefix" "$substring"
else
    python rotate_by_substring.py "${prefix}_clean2.fa" "$prefix"
fi

# Align and convert to uppercase
mafft --auto "${prefix}_rotated.fa" > "${prefix}_aligned.fa"
awk '/^>/ {print; next} {print toupper($0)}' "${prefix}_aligned.fa" > "${prefix}_uppercase.fa"
MSA_FILE="${prefix}_uppercase.fa"

# Calculate base frequencies
A=$(grep -v ">" "$MSA_FILE" | tr -cd 'Aa' | wc -c)
C=$(grep -v ">" "$MSA_FILE" | tr -cd 'Cc' | wc -c)
G=$(grep -v ">" "$MSA_FILE" | tr -cd 'Gg' | wc -c)
T=$(grep -v ">" "$MSA_FILE" | tr -cd 'Tt' | wc -c)

totalLength=$((A + C + G + T))
if [ "$totalLength" -eq 0 ]; then
    echo "Error: Total length of bases in '$MSA_FILE' is zero."
    exit 1
fi

totalA=$(echo "scale=2; $A/$totalLength" | bc)
totalC=$(echo "scale=2; $C/$totalLength" | bc)
totalG=$(echo "scale=2; $G/$totalLength" | bc)
totalT=$(echo "scale=2; $T/$totalLength" | bc)

# Run RAxML and IQ-TREE
nice -19 raxmlHPC-PTHREADS -s "$MSA_FILE" -m GTRCATX -n "$prefix" -p 76 -T 40 -d --HKY85
iqtree -s "$MSA_FILE" -m HKY+G -te "RAxML_bestTree.${prefix}" --prefix "${prefix}.anc" -asr
python extractAncSeq.py "${prefix}.anc.state" "${prefix}.anc.fa"

# Clean ancestral sequences
awk 'BEGIN{FS=""} /^>/{print;next} {for(i=1;i<=NF;i++) {if($i ~ /[ACGTacgt]/) printf toupper($i); else printf "N"}; printf "\n"}' "${prefix}.anc.fa" > "${prefix}.anc.clean.fa"
cat "$MSA_FILE" "${prefix}.anc.clean.fa" | mafft /dev/stdin > "${prefix}.all.fa"
awk '/^>/ {print; next} {print toupper($0)}' "${prefix}.all.fa" > "${prefix}.all2.fa"

# Create tree_dir if needed
mkdir -p tree_dir

# Rename and move tree
sed 's/^>Node\([0-9]\+\)/>N\1'"${prefix}"'/' "${prefix}.all2.fa" > "${prefix}_renamed.anc.fa"
sed 's/Node\([0-9]\+\)/N\1'"${prefix}"'/g' "${prefix}.anc.treefile" > "tree_dir/${prefix}.new.dnd"

# Graph construction and indexing
vg construct -M "${prefix}_renamed.anc.fa" > "${prefix}_pre.vg"
vg mod -X 5 "${prefix}_pre.vg" > "${prefix}.vg"
vg snarls "${prefix}.vg" > "${prefix}.snarls"
vg view "${prefix}.vg" > "${prefix}.gfa"
vg gbwt -o "${prefix}.gbwt" -g "${prefix}.gg" -G "${prefix}.gfa"
vg index -j "${prefix}.dist" "${prefix}.vg"
vg minimizer -g "${prefix}.gbwt" -i "${prefix}.con.min" -k 20 -w 10 "${prefix}.vg"
vg minimizer -g "${prefix}.gbwt" -i "${prefix}.rec.min" -k 17 -w 8 "${prefix}.vg"
vg minimizer -g "${prefix}.gbwt" -i "${prefix}.sen.min" -k 10 -w 3 "${prefix}.vg"
vg convert -g "${prefix}.gfa" -o > "${prefix}.og"

# Save base frequencies
echo "$prefix $totalA $totalC $totalG $totalT" >> soibean_db.baseFreq

echo "Done. You can now use your harvested database. Please use the prefix as --dbprefix input for soibean."
