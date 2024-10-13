from Bio import SeqIO
import sys
import random

def introduce_damage(sequence, proportion):
    sequence = list(sequence)  # Convert to list for mutability
    length = len(sequence)
    
    # Mutate C -> T and G -> A
    for i in range(length):
        if sequence[i] == 'C' and random.random() < proportion:
            sequence[i] = 'T'
        elif sequence[i] == 'G' and random.random() < proportion:
            sequence[i] = 'A'
    
    return ''.join(sequence)

def process_fasta(input_file, proportion):
    output_file = input_file.rsplit('.', 1)[0] + "_damaged.fasta"
    
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            record.seq = introduce_damage(str(record.seq), proportion)
            SeqIO.write(record, outfile, "fasta")

if __name__ == "__main__":
    input_file = sys.argv[1]
    proportion = float(sys.argv[2])
    process_fasta(input_file, proportion)

