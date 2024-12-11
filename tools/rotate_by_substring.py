import argparse
from Bio import SeqIO

def find_longest_common_substring(arr):
    if not arr:
        return ""
    # Initialize the longest common substring to be the first sequence itself
    longest_common_substring = arr[0]
    for sequence in arr:
        temp_substring = ""
        for i in range(len(longest_common_substring)):
            candidate_substring = longest_common_substring[i:]
            if 'N' in candidate_substring:
                continue  # Skip substrings that contain 'N'
            if candidate_substring in sequence:
                temp_substring = candidate_substring
                break
        longest_common_substring = temp_substring
        if not longest_common_substring:
            break  # No common substring exists
    return longest_common_substring


def rotate_sequence(sequence, substring):
    # Find the start of the substring and rotate
    start_index = sequence.find(substring)
    if start_index != -1:
        return sequence[start_index:] + sequence[:start_index]
    return sequence

def process_fasta(input_file, output_prefix):
    sequences = []
    with open(input_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences.append(str(record.seq))

    longest_common_substring = find_longest_common_substring(sequences)
    print("Longest common substring:", longest_common_substring)

    rotated_sequences = []
    with open(input_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            record.seq = rotate_sequence(str(record.seq), longest_common_substring)
            rotated_sequences.append(record)

    output_file = f"{output_prefix}_rotated.fa"
    with open(output_file, "w") as file:
        SeqIO.write(rotated_sequences, file, "fasta")
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rotate sequences in a FASTA file based on the longest common substring.")
    parser.add_argument("input_fasta", help="Input FASTA file path")
    parser.add_argument("output_prefix", help="Output file prefix")
    args = parser.parse_args()

    process_fasta(args.input_fasta, args.output_prefix)
