import argparse
from Bio import SeqIO
from Bio.Seq import Seq

def find_longest_common_substring(arr):
    if not arr:
        return ""
    longest_common_substring = arr[0]
    for sequence in arr:
        temp_substring = ""
        for i in range(len(longest_common_substring)):
            candidate_substring = longest_common_substring[i:]
            if 'N' in candidate_substring:
                continue  # Skip substrings containing 'N'
            if candidate_substring in sequence:
                temp_substring = candidate_substring
                break
        longest_common_substring = temp_substring
        if not longest_common_substring:
            break
    return longest_common_substring

def rotate_sequence(sequence, substring):
    # Convert both to uppercase for case-insensitive matching
    seq_upper = sequence.upper()
    substring_upper = substring.upper()
    start_index = seq_upper.find(substring_upper)
    if start_index != -1:
        return sequence[start_index:] + sequence[:start_index]
    return None  # substring not found

def process_fasta(input_file, output_prefix, user_substring=None):
    sequences = []
    with open(input_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            sequences.append(str(record.seq))

    if user_substring:
        common_substring = user_substring
        print(f"Using user-provided substring: {common_substring}")
    else:
        common_substring = find_longest_common_substring(sequences)
        if not common_substring:
            print("No common substring found among sequences. No rotation performed.")
            return
        print(f"Longest common substring found: {common_substring}")

    rotated_sequences = []
    unrotated_sequences = []

    with open(input_file, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            rotated = rotate_sequence(str(record.seq), common_substring)
            if rotated is not None:
                record.seq = Seq(rotated)
                rotated_sequences.append(record)
            else:
                print(f"Warning: Substring not found in {record.id}, sequence not rotated.")
                unrotated_sequences.append(record)

    # Write rotated sequences
    output_file = f"{output_prefix}_rotated.fa"
    with open(output_file, "w") as file:
        SeqIO.write(rotated_sequences, file, "fasta")
    print(f"Rotated sequences written to {output_file}")

    # Write unrotated sequences if any
    if unrotated_sequences:
        unrotated_file = f"{output_prefix}_unrotated.fa"
        with open(unrotated_file, "w") as file:
            SeqIO.write(unrotated_sequences, file, "fasta")
        print(f"Unrotated sequences written to {unrotated_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rotate sequences in a FASTA file based on a common substring.")
    parser.add_argument("input_fasta", help="Input FASTA file path")
    parser.add_argument("output_prefix", help="Output file prefix")
    parser.add_argument("--substring", help="Optional common substring to use for rotation", default=None)
    args = parser.parse_args()

    process_fasta(args.input_fasta, args.output_prefix, args.substring)
