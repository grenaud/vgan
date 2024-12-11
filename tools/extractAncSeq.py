import sys

def parse_state_file_to_fasta(input_file_path):
    sequences = {}
    with open(input_file_path, 'r') as file:
        for line in file:
            if line.startswith('#') or not line.strip():
                continue  # Skip the header or empty lines
            parts = line.strip().split()
            # Skip lines that do not have at least 3 parts or where the second part is not a digit
            if len(parts) < 3 or not parts[1].isdigit():
                continue
            node, site, state = parts[0], int(parts[1]), parts[2]
            if node not in sequences:
                sequences[node] = []
            sequences[node].append(state)

    # Convert sequences dictionary to FASTA format
    fasta_content = ""
    for node, seq_list in sequences.items():
        fasta_content += f">{node}\n{''.join(seq_list)}\n"

    return fasta_content

def save_fasta_content_to_file(fasta_content, output_file_name):
    with open(output_file_name, 'w') as file:
        file.write(fasta_content)

def main():
    if len(sys.argv) != 3:
        print("Usage: python extractAncSeq.py <input_state_file> <output_fasta_file>")
        sys.exit(1)

    input_file_path = sys.argv[1]
    output_file_name = sys.argv[2]

    fasta_content = parse_state_file_to_fasta(input_file_path)
    save_fasta_content_to_file(fasta_content, output_file_name)

    print(f"FASTA sequences have been saved to '{output_file_name}'")

if __name__ == "__main__":
    main()

