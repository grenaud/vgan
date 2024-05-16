import sys
import pysam

def add_variants_to_fasta(vcf_filename, fasta_filename, chromosome_name):
    vcf_file = pysam.VariantFile(vcf_filename)
    fasta_file = pysam.Fastafile(fasta_filename)

    # Check if there are multiple individuals in the VCF
    if len(vcf_file.header.samples) != 1:
        sys.stderr.write("Error: The VCF file contains more than one individual.\n")
        sys.exit(1)
    else:
        individual_id = vcf_file.header.samples[0]
    
    # Check if there are any records for this chromosome
    records = list(vcf_file.fetch(chromosome_name))
    if len(records) == 0:
        sys.stderr.write("Error: There are no variants for this individual on this chromosome.\n")
        sys.exit(1)
    
    # Create a copy of the fasta sequence for the chromosome
    fasta_sequence = fasta_file.fetch(chromosome_name)

    # Iterate over VCF records
    for record in records:
        # Get position (0-based) and alternate alleles
        pos = record.pos - 1  # Adjust for 1-based VCF indexing
        alts = record.alts

        # Check the genotype
        genotype = record.samples[individual_id]['GT']
        if len(set(genotype)) != 1 or genotype[0] not in range(len(alts)+1):
            sys.stderr.write(f"Error: Unexpected genotype {genotype} at position {pos+1}.\n")
            sys.exit(1)

        # If the genotype is not homozygous reference, proceed to modify the sequence
        #sys.stderr.write(f"{pos} {genotype}\n")

      
        # If the genotype is not homozygous reference, proceed to modify the sequence
        if genotype != (0, ):
            # Fetch the reference segment for this variant
            ref_segment = fasta_file.fetch(chromosome_name, pos, pos+len(record.ref)).upper()

            # Check if the reference in the VCF and FASTA match
            if ref_segment != record.ref:
                sys.stderr.write(f"Error: The reference in the VCF ({record.ref}) and the FASTA ({ref_segment}) at position {pos+1} do not match.\n")
                sys.exit(1)

            # Replace reference base(s) with selected alternate allele
            alt = alts[genotype[0] - 1]
            if len(record.ref) == len(alt):  # SNV or MNP
                fasta_sequence = fasta_sequence[:pos] + alt + fasta_sequence[pos+len(alt):]
            elif len(record.ref) < len(alt):  # insertion
                fasta_sequence = fasta_sequence[:pos] + alt + fasta_sequence[pos+len(record.ref):]
            elif len(record.ref) > len(alt):  # deletion
                fasta_sequence = fasta_sequence[:pos] + alt + fasta_sequence[pos+len(record.ref):]
            
            sys.stderr.write(f"Warning: Variant added at position {pos+1}. Reference {record.ref} replaced with {alt}.\n")

    # Write modified sequence to stdout with the individual's ID as the header
    sys.stdout.write(">" + individual_id + "\n")
    sys.stdout.write(fasta_sequence + "\n")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py [vcf_file] [fasta_file] [chromosome]")
        sys.exit(1)

    vcf_filename = sys.argv[1]
    fasta_filename = sys.argv[2]
    chromosome_name = sys.argv[3]

    add_variants_to_fasta(vcf_filename, fasta_filename, chromosome_name)
