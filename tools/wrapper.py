#!/usr/bin/env python
import argparse
import gzip
from Bio import SeqIO
import subprocess
import os

def is_fasta(filename):
    try:
        with open(filename, "rt") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)
    except:
        return False

def is_zipped_fasta(filename):
    try:
        with gzip.open(filename, "rt") as handle:
            fasta = SeqIO.parse(handle, "fasta")
            return any(fasta)
    except:
        return False

def is_gam(filename):
    try:
        validate_gam = subprocess.check_output("/tools/src/haplocart-1.0/vgan/dep/vg/bin/vg view -a " + filename + " -j > /dev/null 2> /dev/null", shell=True)
        return True
    except subprocess.CalledProcessError:
        return False

def autodetect(input_file):

    if not os.path.exists(input_file):
        print("NONEXISTENT FILE")
        raise Exception("NONEXISTENT FILE")

    # https://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python
    textchars = bytearray({7,8,9,10,12,13,27} | set(range(0x20, 0x100)) - {0x7f})
    is_binary_string = lambda bytes: bool(bytes.translate(None, textchars))
    is_binary = is_binary_string(open(input_file, 'rb').read(1024))

    if not is_binary:
        # The Vgan team fully supports the rights of the non-binary
        with open(input_file, "rt") as f:
            try:
                records = list(SeqIO.parse(f, "fastq"))
                # Make sure it's bona fide FASTQ
                for record in records:
                    score=record.letter_annotations["phred_quality"]
                if records[0].id.split("/")[0] == records[1].id.split("/")[0]:
                    return "Interleaved"
                else:
                    return "FASTQ"
            except:
                if is_fasta(input_file):
                    return "FASTA"
                else:
                    return "ERROR IN FILE"

    else:
        with gzip.open(input_file, "rt") as g:
            if is_gam(input_file):
                return "GAM"
            else:
                try:
                    records = SeqIO.parse(g, "fastq")
                    for record in records:
                        score=record.letter_annotations["phred_quality"]
                    return "FASTQ"
                except:
                    if is_zipped_fasta(input_file):
                       return "FASTA"
                    else:
                        return "ERROR IN FILE"


def run_haplocart():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f1", "--file1", help="First input file")
    parser.add_argument("-f2", "--file2", type=str)
    parser.add_argument("-j", "--jobid", type=str)
    parser.add_argument("-e", "--bg_error", type=float, default=0.0001)
    parser.add_argument("-p", "--posterior", type=str, help="Compute posterior")
    parser.add_argument("-z", "--tempdir", type=str, help="Temporary directory")
    args = parser.parse_args()
    #print("JOB ID: ", args.jobid)
    #print("BG ERROR: ", args.bg_error)
    #print("Posterior: ", args.posterior)
    #print("FQ1: ", args.file1)
    #print("FQ2: ", args.file2)

    #print("ARGS: ", args)
    #print("tempdir: ", args.tempdir)
    format = autodetect(args.file1)
    #print("Format: ", format)

    argslist=[]
    if format == "FASTA":

        arglist = ["/tools/src/haplocart-1.0/vgan/bin/vgan", "haplocart", "-f", args.file1, "-w", "-t", "1", "-e", str(args.bg_error)]

    elif format == "Interleaved":
        arglist = ["/tools/src/haplocart-1.0/vgan/bin/vgan", "haplocart", "-fq1", args.file1, "-i", "-w", "-t", "-1"]

    elif format == "GAM":
        arglist = ["/tools/src/haplocart-1.0/vgan/bin/vgan", "haplocart", "-g", args.file1, "-w", "-t", "-1"]

    elif format == "FASTQ":
        if args.file2 == None:
            arglist = ["/tools/src/haplocart-1.0/vgan/bin/vgan", "haplocart", "-fq1", args.file1, "-w", "-t", "-1"]
        else:
            format2 = autodetect(args.file2)
            if format2 == "FASTQ":
                arglist = ["/tools/src/haplocart-1.0/vgan/bin/vgan", "haplocart", "-fq1", args.file1, "-fq2", args.file2, "-w", "-t", "-1"]
            else:
                return "ERROR IN FILE 2"

    else:
        print("ERROR IN FILE")
        quit()

    if args.posterior=="yes":
        arglist.append("-p")

    ret = subprocess.run(" ".join(arglist), shell=True, capture_output=True, text=True)
    print(ret.stdout)
    #print(ret.stderr)

run_haplocart()

