import argparse
import sys

from Bio import SeqIO

parser = argparse.ArgumentParser(description='Count barcode reads from trimmed FastQ file')
parser.add_argument('-i', '--input', help='Input filename of trimmed barcode FastQ')
parser.add_argument('-o', '--output', help='Output filename for tab-delimited text count')
parser.add_argument('-l', '--length', type=int, help='Length of barcode sequence')
args = parser.parse_args()

with (sys.stdout if args.output is None else open(args.output, 'w')) as hout:
    barcode_counts = {}
    total = 0
    bad_len = 0

    for record in SeqIO.parse(sys.stdin if args.input is None else args.input, "fastq"):
        total += 1
        sequ = str(record.seq)
        if (args.length is None) or len(sequ) == args.length:
            barcode_counts[sequ] = barcode_counts.get(sequ, 0) + 1
        else:
            bad_len += 1
    for barcode,count in barcode_counts.items():
        hout.write("{}\t{}\n".format(barcode, count))
    sys.stderr.write("Processed {:9d} barcodes\n".format(total))
    if args.length is not None:
        sys.stderr.write("          {:9d} were not {:d} long, skipped\n".format(bad_len, args.length))
