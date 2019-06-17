import argparse
import os.path
import sys

parser = argparse.ArgumentParser(description='Tabulate barcode counts across samples')
parser.add_argument('samplefiles', metavar='SAMPLE.TXT', nargs='+', help='Input barcode count (tab-delimited, barcode<TAB>count)')
parser.add_argument('-o', '--output', required=True, help='Output filename for tab-delimited count table')
parser.add_argument('-m', '--mintotal', type=int, help='Minimum total read count across all samples to report barcode')
parser.add_argument('-n', '--minsamples', type=int, help='Minimum number of samples with barcode to report')
parser.add_argument('--mininsample', type=int, help='Minimum read count in at least one sample to report barcode')
parser.add_argument('--omitfile', help='Output filename for omitted barcodes')
args = parser.parse_args()

(root, ext) = os.path.splitext(args.output)
omitfile = '{}-omitted{}'.format(root, ext) if args.omitfile is None else args.omitfile

with open(args.output, 'w') as reportout:
    with open(omitfile, 'w') as omitout:
        barcode_sample_counts = {}
        samples = []
        for samplefile in args.samplefiles:
            sample = os.path.basename(samplefile)
            samples.append(sample)
            for l in open(samplefile, 'r'):
                (barcode, countstr) = l.strip().split('\t')
                count = int(countstr)
                barcode_samples = barcode_sample_counts.get(barcode, {})
                barcode_samples[sample] = count
                barcode_sample_counts[barcode] = barcode_samples

        for out in [reportout, omitout]:
            out.write('barcode')
            for sample in samples:
                out.write('\t{}'.format(sample))
            out.write('\n')
                
        for (barcode, sample_counts) in barcode_sample_counts.items():
            keep_total = ( (args.mintotal is None) or (sum(sample_counts.value() >= args.mintotal)) )
            keep_insample = ( (args.minsamples is None) or (len(sample_counts) >= args.minsamples) )
            keep_minin = ( (args.mininsample is None) or (max(sample_counts.values()) >= args.mininsample) )
            keep = keep_total and keep_insample and keep_minin
            out = reportout if keep else omitout
            out.write('{}'.format(barcode))
            for sample in samples:
                out.write('\t{}'.format(sample_counts.get(sample, 0)))
            out.write('\n')
