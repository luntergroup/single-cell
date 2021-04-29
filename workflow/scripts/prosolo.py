#!/usr/bin/env python3

import argparse
import csv
import concurrent.futures
import subprocess as sp
from pathlib import Path

def read_bed_regions(bed_filename):
    with bed_filename.open() as bed:
        bedreader = csv.reader(bed, delimiter='\t')
        res = []
        for row in bedreader:
            res.append(row[0] + ':' + str(int(row[1]) + 1) + '-' + row[2])
        return res

def split_region(region, window_length=3000000):
    contig, region = region.split(":")
    start, end = region.split("-")
    start, end = int(start), int(end)
    result = []
    for window_start in range(start, end, window_length):
        window_end = min(window_start + window_length, end)
        result.append(contig + ":" + str(window_start) + "-" + str(window_end))
    return result

def generate_candidates(reference, cell_bam, bulk_bam, region, output):
    mpileup = sp.Popen(["bcftools", "mpileup", "-f", reference, "--regions", region, "--annotate", "FORMAT/DP", "--no-BAQ", "-Ou", cell_bam, bulk_bam], stdout=sp.PIPE)
    bcftools_view = sp.Popen(["bcftools", "view", "-i", "N_ALT>1", "-Ou"], stdin=mpileup.stdout, stdout=sp.PIPE)
    bcftools_norm = sp.Popen(["bcftools", "norm", "-m", "-any", "-Ou"], stdin=bcftools_view.stdout, stdout=sp.PIPE)
    report = sp.check_output(["bcftools", "view", "-e", 'ALT==\"<*>\"', "-Ob", "-o", output], stdin=bcftools_norm.stdout)
    mpileup.wait()
    
def prosolo_helper(reference, cell_bam, bulk_bam, region, output):
    candidates = output.with_suffix(".candidates.bcf")
    generate_candidates(reference, cell_bam, bulk_bam, region, candidates)
    sp.call([
        "prosolo", \
        "single-cell-bulk", \
        "--omit-indels", \
        "--candidates", candidates, \
        "--output", output, \
        cell_bam, \
        bulk_bam, \
        reference
    ])
    candidates.unlink()

def main(args):
    regions = read_bed_regions(args.regions)
    tmp_dir = args.output.with_suffix('.tmp')
    tmp_dir.mkdir(exist_ok=True)
    tmp_bcfs = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = []
        for region in regions:
            for window in split_region(region):
                tmp_bcf = tmp_dir / (window.replace(':', '_') + '.bcf')
                futures.append(executor.submit(prosolo_helper, reference=args.reference, cell_bam=args.cell, bulk_bam=args.bulk, region=window, output=tmp_bcf))
                tmp_bcfs.append(tmp_bcf)
        for future in concurrent.futures.as_completed(futures):
            print(future.result())
    bcftools_concat_cmd = ['bcftools', 'concat', '-Ob', '-o', args.output] + tmp_bcfs
    print(" ".join(str(c) for c in bcftools_concat_cmd))
    sp.call(bcftools_concat_cmd)
    bcftools_index_cmd = ['bcftools', 'index', args.output]
    print(" ".join(str(c) for c in bcftools_index_cmd))
    sp.call(bcftools_index_cmd)
    for tmp_bcf in tmp_bcfs:
        tmp_bcf.unlink()
    tmp_dir.rmdir()

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference',
                        type=Path,
                        required=True,
                        help='Reference FASTA')
    parser.add_argument('--cell',
                        type=Path,
                        required=True,
                        help='Cell BAM')
    parser.add_argument('--bulk',
                        type=Path,
                        required=True,
                        help='Bulk BAM')
    parser.add_argument('--regions',
                        type=Path,
                        required=True,
                        help='Regions BED')
    parser.add_argument('--output',
                        type=Path,
                        required=True,
                        help='Output VCF')
    parser.add_argument('--threads',
                        type=int,
                        default=1,
                        help='Threads')
    parsed, unparsed = parser.parse_known_args()
    main(parsed)
