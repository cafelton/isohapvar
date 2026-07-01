#! /usr/bin/env python3

import argparse
import pysam
import vcfpy
import subprocess
from collections import namedtuple

Exon = namedtuple("Exon", ("start", "end"))
Isoform = namedtuple("Isoform", ("chrom", "start", "end", "name", "strand", "exons"))
Variant = namedtuple("Variant", ("pos", "ref", "alt"))

COMPBASE = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
            'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S',
            'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf', required=True, type=str,
                        help="vcf file with haplotyped variants")
    parser.add_argument('-o', '--output', required=True,
                        help="output prefix")
    parser.add_argument('--genome',
                        type=str, required=True,
                        help='genome fasta file')
    parser.add_argument('--isoform_bed', required=True,
                        help='bed12 file of isoforms')
    args = parser.parse_args()
    return args

def get_vars_overlapping_block(block, vars):
    block_vars = []
    for var in vars:
        if block.start <= var.pos <= block.end:
            block_vars.append(var)
    block_vars.sort(key=lambda x: x.pos, reverse=True)
    return block_vars

def check_vcf_line(record):
    if len(record.calls) > 1:
        raise ValueError('Genotypes for multiple samples detected. Please use a vcf with variants called on a single sample')
    if 'GT' not in record.calls[0].data:
        raise ValueError('Genotype not found in vcf. Please use a haplotyped vcf.')

def get_var_haps(alt_index, call):
    var_haps = []
    for hap in range(2):
        if call.gt_alleles[hap] == alt_index:
            var_haps.append(hap)
    return var_haps

def write_alts_with_haps(record, call, out):
    outline = [record.CHROM, record.POS - 1, record.POS]
    for i, alt in enumerate(record.ALT):
        var_haps = get_var_haps(i + 1, call)
        if len(var_haps) > 0:
            outname = [record.REF, alt.value, ','.join([str(x) for x in var_haps])]
            out.write('\t'.join([str(x) for x in outline] + [':'.join([str(x) for x in outname])]) + '\n')

def process_vcf(vcf, temp_prefix):
    with vcfpy.Reader.from_path(vcf) as reader, open(temp_prefix + '.vars.bed', 'w') as out:
        reader = vcfpy.Reader.from_path(vcf)
        for record in reader:
            check_vcf_line(record)
            call = record.calls[0]
            if call.is_variant and (call.is_phased or call.gt_alleles[0] == call.gt_alleles[1]):  # is phased or homozygous
                write_alts_with_haps(record, call, out)        

def get_add_isoform(line, transcript_to_hap_to_vars):
    chrom, name, strand, start, end, esizes, estarts = line[0], line[3], line[5], int(line[1]), int(line[2]), [int(x) for x in line[10].rstrip(',').split(',')], [int(x) for x in line[11].rstrip(',').split(',')]
    exons = tuple([Exon(start + estarts[i], start + estarts[i] + esizes[i]) for i in range(len(esizes))])
    this_iso = Isoform(chrom, start, end, name, strand, exons)
    if this_iso not in transcript_to_hap_to_vars:
        transcript_to_hap_to_vars[this_iso] = {0: [], 1: []}
    return this_iso

def add_var_to_iso(var_info, this_iso, transcript_to_hap_to_vars):
    var_pos = int(var_info[1])
    ref, alt, haps = var_info[3].split(':')
    for hap in haps.split(','):
        transcript_to_hap_to_vars[this_iso][int(hap)].append(Variant(var_pos, ref, alt))

def get_transcript_to_hap_to_vars(temp_prefix):
    transcript_to_hap_to_vars = {}
    for line in open(temp_prefix + '.isovars.txt'):
        line = line.rstrip('\n').split('\t')
        var_info = line[-5:]
        this_iso = get_add_isoform(line, transcript_to_hap_to_vars)
        if var_info[0] != '.':
            add_var_to_iso(var_info, this_iso, transcript_to_hap_to_vars)
    return transcript_to_hap_to_vars

def check_rev_comp(iso_seq, strand):
    if strand == '-':
        iso_seq = iso_seq[::-1]
        for i in range(len(iso_seq)):
            iso_seq[i] = COMPBASE[iso_seq[i]]
    return iso_seq

def adjust_block_seq(block, chrom, block_vars, genome):
    block_seq = list(genome.fetch(chrom, block.start, block.end).upper())
    for var in block_vars:
        before, after = block_seq[:var.pos - block.start], block_seq[(var.pos - block.start) + len(var.ref):]
        block_seq = before + list(var.alt) + after
    return block_seq

def get_iso_seq_with_vars(isoform, genome, hap_vars):
    iso_seq, iso_vars = [], []
    for block in isoform.exons:
        block_vars = get_vars_overlapping_block(block, hap_vars)
        block_seq = adjust_block_seq(block, isoform.chrom, block_vars, genome)
        iso_seq.extend(block_seq)
        iso_vars.extend(block_vars)
    iso_seq = check_rev_comp(iso_seq, isoform.strand)
    return iso_seq, iso_vars

def write_seq_with_vars(transcript_to_hap_to_vars, genome_file, output):
    with open(output + '.isoswithvars.fa', 'w') as out, pysam.FastaFile(genome_file) as genome:
        for isoform in transcript_to_hap_to_vars:
            haps = [transcript_to_hap_to_vars[isoform][x] for x in range(2)]
            # if haplotypes are identical, just report 1 hap
            if haps[0] == haps[1]:
                haps = [haps[0]]
            for h, hap_vars in enumerate(haps):
                iso_seq, iso_vars = get_iso_seq_with_vars(isoform, genome, hap_vars)
                var_strings = [f'{x.pos};{x.ref};{x.alt}' for x in iso_vars]
                out.write(f'>{isoform.name}::{h} {",".join(var_strings)}\n{"".join(iso_seq)}\n')

def addvariants():
    """
    Filter vcf and convert to bed
    use bedtools to overlap bed of variants with bed of transcripts
    for each transcript, separate variants by haplotype
    Generate haplotype-specific transcript
    """
    args = parse_args()
    temp_prefix = 'temp'
    process_vcf(args.vcf, temp_prefix)
    subprocess.check_call(['bedtools', 'intersect', '-split', '-wao', '-a', args.isoform_bed, '-b', temp_prefix + '.vars.bed'], stdout=open(temp_prefix + '.isovars.txt', 'w'))
    transcript_to_hap_to_vars = get_transcript_to_hap_to_vars(temp_prefix)
    write_seq_with_vars(transcript_to_hap_to_vars, args.genome, args.output)
    subprocess.check_call(['rm', temp_prefix + '.vars.bed', temp_prefix + '.isovars.txt'])


if __name__ == "__main__":
    addvariants()
