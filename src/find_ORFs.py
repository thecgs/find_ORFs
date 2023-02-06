#!/usr/bin/env python
# coding: utf-8

import FastaIO
import argparse
import os, re, sys
from CodonTable import getCodonTable

parser = argparse.ArgumentParser(description='Find open Reading Frames (ORFs)', add_help=False, epilog='date:2023/02/05 author:guisen chen email:thecgs001@foxmail.com')
required = parser.add_argument_group('required arguments')
optional = parser.add_argument_group('optional arguments')
required.add_argument('-i', '--input', metavar='[input_file]', help='A file of fasta format', required=True)
optional.add_argument('-o', '--output', metavar='[output_file]', help='A file of output. defualt: stdout.', default=sys.stdout)
optional.add_argument('-outfmt', '--outfmt', metavar='[int]', type=int, choices=[0,1,2], help='Output file format. 0: fasta; 1: tsv; 2: gff. defualt: 0.', default=0)
optional.add_argument('-phase', '--phase', metavar='[int]', type=int, choices=[0,1,2,3], help='Start address of sequence. 0:all; 1: first base; 2: second base; 3: third base. defualt: 0.', default=0)
optional.add_argument('-strand', '--strand', metavar='[int]', type=int, choices=[0,1,2], help='Search strand, 0: both; 1: +; 2: -. defualt: 0.', default=0)
optional.add_argument('-min_len', '--min_len', metavar='[int]', type=int, help='Min Length. defualt: 0.', default=0)
optional.add_argument('-max_len', '--max_len', metavar='[int]', type=int, help='Max Length. defualt: INF.', default=float('inf'))
optional.add_argument('-translate', '--translate', action="store_true", help='Translate DNA to protein sequence. defualt: False.')
optional.add_argument('-codontable', '--codontable', metavar='[int]', type=int, choices=[1,2,3,4,5,6,9,10,11,12,13,14,16,21,22,23,24,25,26,27,28,29,30,31,33],help='Genetic code to use. choices: 1-6, 9-14, 16, 21-31, 33. defualt: 1.', default=1)
optional.add_argument('-stop_codons', '--stop_codons', metavar='[str]', nargs= '+', type=str, help='Stop codons. defualt: None.', default=None)
optional.add_argument('-start_codons', '--start_codons', metavar='[str]', nargs='+', type=str, help='Start codons. defualt: None.', default=None)
optional.add_argument('-start_codon_model', '--start_codon_model', metavar='[int]', type=int, choices=[0,1], help='ORF start codon to use. 0: only "ATG"; 1: "ATG" and alternative initiation codons. defualt: 0.', default=0)
optional.add_argument('-remove_stop_codon', '--remove_stop_codon', action="store_true", help='Remove portein seqence stop codon. defualt: False.')
optional.add_argument('-remove_nested', '--remove_nested', action="store_true", help='Ignore nested ORFs:. defualt: False.')
optional.add_argument('-h', '--help', action='help', help='show this help message and exit')
optional.add_argument('-v', '--version', action='version', version='v1.00')
args = parser.parse_args()

def find_ORFs_Pos(seqence, phase=0, codontable=1, start_codon_model=0, start_codons=None, stop_codons = None, min_len=0, max_len=float('inf'), remove_nested=True):
    if start_codons == None:
        if start_codon_model == 0:
            start_codons = ['ATG']
        else:
            start_codons = getCodonTable(codontable)[1]
    if stop_codons  == None:
        stop_codons = getCodonTable(codontable)[2]
    
    patterns = []
    for start_codon in start_codons:
        for stop_codon in stop_codons:
            patterns.append(start_codon + '(...)*' + stop_codon)
            patterns.append(start_codon + '(...)*?' + stop_codon)
    
    positions = []
    for pattern in patterns:
        positions.extend(list(re.finditer(pattern, seqence)))

    positions = sorted(list({i.span() for i in positions})) #去重并排序
    positions = [i for i in positions if ((i[1] - i[0]) >= min_len and (i[1] - i[0]) <= max_len)] #ORF长度限制
    
    for i in positions:
        for j in positions:
            #去除中间有终止密码子的ORF
            if i[0] == j[0]:
                if i[1] < j[1]:
                    if j in positions:
                        positions.remove(j)
                elif i[1] > j[1]:
                    if i in positions:
                        positions.remove(i)
            #仅保留最长ORF
            if remove_nested==True:
                if i[0] < j[0] and j[1] <= i[1]:
                    if j in positions:
                        positions.remove(j)
                elif j[0] < i[0] and i[1] <= i[1]:
                    if i in positions:
                        positions.remove(i)
            #相位0, 1, 2, 3
            if phase == 1:
                if i[0]%3 != 0 and i in positions:
                    positions.remove(i)
            if phase == 2:
                if i[0]%3 != 1 and i in positions:
                    positions.remove(i)
            if phase == 3:
                if i[0]%3 != 2 and i in positions:
                    positions.remove(i)
            else:
                pass
    
    return positions



def main(input_file, outfmt=0, output_file=sys.stdout, remove_stop_codon=True, strand=0, phase=0, codontable=1, start_codon_model=0, start_codons=None, stop_codons = None, min_len=0, max_len=float('inf'), remove_nested=True, translate=False):
    
    num = 1
    if output_file==sys.stdout:
        out = output_file
    else:
        out = open(output_file,'w')
    
    if outfmt==1:
        out.write('ID\tStrand\tSource\tStart\tEnd\tLength\tSeqence\n')
            
    for record in FastaIO.FastaIO(input_file).parse():
        # +
        if strand == 1 or strand == 0:
            seqence = record.Seq._Seq
            positions = find_ORFs_Pos(seqence, phase, codontable, start_codon_model, start_codons, stop_codons, min_len, max_len, remove_nested)
            for position in positions:
                if translate == True:
                    if remove_stop_codon == True:
                        if outfmt == 0:
                            out.write('>ORF'+ str(num) + '\tLength: ' + str(position[1] - position[0]) + '\t' + 'Strand: +' + '\t' + 'Source: ' + record.Name + ':' + str(position[0]+1) + '-' + str(position[1]) + '\n' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable).remove_StopCodon())+'\n')
                        elif outfmt == 1:
                            out.write('ORF'+ str(num) + '\t+\t' + record.Name + '\t' + str(position[0]+1) + '\t' + str(position[1]) + '\t' + str(position[1] - position[0]) + '\t' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable).remove_StopCodon()) + '\n')
                        elif outfmt == 2:
                            out.write(record.Name + '\tfind_ORFs\tORF\t' + str(position[0]+1) + '\t' + str(position[1]) + '\t.\t+\t.\tID=ORF' + str(num) + '; Length=' + str(position[1] - position[0]) + '; Seqence=' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable).remove_StopCodon()) + ';\n')
                    else:
                        if outfmt == 0:
                            out.write('>ORF'+ str(num) + '\tLength: ' + str(position[1] - position[0]) + '\t' + 'Strand: +' + '\t' + 'Source: ' + record.Name + ':' + str(position[0]+1) + '-' + str(position[1]) + '\n' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable))+'\n')
                        elif outfmt == 1:
                            out.write('ORF'+ str(num) + '\t+\t' + record.Name + '\t' + str(position[0]+1) + '\t' + str(position[1]) + '\t' + str(position[1] - position[0]) + '\t' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable)) + '\n')
                        elif outfmt == 2:
                            out.write(record.Name + '\tfind_ORFs\tORF\t' + str(position[0]+1) + '\t' + str(position[1]) + '\t.\t+\t.\tID=ORF' + str(num) + '; Length=' + str(position[1] - position[0]) + '; Seqence=' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable)) + ';\n')      
                else:
                    if outfmt == 0:
                        out.write('>ORF'+ str(num) + '\tLength: ' + str(position[1] - position[0]) + '\t' + 'Strand: +' + '\t' + 'Source: ' + record.Name + ':' + str(position[0]+1) + '-' + str(position[1]) + '\n' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]))+'\n')
                    elif outfmt == 1:
                        out.write('ORF'+ str(num) + '\t+\t' + record.Name + '\t' + str(position[0]+1) + '\t' + str(position[1]) + '\t' + str(position[1] - position[0]) + '\t' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]])) + '\n')
                    elif outfmt == 2:
                        out.write(record.Name + '\tfind_ORFs\tORF\t' + str(position[0]+1) + '\t' + str(position[1]) + '\t.\t+\t.\tID=ORF' + str(num) + '; Length=' + str(position[1] - position[0]) + '; Seqence=' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]])) + ';\n')     
                num += 1
        # -
        if strand == 2 or strand == 0:
            seqence = str(record.Seq.reverse_compliment())
            positions = find_ORFs_Pos(seqence, phase, codontable, start_codon_model, start_codons, stop_codons, min_len, max_len, remove_nested)
            
            l = len(seqence)
            for position in positions:
                if translate == True:
                    if remove_stop_codon == True:
                        if outfmt == 0:
                            out.write('>ORF'+ str(num) + '\tLength: ' + str(position[1] - position[0]) + '\t' + 'Strand: -' + '\t' + 'Source: ' + record.Name + ':' + str(l - position[1]) + '-' + str(l - position[0]+1) + '\n' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable).remove_StopCodon())+'\n')
                        elif outfmt == 1:
                            out.write('ORF'+ str(num) + '\t-\t' + record.Name + '\t' + str(l - position[1]) + '\t' + str(l - position[0]+1) + '\t' + str(position[1] - position[0]) + '\t' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable).remove_StopCodon()) + '\n')
                        elif outfmt == 2:
                            out.write(record.Name + '\tfind_ORFs\tORF\t' + str(l - position[1]) + '\t' + str(l - position[0]+1) + '\t.\t-\t.\tID=ORF' + str(num) + '; Length=' + str(position[1] - position[0]) + '; Seqence=' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable).remove_StopCodon()) + ';\n')
                    else:
                        if outfmt == 0:
                            out.write('>ORF'+ str(num) + '\tLength: ' + str(position[1] - position[0]) + '\t' + 'Strand: -' + '\t' + 'Source: ' + record.Name + ':' + str(l - position[1]) + '-' + str(l - position[0]+1) + '\n' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable))+'\n')
                        elif outfmt == 1:
                            out.write('ORF'+ str(num) + '\t-\t' + record.Name + '\t' + str(l - position[1]) + '\t' + str(l - position[0]+1) + '\t' + str(position[1] - position[0]) + '\t' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable)) + '\n')
                        elif outfmt == 2:
                            out.write(record.Name + '\tfind_ORFs\tORF\t' + str(l - position[1]) + '\t' + str(l - position[0]+1) + '\t.\t-\t.\tID=ORF' + str(num) + '; Length=' + str(position[1] - position[0]) + '; Seqence=' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]).translate(codontable=codontable)) + ';\n')      
                else:
                    if outfmt == 0:
                        out.write('>ORF'+ str(num) + '\tLength: ' + str(position[1] - position[0]) + '\t' + 'Strand: -' + '\t' + 'Source: ' + record.Name + ':' + str(l - position[1]) + '-' + str(l - position[0]+1) + '\n' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]]))+'\n')
                    elif outfmt == 1:
                        out.write('ORF'+ str(num) + '\t-\t' + record.Name + '\t' + str(l - position[1]) + '\t' + str(l - position[0]+1) + '\t' + str(position[1] - position[0]) + '\t' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]])) + '\n')
                    elif outfmt == 2:
                        out.write(record.Name + '\tfind_ORFs\tORF\t' + str(l - position[1]) + '\t' + str(l - position[0]+1) + '\t.\t-\t.\tID=ORF' + str(num) + '; Length=' + str(position[1] - position[0]) + '; Seqence=' + str(FastaIO.NucleicSeqence(seqence[position[0]:position[1]])) + ';\n')
                num += 1

main(input_file=args.input, \
     outfmt=args.outfmt, \
     output_file=args.output, \
     remove_stop_codon=args.remove_stop_codon, \
     strand=args.strand, \
     phase=args.phase, \
     codontable=args.codontable, \
     start_codon_model=args.start_codon_model, \
     start_codons=args.start_codons, \
     stop_codons=args.stop_codons, \
     min_len=args.min_len, \
     max_len=args.max_len, \
     remove_nested=args.remove_nested, \
     translate=args.translate)