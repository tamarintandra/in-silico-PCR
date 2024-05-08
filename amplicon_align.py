#!/usr/bin/env python3

import magnumopus
import argparse

def main():
    parser = argparse.ArgumentParser(description="Perform in-silico PCR on two assemblies and align the amplicons")

    parser.add_argument("-1", dest="ASSEMBLY1", required=True, help="Path to the first assembly file")
    parser.add_argument("-2", dest="ASSEMBLY2", required=True, help="Path to the second assembly file")
    parser.add_argument("-p", dest="PRIMERS", required=True, help="Path to the primer file")
    parser.add_argument("-m", dest="MAX_AMPLICON_SIZE", required=True, type=int, help="Maximum amplicon size for isPCR")
    parser.add_argument("--match", dest="match", required=True, type=int, help="Match score to use in alignment")
    parser.add_argument("--mismatch", dest="mismatch", required=True, type=int, help="Mismatch penalty to use in alignment")
    parser.add_argument("--gap", dest="gap", required=True, type=int, help="Gap penalty to use in alignment")

    args = parser.parse_args()

    #Perform isPCR on both assembly files
    amplicon1 = magnumopus.ispcr(args.PRIMERS, args.ASSEMBLY1, args.MAX_AMPLICON_SIZE)
    amplicon2 = magnumopus.ispcr(args.PRIMERS, args.ASSEMBLY2, args.MAX_AMPLICON_SIZE)

    #Check the best alignment orientation
    amplicon1 = amplicon1.split('\n', 1)[1].strip()
    amplicon2 = amplicon2.split('\n', 1)[1].strip()
    rc_amplicon2 = reverse_complement(amplicon2)
    alignment1, score1 = magnumopus.needleman_wunsch(amplicon1, amplicon2, args.match, args.mismatch, args.gap)
    alignment2, score2 = magnumopus.needleman_wunsch(amplicon1, rc_amplicon2, args.match, args.mismatch, args.gap)

    #Choose the best alignment based on score
    if score1 >= score2:
        best_alignment = alignment1
        best_score = score1
    else:
        best_alignment = alignment2
        best_score = score2
    
    #Print the alignment and alignment score
    print(best_alignment[0])
    print(best_alignment[1])
    print(best_score)

def reverse_complement(sequence):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    reversed_sequence = sequence[::-1]
    complemented_sequence = ''.join(complement_dict[base] for base in reversed_sequence)

    return complemented_sequence

if __name__ == "__main__":
    main()