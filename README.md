# *In-silico*-PCR and Needleman-Wunsch Alignment
This Python module performs in-silico PCR and Needleman-Wunsch alignment sequentially. Some example data files are provided in the **data** folder.

## *In-Silico* PCR
*is*PCR is a computational simulation of the Polymerase Chain Reaction (PCR), a widely used technique in laboratory settings to amplify target DNA sequences. It employs similar principles as traditional PCR but leverages computational methods to model the chemical processes inherent in PCR, thereby simulating its outcomes.

```ispcr.py``` script employs three distinct steps. The stepwise development process involves:
1. Identify locations where primers would anneal to the target sequence
2. Identify pairs of locations where two primers anneal close enough together and in the correct orientation for amplification to occur
3. Extract the amplified sequence

## Needleman-Wunsch Global Sequence Alignment
The Needleman-Wunsch algorithm uses dynamic programming to find optimal alignment between two sequences, which may or may not contain gaps. The ```nw.py``` script returns one of the optimal alignments as well as its alignment score in the output.

## Integration
The ```amplicon_align.py``` script performs *in-silico* PCR on two assemblies and align the amplicons using the needleman-wunsch algorithm. It takes two assembly files, a primer file, a maximum amplicon size, and a match, mismatch, and gap score as command-line input.

```
An example usage of the script is:

    ./amplicon_align.py -h

usage: amplicon_align.py [-h] -1 ASSEMBLY1 -2 ASSEMBLY2 -p PRIMERS -m MAX_AMPLICON_SIZE
--match MATCH --mismatch MISMATCH --gap GAP

Perform in-silico PCR on two assemblies and align the amplicons

options:
    -h, --help            show this help message and exit
    -1 ASSEMBLY1          Path to the first assembly file
    -2 ASSEMBLY2          Path to the second assembly file
    -p PRIMERS            Path to the primer file
    -m MAX_AMPLICON_SIZE  maximum amplicon size for isPCR
    --match MATCH         match score to use in alignment
    --mismatch MISMATCH   mismatch penalty to use in alignment
    --gap GAP             gap penalty to use in alignment
```

**Note**: Providing negative numbers as command-line inputs can lead argparse to interpret them as options. You can work around this in two ways: 
1) Enclose the number with a leading space (```--gap ' -1'```)
2) Link an option and its value using the "=" symbol (```--gap=-1```)
