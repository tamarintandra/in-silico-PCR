import subprocess
from typing import List, Tuple
from tempfile import NamedTemporaryFile

def step_one(primer_file: str, assembly_file: str) -> List[List[str]]:
    """Identify primer annealing location

    Args:
        primer_file: Path to the primer file containing the primer sequences
        assembly_file: Path to assembly file containing the template sequences
    
    Return:
        Filtered BLAST hits sorted by the start position of alignment
    """
    cmd = [
        "blastn",
        "-task", "blastn-short",
        "-query", primer_file,
        "-subject", assembly_file,
        "-outfmt", "6 std qlen"
    ]
    result = subprocess.run(
        cmd, 
        capture_output=True,
        text=True,
        shell=False,
    )
    filtered_results = []
    for line in result.stdout.strip().split('\n'):
        fields = line.strip().split('\t')
        percent_identity = float(fields[2])
        qlen = int(fields[12])
        if percent_identity >= 80.0 and qlen == int(fields[3]):
            filtered_results.append(fields)
    
    sorted_results = sorted(filtered_results, key=lambda x: (int(x[8])))
    
    return sorted_results

def step_two(sorted_hits: List[str], max_amplicon_size: int) -> List[Tuple[List[str]]]:
    """Identify pairs of primer annealing sites which would yield an amplicon

    Args:
        sorted_hits: sorted BLAST hits from step one, indicating the primer annealing location
        max_amplicon_size: maximum length of amplicon
    
    Return:
        Pairs of two primers that can make an amplicon
    """
    pairs = []
    for i in range(len(sorted_hits)):
        for j in range(i+1, len(sorted_hits)):
            hit1, hit2 = sorted_hits[i], sorted_hits[j]
            start1, end1 = int(hit1[8]), int(hit1[9])
            start2, end2 = int(hit2[8]), int(hit2[9])
            distance = abs(end2-end1)
            if distance <= max_amplicon_size and (
                (start1 < end1 and start2 > end2 and end1 < end2) or
                (start1 > end1 and start2 < end2 and end2 < end1)
                ):
                pairs.append((hit1,hit2))
    return pairs

def step_three(hit_pairs: List[Tuple[List[str]]], assembly_file: str) -> str:
    """Extract amplified sequences

    Args:
        hit_pairs: pairs of primer that would yield an amplicon identified from step two
        assembly_file: path to the assembly file which contains amplicon sequences to be extracted
    
    Return:
        Extracted amplicon sequences
    """
    with NamedTemporaryFile(mode='w', delete=False) as bed_file:
        for hit_pair in hit_pairs:
            contig, start, stop = hit_pair[0][1], int(hit_pair[0][9]), int(hit_pair[1][9])-1
            bed_file.write(f"{contig}\t{start}\t{stop}\n")

    cmd = [
        "seqtk",
        "subseq",
        assembly_file,
        bed_file.name
    ]
    result =  subprocess.run(
        cmd,
        capture_output=True,
        text=True
    )
    amplicon_sequences = result.stdout
    return amplicon_sequences

def ispcr(primer_file: str, assembly_file: str, max_amplicon_size: int) -> str:
    """Simulation of the Polymerase Chain Reaction (PCR) to amplify target DNA sequences

    Args: 
        primer file: path to the primer file
        assembly file: path to the assembly file
        max_amplicon_size: maximum length of amplicon
    
    Return:
        Amplicon sequences between the two primers, produced by the PCR reaction
    """
    sorted_results = step_one(primer_file, assembly_file)
    pairs = step_two(sorted_results, max_amplicon_size)
    amplicon_sequences = step_three(pairs, assembly_file)

    return amplicon_sequences