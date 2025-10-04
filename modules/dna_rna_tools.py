def is_nucleic_acid(seq):
    if not seq:
        return False
    seq = set(seq.upper())
    dna_set = set('ATGC')
    rna_set = set('AUGC')
    return seq <= dna_set or seq <= rna_set


def transcribe(seq):
    return seq.replace('T', 'U').replace('t', 'u')


def reverse(seq):
    return seq[::-1]


def complement(seq):
    dna_nucs = {
        'T': 'A',
        'A': 'T',
        'C': 'G',
        'G': 'C',
        'U': 'A',
        't': 'a',
        'a': 't',
        'c': 'g',
        'g': 'c',
        'u': 'a',
    }

    comp_seq = ''
    for nuc in seq:
        comp_seq += dna_nucs[nuc]
    if 'U' in seq.upper():
        comp_seq = comp_seq.replace('T', 'U').replace('t', 'u')
    return comp_seq


def reverse_complement(seq):
    return complement(seq)[::-1]
