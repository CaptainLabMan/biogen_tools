from collections import Counter


def is_valid_gc(seq: str, gc_bounds: int | float | tuple | list) -> bool:
    '''
    Performs a GC content check of the nucleotide sequence.

    Arguments:
     seq: str
     gc_bounds: int | tuple | list

     Returns bool
    '''

    nucs_count = Counter(seq)
    seq_gc = (nucs_count['G'] + nucs_count['C']) / len(seq) * 100
    if isinstance(gc_bounds, (int, float)):
        return seq_gc <= gc_bounds
    return gc_bounds[0] <= seq_gc <= gc_bounds[1]


def is_valid_length(seq: str, length_bounds: int | tuple | list) -> bool:
    '''
    Performs a length check of the nucleotide sequence.

    Arguments:
     seq: str
     gc_bounds: int | tuple | list

     Returns bool
    '''

    if isinstance(length_bounds, (int, float)):
        return len(seq) <= length_bounds
    return length_bounds[0] <= len(seq) <= length_bounds[1]


def is_valid_quality(quality: str, quality_threshold: int) -> bool:
    '''
    Performs a quality check of the nucleotide sequence.

    Arguments:
     seq: str
     quality_threshold: int

     Returns bool
    '''

    nuc_quals = [ord(x) - 33 for x in quality]
    quals_mean = sum(nuc_quals) / len(nuc_quals)
    return quals_mean >= quality_threshold
