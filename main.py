import modules.dna_rna_tools as drt
from modules.fastq_tools import is_valid_gc, is_valid_length, is_valid_quality


# HW 3 https://github.com/Python-BI-2025-26/hw3-functions-CaptainLabMan
# HW 4 https://github.com/Python-BI-2025-26/course_materials/blob/main/homeworks/HW4_Modules/HW4_Modules.md


processes = {
    'is_nucleic_acid': drt.is_nucleic_acid,
    'transcribe': drt.transcribe,
    'reverse': drt.reverse,
    'complement': drt.complement,
    'reverse_complement': drt.reverse_complement,
}


def run_dna_rna_tools(*args: str) -> str | list:
    '''
    Performs validation and various operations on a nucleotide sequence (DNA or RNA).

    Available operations:
     - is_nucleic_acid
     - transcribe
     - reverse
     - complement
     - reverse_complement

     Arguments:
     *args: str

     Returns: str | list
    '''

    *seqs, proc = args
    single_seq = len(args) == 2
    results = []
    for seq in seqs:
        if proc == 'is_nucleic_acid':
            result = processes['is_nucleic_acid'](seq)
        elif not processes['is_nucleic_acid'](seq):
            result = None
        else:
            result = processes[proc](seq)

        if single_seq:
            return result
        results.append(result)
    return results


def filter_fastq(seqs: dict, gc_bounds: int | float | tuple = (0, 100), length_bounds: int | float | tuple = (0, 2**32), quality_threshold: int = 0) -> dict:
    '''
    Performs evaluation of various parameters of reads from a FASTQ file.

    Available operations:
     - is_valid_gc
     - length_bounds
     - is_valid_quality

    Arguments:
    seqs: dict
    gc_bounds: int | float | tuple
    length_bounds: int | float | tuple
    quality_threshold: int

    Returns dict
    '''

    output = {}
    for seq_id in seqs:
        seq, qual = seqs[seq_id]
        is_valid_gc_ = is_valid_gc(seq, gc_bounds)
        is_valid_length_ = is_valid_length(seq, length_bounds)
        is_valid_quality_ = is_valid_quality(qual, quality_threshold)
        if is_valid_gc_ and is_valid_length_ and is_valid_quality_:
            output[seq_id] = seqs[seq_id]
    return output
