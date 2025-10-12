import modules.dna_rna_tools as drt
from modules.fastq_tools import is_valid_gc, is_valid_length, is_valid_quality

import os
from pathlib import Path


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


def filter_fastq(input_fastq: str, output_fastq: str, gc_bounds: int | float | tuple = (0, 100), length_bounds: int | float | tuple = (0, 2**32), quality_threshold: int = 0, overwrite=False) -> str:
    '''
    Filter FASTQ reads by GC%, length, and mean Phred quality; write passing reads to output.

    Args:
        input_fastq: Path to input FASTQ.
        output_fastq: Path to output FASTQ (overwrites).
        gc_bounds: Allowed GC% (max or (min, max)).
        length_bounds: Allowed read length (max or (min, max)).
        quality_threshold: Minimum mean Phred score.

    Returns:
        Path to output_fastq.
    '''

    input_filepath = os.path.abspath(input_fastq)
    output_dir = '/'.join(input_filepath.split('/')[0:-1]) + '/filtered'
    output_filename = output_fastq.split('/')[-1]
    output_filepath = f'{output_dir}/{output_filename}'
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if os.path.exists(output_filepath) and not overwrite:
        print(f'ПРЕДУПРЕЖДЕНИЕ: такой файл уже существует: {output_filepath}')
        return None

    with open(input_fastq, 'r') as input_fastq_file, open(output_filepath, 'w') as output_fastq_file:
        fastq_data = []
        for row in input_fastq_file:
            fastq_data.append(row.strip())
            if len(fastq_data) == 4:
                header, seq, sep, qual = fastq_data
                is_valid_gc_ = is_valid_gc(seq, gc_bounds)
                is_valid_length_ = is_valid_length(seq, length_bounds)
                is_valid_quality_ = is_valid_quality(qual, quality_threshold)
                if is_valid_gc_ and is_valid_length_ and is_valid_quality_:
                    output_fastq_file.write(f'{header}\n{seq}\n{sep}\n{qual}\n')
                fastq_data = []
    return output_filepath
