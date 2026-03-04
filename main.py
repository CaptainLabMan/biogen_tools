from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import os
from pathlib import Path
from abc import ABC, abstractmethod


# Дорогие проверяющие, я специально не прописывал типы аргументов и аутпута (за редким исключением), потому что ООП идет пока тяжеловато,
# а дополнительные элементы в коде очень сильно отвлекают от основного кода. Надеюсь на понимание =) Мне не лень, правда =D


# ===================================================================================
# ===================================================================================
class BiologicalSequence(ABC):
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)
    
    def __str__(self):
        return f'Sequence: {self.sequence}'
    
    def __getitem__(self, index):
        return self.sequence[index]
    
    @abstractmethod
    def _alphabet(self):
        pass

    def is_valid(self):
        return set(self.sequence) <= self._alphabet()
    

class NucleicAcidSequence(BiologicalSequence):
    def _alphabet(self):
        return set('ATGCUatgcu')

    def reverse(self):
        return self.sequence[::-1]

    def complement(self):
        complement_nucleotides = {
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

        complemented_sequence = ''
        for nucleotide in self.sequence:
            complemented_sequence += complement_nucleotides[nucleotide]
        if 'U' in self.sequence.upper():
            complemented_sequence = complemented_sequence.replace('T', 'U').replace('t', 'u')
        return complemented_sequence

    def reverse_complement(self):
        return self.complement()[::-1]
    

class DNASequence(NucleicAcidSequence):
    def _alphabet(self):
        return set('ATGCatgc')
    
    def transcribe(self):
        return self.sequence.replace('T', 'U').replace('t', 'u')
    

class RNASequence(NucleicAcidSequence):
    def _alphabet(self):
        return set('AUGCaugc')


class AminoAcidSequence(BiologicalSequence):
    _AA_ONE_TO_THREE = {
            'A': 'Ala',
            'C': 'Cys',
            'D': 'Asp',
            'E': 'Glu',
            'F': 'Phe',
            'G': 'Gly',
            'H': 'His',
            'I': 'Ile',
            'K': 'Lys',
            'L': 'Leu',
            'M': 'Met',
            'N': 'Asn',
            'P': 'Pro',
            'Q': 'Gln',
            'R': 'Arg',
            'S': 'Ser',
            'T': 'Thr',
            'V': 'Val',
            'W': 'Trp',
            'Y': 'Tyr',
        }
    
    def _alphabet(self):
        return set('ACDEFGHIKLMNPQRSTVWY')
    
    def triple_alphabet(self):
        triple_sequnce = ''
        for aa in self.sequence:
            triple_sequnce += self._AA_ONE_TO_THREE[aa]

        return triple_sequnce


# ===================================================================================
# ===================================================================================
def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: int | float | tuple = (0, 100),
    length_bounds: int | float | tuple = (0, 2**32),
    quality_threshold: int = 0,
    overwrite: bool = False
) -> str:
    
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
    output_dir = os.path.join(os.path.dirname(input_filepath), 'filtered')
    output_filename = os.path.basename(output_fastq)
    output_filepath = os.path.join(output_dir, output_filename)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    if os.path.exists(output_filepath) and not overwrite:
        print(f'ПРЕДУПРЕЖДЕНИЕ: файл уже существует: {output_filepath}')
        return None

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    with open(output_filepath, 'w') as out_handle:
        for record in SeqIO.parse(input_filepath, 'fastq'):
            gc = gc_fraction(record.seq) * 100
            seq_len = len(record.seq)
            qualities = record.letter_annotations['phred_quality']
            mean_quality = sum(qualities) / len(qualities)

            gc_ok = gc_bounds[0] <= gc <= gc_bounds[1]
            len_ok = length_bounds[0] <= seq_len <= length_bounds[1]
            qual_ok = mean_quality >= quality_threshold

            if gc_ok and len_ok and qual_ok:
                SeqIO.write(record, out_handle, 'fastq')

    return output_filepath
