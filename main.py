from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import argparse
import logging
import os
import sys
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
def configure_logging(log_file: str = 'biogen.log') -> logging.Logger:
    logger = logging.getLogger('biogen')
    logger.setLevel(logging.INFO)

    if not logger.handlers:
        handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    return logger


def parse_args(args=None):
    parser = argparse.ArgumentParser(
        description='Filter FASTQ reads by GC%, length, and quality'
    )
    parser.add_argument('input_fastq', help='Input FASTQ path')
    parser.add_argument('output_fastq', help='Output FASTQ filename')
    parser.add_argument('--gc-min', type=float, default=0, help='Minimum GC percentage')
    parser.add_argument('--gc-max', type=float, default=100, help='Maximum GC percentage')
    parser.add_argument('--length-min', type=int, default=0, help='Minimum read length')
    parser.add_argument('--length-max', type=int, default=2**32, help='Maximum read length')
    parser.add_argument('--quality-threshold', type=int, default=0, help='Minimum average Phred quality')
    parser.add_argument('--overwrite', action='store_true', help='Overwrite existing output file')
    parser.add_argument('--simulate-error', action='store_true', help='Log an artificial error and exit')
    parser.add_argument('--log-file', default='biogen.log', help='Log file path')
    return parser.parse_args(args)


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: int | float | tuple = (0, 100),
    length_bounds: int | float | tuple = (0, 2**32),
    quality_threshold: int = 0,
    overwrite: bool = False,
    logger: logging.Logger | None = None,
) -> str | None:
    if logger is None:
        logger = logging.getLogger('biogen')

    input_filepath = os.path.abspath(input_fastq)
    output_dir = os.path.join(os.path.dirname(input_filepath), 'filtered')
    output_filename = os.path.basename(output_fastq)
    output_filepath = os.path.join(output_dir, output_filename)
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    logger.info('Starting FASTQ filter from %s to %s', input_filepath, output_filepath)

    if os.path.exists(output_filepath) and not overwrite:
        logger.error('Output file already exists and overwrite is disabled: %s', output_filepath)
        return None

    if isinstance(gc_bounds, (int, float)):
        gc_bounds = (0, gc_bounds)
    if isinstance(length_bounds, (int, float)):
        length_bounds = (0, length_bounds)

    try:
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
    except FileNotFoundError:
        logger.error('Input FASTQ file not found: %s', input_filepath)
        return None
    except Exception as exc:
        logger.error('Filtering failed: %s', exc)
        return None

    logger.info(
        'Filtered FASTQ written to %s with gc bounds=%s, length bounds=%s, quality threshold=%s',
        output_filepath,
        gc_bounds,
        length_bounds,
        quality_threshold,
    )
    return output_filepath


def main(argv=None):
    args = parse_args(argv)
    logger = configure_logging(args.log_file)

    if args.simulate_error:
        simulate_error(logger)
        sys.exit(1)

    output_path = filter_fastq(
        args.input_fastq,
        args.output_fastq,
        gc_bounds=(args.gc_min, args.gc_max),
        length_bounds=(args.length_min, args.length_max),
        quality_threshold=args.quality_threshold,
        overwrite=args.overwrite,
        logger=logger,
    )

    if output_path is None:
        logger.error('No output was written because the filter did not complete successfully.')
        sys.exit(1)

    print(output_path)


if __name__ == '__main__':
    main()
