![Описание картинки](biogen_logo.png)

# BioGen

BioGen is a toolkit for genetic data processing. It allows you to evaluate, validate, and prepare your nucleotide sequences for downstream analysis, including FASTQ filtering by GC%, length, and quality.

Authors:
* Andrey Nekrasov
* Gregor Mendel
* James D. Watson
* Francis H.C. Crick

*BioGen* documentation is available at https://site.lol/this/link/does/not/exists

## Installation

You don't need that...

## Usage and examples

### main.py

#### Biological sequence classes

`BiologicalSequence` — abstract base class for all biological sequences. Supports `len()`, indexing, slicing, and alphabet validation.

`NucleicAcidSequence(BiologicalSequence)` — base class for nucleic acids. Provides `reverse()`, `complement()`, `reverse_complement()`.
```python
seq = NucleicAcidSequence('ATGC')
seq.reverse()            # 'CGTA'
seq.complement()         # 'TACG'
seq.reverse_complement() # 'GCAT'
seq.is_valid()           # True
```

`DNASequence(NucleicAcidSequence)` — DNA sequence. Adds `transcribe()`.
```python
dna = DNASequence('ATGCGT')
dna.transcribe()         # 'AUGCGU'
dna.is_valid()           # True
DNASequence('AUGC').is_valid()  # False — U not allowed in DNA
```

`RNASequence(NucleicAcidSequence)` — RNA sequence.
```python
rna = RNASequence('AUGCGU')
rna.complement()         # 'UACGCA'
rna.is_valid()           # True
```

`AminoAcidSequence(BiologicalSequence)` — protein sequence. Adds `triple_alphabet()` for one-to-three letter conversion.
```python
prot = AminoAcidSequence('MKTLL')
prot.triple_alphabet()   # 'MetLysThrLeuLeu'
prot.is_valid()          # True
```

#### filter_fastq

Filters FASTQ reads by GC%, length, and mean Phred quality using Biopython (`SeqIO`, `gc_fraction`). Writes results directly to file without accumulating reads in memory.
```python
filter_fastq(
    input_fastq='example_data/example_fastq.fastq',
    output_fastq='output.fastq',
    gc_bounds=(20, 80),
    length_bounds=(10, 100),
    quality_threshold=20
)
```

#### Command-line usage

Run the script from the `biogen_tools` directory:
```bash
python3 main.py example_data/example_fastq.fastq filtered.fastq --gc-min 10 --gc-max 90 \
    --length-min 50 --length-max 150 --quality-threshold 20 --overwrite
```

This creates the output file in `filtered/filtered.fastq` and appends log messages to `biogen.log` by default.

You can also set a custom log file:
```bash
python3 main.py example_data/example_fastq.fastq filtered.fastq --log-file mylog.log
```

The script logs informational messages when filtering starts and finishes, and logs errors if the output file already exists or the input file is missing.

### Testing

Tests are stored in `tests/test_main.py` and cover:
- command-line parsing
- FASTQ output creation
- GC, length, and quality filtering
- overwrite protection
- error handling for missing input files
- logging of informational and error messages

Run the suite with:
```bash
python3 -m pytest -q
```


## License

The MIT License (Made In Tears)


### Disclaimer: The author takes no responsibility for any data loss, mental breakdowns, or existential crises caused by using this toolkit. You chose this path.