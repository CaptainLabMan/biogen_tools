![Описание картинки](biogen_logo.png)

# BioGen

BioGen is a toolkit for genetic data processing. It allows you to evaluate, validate, and prepare your nucleotide sequences for downstream analysis, including FASTQ filtering by GC%, length, and quality; FASTA reformatting; BLAST best-hit extraction; and retrieval of flanking genes from GenBank files.

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

### bio_files_processor.py
convert_multiline_fasta_to_oneline – Converts multi-line FASTA sequences into a single-line format and saves the result.
```python
convert_multiline_fasta_to_oneline(input_fasta: "input.fasta", output_fasta: "output.fasta")
```

parse_blast_output - Parses a BLAST output file and keeps only the best match for each query sequence.
```python
parse_blast_output(input_file: "input.txt", output_file: "output.txt")
```

select_genes_from_gbk_to_fasta - Extracts the nearest flanking genes of a target gene from a GenBank file and saves them in FASTA format.
```python
select_genes_from_gbk_to_fasta(input_gbk: "input.gbk", genes: ["gene1", "gene2"], n_before: 20, n_after: 25, output_fasta: "output.fasta")
```


## License

The MIT License (Made In Tears)


### Disclaimer: The author takes no responsibility for any data loss, mental breakdowns, or existential crises caused by using this toolkit. You chose this path.