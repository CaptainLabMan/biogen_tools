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
run_dna_rna_tools - Performs validation and various operations on a nucleotide sequence (DNA or RNA).
```python
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False !!
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```

filter_fastq - Filter FASTQ reads by GC%, length, and mean Phred quality; write passing reads to output.
```python
filter_fastq(seqs = EXAMPLE_FASTQ, gc_bounds = (20, 80), length_bounds = (10, 30), quality_threshold = 10)
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