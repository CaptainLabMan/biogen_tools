![Описание картинки](biogen_logo.png)

# BioGen

BioGen is a toolkit for genetic data processing. It allows you to evaluate, validate and prepare your nucleotide sequences for downstream analysis.

Authors:
* Andrey Nekrasov
* Gregor Mendel
* James D. Watson
* Francis H.C. Crick

*BioGen* documentation is available at https://site.lol/this/link/does/not/exists

## Installation

Right after publication of the article in Nature. Maybe never...

## Usage and examples

run_dna_rna_tools - performs validation and various operations on a nucleotide sequence (DNA or RNA).
```python
run_dna_rna_tools('TTUU', 'is_nucleic_acid') # False !!
run_dna_rna_tools('ATG', 'transcribe') # 'AUG'
run_dna_rna_tools('ATG', 'reverse') # 'GTA'
run_dna_rna_tools('AtG', 'complement') # 'TaC'
run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'
run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
```

filter_fastq - performs evaluation of various parameters of reads from a FASTQ file.
```python
EXAMPLE_FASTQ = {
    # 'name' : ('sequence', 'quality')
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
    'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    ...
}

filter_fastq(seqs = EXAMPLE_FASTQ, gc_bounds = (20, 80), length_bounds = (10, 30), quality_threshold = 10)
```

## License

The MIT License (Made In Tears)


### Disclaimer: The author takes no responsibility for any data loss, mental breakdowns, or existential crises caused by using this toolkit. You chose this path.