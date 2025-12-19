# 🧬 BioSeq Toolkit

A lightweight Python toolkit for basic genomic and proteomic sequence analysis.
Designed for teaching, rapid prototyping, and small-scale analyses, with minimal dependencies.


### Features

	•	DNA → amino acid translation (6 reading frames)
	•	Open reading frame (ORF) detection
	•	FASTA file reading & writing
	•	Amino acid composition statistics
	•	Distance metrics & distance matrices
	•	Simple pairwise sequence alignment
	•	Support for Biopython substitution matrices (e.g. BLOSUM)


### Installation
```angular2html
git clone https://github.com/yourusername/bioseq-toolkit.git
cd bioseq-toolkit
pip install -r requirements.txt
```

### Quick examples
Translate DNA:
```angular2html
from bioseq.translation import translate

frames = translate("ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG")
print(frames["f1"])

```
Find longest ORF:
```angular2html
from bioseq.orf import candidateProtein

protein = candidateProtein(dna_sequence)

```
Simple alignment:
```angular2html
from bioseq.alignment import simple_align
from bioseq.scoring import Scoring

scoring = Scoring("BLOSUM62")
result = simple_align(("MKT", "MNT"), scoring)
```

Design philosophy

This project intentionally avoids complex optimisations and focuses on:

	•	Readability
	•	Clear algorithms
	•	Educational value

For large-scale or production analyses, consider libraries such as Biopython, scikit-bio, or SeqAn.


⚠️ Limitations

	•	ORF detection currently uses only the first forward reading frame
	•	Alignment algorithms are heuristic and not guaranteed optimal
	•	Designed for short to moderate-length sequences