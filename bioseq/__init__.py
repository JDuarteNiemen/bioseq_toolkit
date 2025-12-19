"""
bioseq: lightweight utilities for genomic and proteomic sequence analysis
"""

from .fasta import readFASTAseq, writeFASTA
from .orf import candidateProtein, maximalORF
from .translation import translate
from .alignment import simple_align, seeded_simple_align
from .scoring import Scoring