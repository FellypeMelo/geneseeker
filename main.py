"""
GeneSeeker - Identificador de ORFs

Propósito: Identificar Open Reading Frames (ORFs) em sequências de DNA
analisando os 3 quadros de leitura possíveis.

Um ORF é uma sequência de DNA que começa com um códon de início (ATG)
e termina com um códon de parada (TAA, TAG ou TGA).
"""

import argparse
from Bio.Seq import Seq
from Bio import SeqIO



