"""
Módulo para leitura e validação de sequências de DNA.
"""

from Bio import SeqIO
from typing import Dict

def read_fasta(file_path: str) -> Dict[str, str]:
    """
    Lê um arquivo FASTA usando Biopython.
    """
    sequences = {}
    try:
        for record in SeqIO.parse(file_path, "fasta"):
            sequences[record.id] = str(record.seq)
    except Exception as e:
        print(f"Erro ao ler arquivo FASTA: {e}")
    return sequences

def validate_sequence(sequence: str) -> bool:
    """
    Verifica se a sequência contém apenas caracteres válidos (ACGTN).
    """
    valid_bases = set("ACGTN")
    return all(base.upper() in valid_bases for base in sequence)
