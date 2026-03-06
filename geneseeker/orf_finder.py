"""
Módulo contendo a lógica central de busca de Open Reading Frames (ORFs).
"""

from Bio.Seq import Seq
from typing import List, Dict, Any

START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]

def find_orfs_in_frame(sequence: str, frame: int, min_length: int = 100) -> List[Dict[str, Any]]:
    """
    Busca ORFs em um quadro de leitura específico.
    
    Args:
        sequence: Sequência de DNA (string).
        frame: Quadro de leitura (0, 1 ou 2).
        min_length: Tamanho mínimo em pares de bases.
        
    Returns:
        List[Dict]: Informações das ORFs encontradas.
    """
    orfs = []
    seq = str(sequence).upper()
    i = frame
    
    while i < len(seq) - 2:
        if seq[i : i + 3] == START_CODON:
            start_pos = i
            j = i + 3
            while j < len(seq) - 2:
                if seq[j : j + 3] in STOP_CODONS:
                    orf_len = (j + 3) - start_pos
                    if orf_len >= min_length:
                        orf_seq = seq[start_pos : j + 3]
                        orfs.append({
                            "start": start_pos,
                            "end": j + 3,
                            "sequence": orf_seq,
                            "length": orf_len
                        })
                    i = j
                    break
                j += 3
        i += 3
    return orfs

def find_orfs(sequence: str, min_length: int = 100) -> List[Dict[str, Any]]:
    """
    Encontra ORFs em todos os 6 quadros de leitura.
    
    Args:
        sequence: Sequência de DNA.
        min_length: Comprimento mínimo em pares de bases.
        
    Returns:
        List[Dict]: Lista consolidada de ORFs.
    """
    all_orfs = []
    seq_obj = Seq(sequence)
    rev_comp = str(seq_obj.reverse_complement())
    seq_len = len(sequence)
    
    for frame in range(6):
        if frame < 3:
            # Fita Direta
            orfs = find_orfs_in_frame(sequence, frame, min_length)
            for orf in orfs:
                orf.update({"frame": frame, "strand": "+"})
                all_orfs.append(orf)
        else:
            # Fita Reversa
            orfs = find_orfs_in_frame(rev_comp, frame - 3, min_length)
            for orf in orfs:
                # Ajusta coordenadas para a fita original
                real_start = seq_len - orf["end"]
                real_end = seq_len - orf["start"]
                orf.update({
                    "start": real_start,
                    "end": real_end,
                    "frame": frame - 3,
                    "strand": "-"
                })
                all_orfs.append(orf)
                
    return all_orfs

def translate_orf(orf_sequence: str) -> str:
    """
    Traduz a sequência de DNA de uma ORF para aminoácidos.
    """
    return str(Seq(orf_sequence).translate(to_stop=True))
