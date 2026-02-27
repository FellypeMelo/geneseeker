import pytest
from Bio.Seq import Seq
from main import find_orfs_in_frame, analyze_all_frames

def test_find_orfs_with_min_len():
    """Testa a filtragem de ORFs por tamanho mínimo."""
    # ORF 1: ATG CGA TAC TGA (12 bp)
    # ORF 2: ATG CCC GGG AAA CCC TTT TAA (18 bp)
    sequence = "ATGCGATACTGAATGCCCGGGAAACCCTTTTAA"
    
    # Teste sem filtro (padrão deve ser mantido ou definido)
    results_all = analyze_all_frames(sequence)
    total_orfs_all = sum(len(orfs) for orfs in results_all.values())
    assert total_orfs_all >= 2
    
    # Teste com filtro de 15 bp (deve encontrar apenas o de 18bp)
    results_filtered = analyze_all_frames(sequence, min_length=15)
    total_orfs_filtered = sum(len(orfs) for orfs in results_filtered.values())
    assert total_orfs_filtered == 1
    
    # Verifica se o ORF encontrado é o correto
    found = False
    for frame_orfs in results_filtered.values():
        for start, end, seq, prot in frame_orfs:
            if len(seq) == 21:
                found = True
    assert found

def test_find_orfs_with_very_high_min_len():
    """Testa filtro que remove todos os ORFs."""
    sequence = "ATGCGATACTGAATGCCCGGGAAATAA"
    results = analyze_all_frames(sequence, min_length=100)
    total_orfs = sum(len(orfs) for orfs in results.values())
    assert total_orfs == 0
