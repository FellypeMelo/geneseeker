import pytest
from Bio.Seq import Seq
from geneseeker.domain import orf_finder, analysis

def test_identify_protein_domains():
    """Testa a identificação de domínios proteicos."""
    # Sequência com motivo de dedo de zinco (C-x(2,4)-C-x(12)-H-x(3,5)-H)
    protein_seq = "ACCCCDDDDDDDDDDDDHAAAH"
    domains = analysis.identify_protein_domains(protein_seq)
    assert domains
    assert "Zinc Finger" in domains

def test_identify_protein_domains_none():
    """Testa quando nenhum domínio é encontrado."""
    protein_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    domains = analysis.identify_protein_domains(protein_seq)
    assert not domains

def test_predict_splice_sites():
    """Testa a predição de sítios de splicing (GT-AG)."""
    # GT na pos 6, AG na pos 12
    sequence = "ATG" + "CCC" + "GTA" + "AAA" + "AGG" + "TAA"
    splice_info = analysis.predict_splice_sites(sequence)
    
    assert 6 in splice_info["donor_sites"] # GT na pos 6
    assert 12 in splice_info["acceptor_sites"] # AG na pos 12

def test_analyze_promoters():
    """Testa a análise de motivos em regiões upstream."""
    # TATA box (TATAAT) na região upstream
    sequence = "CCCGGGTATAATGGGCCCATGAAATAA"
    # ATG começa na pos 18
    promoter_info = analysis.analyze_promoters(sequence, 18)
    assert "TATAAT" in promoter_info["found_motifs"]

def test_find_orfs_min_length():
    """Testa o filtro de tamanho mínimo."""
    sequence = "ATG" + "CCC" * 30 + "TAA" # 96 bp
    # Busca com min_length 100
    orfs = orf_finder.find_orfs(sequence, min_length=100)
    assert len(orfs) == 0
    
    # Busca com min_length 90
    orfs = orf_finder.find_orfs(sequence, min_length=90)
    assert len(orfs) == 1

def test_analyze_all_6_frames():
    """Testa se encontra ORFs em ambas as fitas."""
    # Forward ORF
    forward_seq = "ATG" + "CCC" * 10 + "TAA"
    # Reverse ORF (será o complemento reverso da fita direta)
    rev_seq_obj = Seq("ATG" + "GGG" * 10 + "TAG")
    reverse_seq_in_forward = str(rev_seq_obj.reverse_complement())
    
    # Combina sequências separadas por Ns para não criar ORFs artificiais
    combined = forward_seq + "NNNNNN" + reverse_seq_in_forward
    
    orfs = orf_finder.find_orfs(combined, min_length=30)
    
    strands = [o["strand"] for o in orfs]
    assert "+" in strands
    assert "-" in strands
