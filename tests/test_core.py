import pytest
from geneseeker import sequence
from geneseeker.domain import orf_finder

def test_read_fasta_empty(tmp_path):
    """Testa leitura de arquivo vazio."""
    d = tmp_path / "sub"
    d.mkdir()
    p = d / "empty.fasta"
    p.write_text("")
    
    seqs = sequence.read_fasta(str(p))
    assert seqs == {}

def test_validate_sequence():
    """Testa validação de DNA."""
    assert sequence.validate_sequence("ATGC") == True
    assert sequence.validate_sequence("atgc") == True
    assert sequence.validate_sequence("ATGCN") == True
    assert sequence.validate_sequence("ATGCX") == False

def test_orf_finder_basic():
    """Testa busca básica de ORF."""
    seq = "ATG" + "CCC" * 5 + "TAA" # 21 bp
    orfs = orf_finder.find_orfs(seq, min_length=15)
    assert len(orfs) == 1
    assert orfs[0]["start"] == 0
    assert orfs[0]["end"] == 21
    assert orfs[0]["strand"] == "+"

def test_orf_finder_reverse():
    """Testa busca de ORF na fita reversa."""
    # ORF na fita reversa: ATG CCC TAA
    # Na fita direta (forward): TTA GGG CAT
    seq = "TTAGGGCAT"
    orfs = orf_finder.find_orfs(seq, min_length=9)
    assert len(orfs) == 1
    assert orfs[0]["strand"] == "-"
    assert orfs[0]["start"] == 0
    assert orfs[0]["end"] == 9
