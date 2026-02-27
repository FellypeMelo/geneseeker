import os
import pytest
from main import find_orfs_in_frame, analyze_all_frames, START_CODON, STOP_CODONS, read_fasta_file

# --- Fixtures ---

@pytest.fixture
def sample_fasta_file(tmp_path):
    """Cria um arquivo FASTA temporário para testes."""
    content = ">seq1\nATGCGATACTGAATGCCCTAGATGAAATAA"
    file_path = tmp_path / "test.fasta"
    file_path.write_text(content)
    return str(file_path)

# --- Tests ---

def test_read_fasta_file(sample_fasta_file):
    """Testa se a função de leitura de FASTA retorna a sequência correta."""
    sequence = read_fasta_file(sample_fasta_file)
    assert sequence == "ATGCGATACTGAATGCCCTAGATGAAATAA"

def test_analyze_all_frames_six_frames():
    """Verifica se analyze_all_frames retorna resultados para 6 quadros de leitura (0-5)."""
    # A sequence with a forward ORF and a reverse ORF
    # Forward: ATG CCC TAA -> Frame 0
    # Reverse compl: TTA GGG CAT => ATG CCC TAA (Frame 0 reversed)
    sequence = "ATGCCCTAATTAGGGCAT"
    results = analyze_all_frames(sequence)
    # The dictionary should have keys 0, 1, 2, 3, 4, 5
    assert set(results.keys()) == {0, 1, 2, 3, 4, 5}
    
    # Check if protein sequence is present in the output
    assert len(results[0]) > 0
    # Expected structure: (start, end, orf_seq, protein_seq)
    for orf in results[0]:
        assert len(orf) == 4, "A tradução proteica não está sendo retornada."
        assert orf[3] == "MP*", "Tradução incorreta"
