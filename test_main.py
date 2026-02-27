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
