import os
import pytest
from Bio.Seq import Seq
from main import find_orfs_in_frame, analyze_all_frames, read_fasta_file, generate_report, analyze_promoters

def test_analyze_promoters():
    """Testa a análise de motivos em regiões upstream."""
    # Sequência com TATA box (TATAAT) na região upstream
    # TATA box geralmente em -10 (Pribnow box em procariotos)
    # Upstream: ...TATAAT... [ATG]
    sequence = "CCCGGG TATAAT GGGCCC " + "ATG" + "AAATAA"
    # sequence = "CCCGGGTATAATGGGCCCATGAAATAA"
    sequence = sequence.replace(" ", "")
    
    # Encontra o ORF
    results = analyze_all_frames(sequence)
    orf_info = results[0][0] # Frame 0, primeiro ORF
    
    # Analisa promotores
    promoter_info = analyze_promoters(sequence, orf_info)
    
    assert promoter_info["found_motifs"]
    assert "TATAAT" in promoter_info["found_motifs"]
    assert promoter_info["upstream_seq"] != ""

def test_analyze_promoters_no_motif():
    """Testa análise em região sem motivos conhecidos."""
    sequence = "GGGGGGGGGGGGGGGGGGGG" + "ATG" + "AAATAA"
    results = analyze_all_frames(sequence)
    orf_info = results[2][0] # Está no Quadro 2
    
    promoter_info = analyze_promoters(sequence, orf_info)
    assert not promoter_info["found_motifs"]

def test_find_orfs_with_min_len():
    """Testa a filtragem de ORFs por tamanho mínimo."""
    # ORF 1: ATG CGA TAC TGA (12 bp)
    # ORF 2: ATG CCC GGG AAA CCC TTT TAA (21 bp)
    sequence = "ATGCGATACTGAATGCCCGGGAAACCCTTTTAA"
    
    # Teste sem filtro
    results_all = analyze_all_frames(sequence)
    total_orfs_all = sum(len(orfs) for orfs in results_all.values())
    assert total_orfs_all >= 2
    
    # Teste com filtro de 15 bp
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

def test_read_fasta_file(tmp_path):
    """Testa a leitura de arquivo FASTA."""
    fasta_content = ">test_seq\nATGCGT\n"
    fasta_file = tmp_path / "test.fasta"
    fasta_file.write_text(fasta_content)
    
    seq = read_fasta_file(str(fasta_file))
    assert seq == "ATGCGT"

def test_read_fasta_file_error():
    """Testa erro na leitura de arquivo FASTA inexistente."""
    seq = read_fasta_file("non_existent.fasta")
    assert seq == ""

def test_generate_report(tmp_path):
    """Testa a geração do relatório."""
    results = {
        0: [(0, 12, "ATGCGATACTGA", "MRY*")],
        1: []
    }
    report_file = tmp_path / "report.txt"
    generate_report(results, output_file=str(report_file))
    
    assert os.path.exists(report_file)
    content = report_file.read_text(encoding="utf-8")
    assert "RELATÓRIO GENESEEKER" in content
    assert "Total de ORFs encontrados: 1" in content
    assert "Quadro de Leitura 0" in content

def test_find_orfs_reverse_strand():
    """Testa a busca de ORFs na fita reversa (quadros 3-5)."""
    # Sequência: ATGCGATACTGA
    # RevComp: TCAGTATCGCAT (tem ATG na posição 7-10 em revcomp)
    # Start em revcomp: ATG -> pos 7 em revcomp
    # Stop em revcomp: TGA na pos 0-3? Não.
    
    # Vamos usar uma sequência conhecida
    sequence = "ATGCGATACTGATCAT" # ATGCGATACTGA TCAT
    # RevComp: ATGATCAGTATCGCAT
    # No RevComp (Frame 0): ATG ATC AGT ATC GCA T -> sem stop?
    # Vamos simplificar: ATG AAA TAA
    # RevComp: TTA TTT CAT -> sem start?
    
    # ATG AAA TAA (Forward)
    # TTATTT CAT (Reverse)
    
    # ATG CGA TAG (Forward)
    # CTA TCG CAT (Reverse) -> ATG na pos 9?
    
    sequence = "ATGCGA" + str(Seq("ATGCCCTAA").reverse_complement())
    # sequence = ATGCGA + TTAGGGCAT = ATGCGATTAGGGCAT
    # RevComp: ATGCCCTAATCGCAT
    # Frame 0 em RevComp: ATG CCC TAA (9 bp)
    
    results = analyze_all_frames(sequence)
    # Quadros 3, 4, 5 são reverse
    total_rev_orfs = sum(len(results[f]) for f in [3, 4, 5])
    assert total_rev_orfs >= 1
