"""
Módulo de análise biológica avançada para o GeneSeeker.
Contém lógica para promotores, splicing e domínios proteicos.
"""

import re
from typing import List, Dict, Any

def identify_protein_domains(protein_seq: str) -> List[str]:
    """
    Identifica domínios proteicos conhecidos na sequência traduzida.
    
    Args:
        protein_seq: Sequência de aminoácidos.
        
    Returns:
        List[str]: Nomes dos domínios encontrados.
    """
    found_domains = []
    
    # Padrão Prosite simplificado para Zinc Finger C2H2
    # C-x(2,4)-C-x(12)-H-x(3,5)-H
    ZINC_FINGER_REGEX = r"C.{2,4}C.{12}H.{3,5}H"
    
    if re.search(ZINC_FINGER_REGEX, protein_seq):
        found_domains.append("Zinc Finger")
        
    return found_domains

def analyze_promoters(sequence: str, start: int, upstream_len: int = 50) -> Dict[str, Any]:
    """
    Analisa a região upstream de uma ORF em busca de motivos de promotores.
    
    Args:
        sequence: Sequência de DNA original (fita correspondente).
        start: Posição de início da ORF.
        upstream_len: Tamanho da região a analisar.
        
    Returns:
        Dict: Informações sobre a região e motivos encontrados.
    """
    # Motivos biológicos (Pribnow box, -35 box, TATA box)
    KNOWN_MOTIFS = ["TATAAT", "TTGACA", "TATAAA"]
    
    upstream_start = max(0, start - upstream_len)
    upstream_seq = sequence[upstream_start:start].upper()
    
    found_motifs = [motif for motif in KNOWN_MOTIFS if motif in upstream_seq]
    
    return {
        "upstream_seq": upstream_seq,
        "found_motifs": found_motifs,
        "start_pos": upstream_start
    }

def predict_splice_sites(dna_seq: str) -> Dict[str, List[int]]:
    """
    Prediz potenciais sítios de splicing (GT-AG) dentro de um DNA.
    
    Args:
        dna_seq: Sequência de DNA da ORF.
        
    Returns:
        Dict: Listas de posições de doadores e aceitadores.
    """
    donor_sites = []
    acceptor_sites = []
    
    for i in range(len(dna_seq) - 1):
        dinucleotide = dna_seq[i : i + 2].upper()
        if dinucleotide == "GT":
            donor_sites.append(i)
        elif dinucleotide == "AG":
            acceptor_sites.append(i)
            
    return {
        "donor_sites": donor_sites,
        "acceptor_sites": acceptor_sites
    }
