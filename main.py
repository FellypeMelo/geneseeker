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

# Constantes
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]


def read_fasta_file(file_path):
    """
    Lê uma sequência de um arquivo FASTA.

    Args:
        file_path: Caminho para o arquivo FASTA.

    Returns:
        str: Sequência de DNA.
    """
    try:
        record = SeqIO.read(file_path, "fasta")
        return str(record.seq)
    except Exception as e:
        print(f"Erro ao ler arquivo FASTA: {e}")
        return ""


import re

def identify_protein_domains(orf_info):
    """
    Identifica domínios proteicos conhecidos na sequência traduzida de um ORF.

    Args:
        orf_info: Tupla (start, end, seq, prot).

    Returns:
        list: Nomes dos domínios encontrados.
    """
    start, end, orf_seq, prot = orf_info
    found_domains = []
    
    # Padrões Prosite (Simplificados para o exemplo)
    # Zinc Finger C2H2: C-x(2,4)-C-x(12)-H-x(3,5)-H
    ZINC_FINGER_REGEX = r"C.{2,4}C.{12}H.{3,5}H"
    
    if re.search(ZINC_FINGER_REGEX, prot):
        found_domains.append("Zinc Finger")
        
    return found_domains


def analyze_promoters(sequence, orf_info, upstream_len=50):
    """
    Analisa a região upstream de um ORF em busca de motivos de promotores.

    Args:
        sequence: Sequência de DNA original.
        orf_info: Tupla (start, end, seq, prot).
        upstream_len: Tamanho da região upstream a ser analisada.

    Returns:
        dict: Informações sobre a região upstream e motivos encontrados.
    """
    start, end, orf_seq, prot = orf_info
    
    # Motivos conhecidos (Simplificado)
    # TATAAT: Pribnow box (Procariotos)
    # TTGACA: -35 box (Procariotos)
    # TATAAA: TATA box (Eucariotos)
    KNOWN_MOTIFS = ["TATAAT", "TTGACA", "TATAAA"]
    
    # Extrai região upstream
    upstream_start = max(0, start - upstream_len)
    upstream_seq = sequence[upstream_start:start].upper()
    
    found_motifs = []
    for motif in KNOWN_MOTIFS:
        if motif in upstream_seq:
            found_motifs.append(motif)
            
    return {
        "upstream_seq": upstream_seq,
        "found_motifs": found_motifs,
        "start_pos": upstream_start
    }


def predict_splice_sites(orf_info):
    """
    Prediz potenciais sítios de splicing (GT-AG) dentro de um ORF.

    Args:
        orf_info: Tupla (start, end, seq, prot).

    Returns:
        dict: Listas de posições de doadores e aceitadores.
    """
    start, end, orf_seq, prot = orf_info
    
    donor_sites = []
    acceptor_sites = []
    
    # Busca por GT (Doador) e AG (Aceitador)
    for i in range(len(orf_seq) - 1):
        dinucleotide = orf_seq[i : i + 2]
        if dinucleotide == "GT":
            donor_sites.append(i)
        elif dinucleotide == "AG":
            acceptor_sites.append(i)
            
    return {
        "donor_sites": donor_sites,
        "acceptor_sites": acceptor_sites
    }


def find_orfs_in_frame(sequence, frame, min_length=0):
    """
    Busca ORFs em um quadro de leitura específico.

    Args:
        sequence: Sequência de DNA (string)
        frame: Quadro de leitura (0, 1 ou 2)
        min_length: Tamanho mínimo da ORF em pb

    Returns:
        list: Lista de tuplas (posição_inicio, posição_fim, sequência_orf, sequência_proteina)
    """
    orfs = []
    seq = str(sequence).upper()
    i = frame

    while i < len(seq) - 2:
        codon = seq[i : i + 3]

        # Encontrou códon de início
        if codon == START_CODON:
            start_pos = i
            j = i + 3

            # Procura códon de parada
            while j < len(seq) - 2:
                stop_codon = seq[j : j + 3]
                if stop_codon in STOP_CODONS:
                    orf_seq = seq[start_pos : j + 3]
                    if len(orf_seq) >= min_length:
                        # Tradução usando Biopython
                        protein_seq = str(Seq(orf_seq).translate())
                        orfs.append((start_pos, j + 3, orf_seq, protein_seq))
                    i = j  # Continua busca após este ORF
                    break
                j += 3

        i += 3

    return orfs


def analyze_all_frames(sequence, min_length=0):
    """
    Analisa todos os 6 quadros de leitura (3 forward, 3 reverse complement).

    Args:
        sequence: Sequência de DNA a ser analisada
        min_length: Tamanho mínimo do ORF

    Returns:
        dict: Dicionário com resultados por quadro (0 a 5)
    """
    results = {}
    
    seq_obj = Seq(sequence)
    rev_comp = str(seq_obj.reverse_complement())
    seq_len = len(sequence)

    print("=" * 60)
    print("GeneSeeker - Identificador de ORFs")
    print("=" * 60)
    print(f"\nSequência analisada: {sequence[:30]}{'...' if len(sequence)>30 else ''} ({seq_len} bp)")

    for frame in range(6):
        if frame < 3:
            orfs = find_orfs_in_frame(sequence, frame, min_length)
            resultados_ajustados = orfs
            tipo_fita = "Forward"
        else:
            # Quadros 3, 4 e 5 são a fita complementar reversa
            orfs = find_orfs_in_frame(rev_comp, frame - 3, min_length)
            resultados_ajustados = []
            for start, end, orf_seq, protein_seq in orfs:
                # Na fita reversa, o end (j+3) de rev complementa o start original inversamente
                real_start = seq_len - end
                real_end = seq_len - start
                resultados_ajustados.append((real_start, real_end, orf_seq, protein_seq))
            tipo_fita = "Reverse"

        results[frame] = resultados_ajustados

        print(f"\nQuadro de Leitura {frame} ({tipo_fita}):")
        print("-" * 40)

        if resultados_ajustados:
            for start, end, orf_seq, protein_seq in resultados_ajustados:
                # Análises adicionais para verificação
                promoter_info = analyze_promoters(sequence if frame < 3 else rev_comp, (start, end, orf_seq, protein_seq))
                splice_info = predict_splice_sites((start, end, orf_seq, protein_seq))
                
                print(f"  ORF encontrado:")
                print(f"    Posição: {start} - {end}")
                print(f"    Comprimento: {len(orf_seq)} bp")
                print(f"    Seq. DNA: {orf_seq[:30]}{'...' if len(orf_seq) > 30 else ''}")
                print(f"    Seq. Prot: {protein_seq[:30]}{'...' if len(protein_seq) > 30 else ''}")
                
                if promoter_info["found_motifs"]:
                    print(f"    Motivos Upstream: {', '.join(promoter_info['found_motifs'])}")
                
                if splice_info["donor_sites"] or splice_info["acceptor_sites"]:
                    print(f"    Sítios Splicing: Doador(GT) em {splice_info['donor_sites']}, Aceitador(AG) em {splice_info['acceptor_sites']}")
        else:
            print("  Nenhum ORF completo encontrado")

    return results


def generate_report(results, output_file="orf_report.txt"):
    """
    Gera relatório estruturado em arquivo texto.

    Args:
        results: Dicionário com resultados dos ORFs (todos os 6 quadros)
        output_file: Nome do arquivo de saída
    """
    with open(output_file, "w", encoding="utf-8") as f:
        f.write("RELATÓRIO GENESEEKER - IDENTIFICAÇÃO DE ORFs\n")
        f.write("=" * 60 + "\n\n")

        total_orfs = sum(len(orfs) for orfs in results.values())
        f.write(f"Total de ORFs encontrados: {total_orfs}\n\n")

        for frame, orfs in results.items():
            f.write(f"\nQuadro de Leitura {frame}:\n")
            f.write("-" * 40 + "\n")

            if orfs:
                for start, end, seq, prot in orfs:
                    f.write(f"  Posição {start}-{end} ({len(seq)} bp) -> Prot: {prot[:20]}{'...' if len(prot)>20 else ''}\n")
            else:
                f.write("  Nenhum ORF encontrado\n")

    print(f"\nRelatório salvo em: {output_file}")


def main():
    """Função principal do programa."""
    parser = argparse.ArgumentParser(description="GeneSeeker - Identificador de ORFs em sequências de DNA")
    parser.add_argument("input_file", help="Caminho para o arquivo FASTA contendo a sequência de DNA.")
    parser.add_argument("-o", "--output", default="orf_report.txt", help="Nome do arquivo de saída para o relatório (padrão: orf_report.txt)")
    parser.add_argument("--min-length", type=int, default=0, help="Tamanho mínimo da ORF em pares de bases (bp)")
    
    args = parser.parse_args()

    print(f"Lendo sequência do arquivo: {args.input_file}")
    sequence = read_fasta_file(args.input_file)
    
    if not sequence:
        print("Erro: Sequência inválida ou arquivo vazio.")
        return

    # Analisa todos os quadros
    results = analyze_all_frames(sequence, min_length=args.min_length)

    # Gera relatório
    generate_report(results, args.output)

    print("\nAnálise concluída!")


if __name__ == "__main__":
    main()
