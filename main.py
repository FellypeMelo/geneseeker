"""
GeneSeeker - Identificador de ORFs

Propósito: Identificar Open Reading Frames (ORFs) em sequências de DNA
analisando os 3 quadros de leitura possíveis.

Um ORF é uma sequência de DNA que começa com um códon de início (ATG)
e termina com um códon de parada (TAA, TAG ou TGA).
"""

from Bio.Seq import Seq


# Constantes
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]


def find_orfs_in_frame(sequence, frame):
    """
    Busca ORFs em um quadro de leitura específico.

    Args:
        sequence: Sequência de DNA (string)
        frame: Quadro de leitura (0, 1 ou 2)

    Returns:
        list: Lista de tuplas (posição_inicio, posição_fim, sequência_orf)
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
                    orfs.append((start_pos, j + 3, orf_seq))
                    i = j  # Continua busca após este ORF
                    break
                j += 3

        i += 3

    return orfs


def analyze_all_frames(sequence):
    """
    Analisa todos os 3 quadros de leitura.

    Args:
        sequence: Sequência de DNA a ser analisada

    Returns:
        dict: Dicionário com resultados por quadro
    """
    results = {}

    print("=" * 60)
    print("GeneSeeker - Identificador de ORFs")
    print("=" * 60)
    print(f"\nSequência analisada: {sequence}")
    print(f"Comprimento: {len(sequence)} bp\n")

    for frame in range(3):
        orfs = find_orfs_in_frame(sequence, frame)
        results[frame] = orfs

        print(f"\nQuadro de Leitura {frame}:")
        print("-" * 40)

        if orfs:
            for start, end, orf_seq in orfs:
                print(f"  ORF encontrado:")
                print(f"    Posição: {start} - {end}")
                print(f"    Comprimento: {len(orf_seq)} bp")
                print(
                    f"    Sequência: {orf_seq[:50]}{'...' if len(orf_seq) > 50 else ''}"
                )
        else:
            print("  Nenhum ORF completo encontrado")

    return results


def generate_report(results, output_file="orf_report.txt"):
    """
    Gera relatório estruturado em arquivo texto.

    Args:
        results: Dicionário com resultados dos ORFs
        output_file: Nome do arquivo de saída
    """
    with open(output_file, "w") as f:
        f.write("RELATÓRIO GENESEEKER - IDENTIFICAÇÃO DE ORFs\n")
        f.write("=" * 60 + "\n\n")

        total_orfs = sum(len(orfs) for orfs in results.values())
        f.write(f"Total de ORFs encontrados: {total_orfs}\n\n")

        for frame, orfs in results.items():
            f.write(f"\nQuadro de Leitura {frame}:\n")
            f.write("-" * 40 + "\n")

            if orfs:
                for start, end, seq in orfs:
                    f.write(f"  Posição {start}-{end} ({len(seq)} bp)\n")
            else:
                f.write("  Nenhum ORF encontrado\n")

    print(f"\nRelatório salvo em: {output_file}")


def main():
    """Função principal do programa."""
    # Sequência de exemplo com ORFs
    test_sequence = "ATGCGATACTGAATGCCCTAGATGAAATAA"

    # Analisa todos os quadros
    results = analyze_all_frames(test_sequence)

    # Gera relatório
    generate_report(results)

    print("\nAnálise concluída!")


if __name__ == "__main__":
    main()
