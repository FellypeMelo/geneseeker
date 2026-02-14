#!/usr/bin/env python3
"""
GeneSeeker - Gerador de Dados de Teste

Este script gera 50+ conjuntos de dados sintéticos com ORFs conhecidos.
Útil para validar se o GeneSeeker está identificando corretamente os frames.

Os dados de teste são COMMITADOS no GitHub.
Para dados reais, use a pasta data/ (gitignored)
"""

import random
import os
from datetime import datetime

TEST_DATA_DIR = "test_data"
NUM_DATASETS = 55

# Códons
START_CODON = "ATG"
STOP_CODONS = ["TAA", "TAG", "TGA"]
CODONS = [
    "ATG",
    "GCT",
    "GCC",
    "GCA",
    "GCG",
    "CGT",
    "CGC",
    "CGA",
    "CGG",
    "AGA",
    "AGG",
    "AAT",
    "AAC",
    "GAT",
    "GAC",
    "TGT",
    "TGC",
    "CAA",
    "CAG",
    "GAA",
    "GAG",
    "GGT",
    "GGC",
    "GGA",
    "GGG",
    "CAT",
    "CAC",
    "ATT",
    "ATC",
    "ATA",
    "TTA",
    "TTG",
    "CTT",
    "CTC",
    "CTA",
    "CTG",
    "AAA",
    "AAG",
    "TTT",
    "TTC",
    "CCT",
    "CCC",
    "CCA",
    "CCG",
    "TCT",
    "TCC",
    "TCA",
    "TCG",
    "AGT",
    "AGC",
    "ACT",
    "ACC",
    "ACA",
    "ACG",
    "GTT",
    "GTC",
    "GTA",
    "GTG",
    "TGG",
    "TAT",
    "TAC",
    "TAA",
    "TAG",
    "TGA",
]


def generate_random_codon():
    """Gera um códon aleatório válido."""
    return random.choice(CODONS)


def generate_sequence_with_orfs(num_orfs=3, min_orf_length=5, max_orf_length=15):
    """
    Gera uma sequência com ORFs específicos em frames diferentes.

    Returns:
        tuple: (sequência, lista de ORFs encontrados)
    """
    sequence = []
    orfs_info = []

    # Adiciona lixo inicial (frame 0)
    initial_junk_len = random.randint(0, 5) * 3
    for _ in range(initial_junk_len // 3):
        codon = generate_random_codon()
        while codon == START_CODON or codon in STOP_CODONS:
            codon = generate_random_codon()
        sequence.append(codon)

    # Adiciona ORF no frame 0
    orf0_len = random.randint(min_orf_length, max_orf_length)
    sequence.append(START_CODON)
    for _ in range(orf0_len - 2):
        codon = generate_random_codon()
        while codon in STOP_CODONS:
            codon = generate_random_codon()
        sequence.append(codon)
    sequence.append(random.choice(STOP_CODONS))
    orfs_info.append((0, len(sequence) * 3 - 3, orf0_len))

    # Adiciona ORF no frame 1 (começa na posição 1)
    # Primeiro adicionamos 1 base
    sequence.append(random.choice(["A", "T", "G", "C"]))
    start_pos_frame1 = len(sequence) - 1
    sequence.append(START_CODON)
    orf1_len = random.randint(min_orf_length, max_orf_length)
    for _ in range(orf1_len - 2):
        codon = generate_random_codon()
        while codon in STOP_CODONS:
            codon = generate_random_codon()
        sequence.append(codon)
    sequence.append(random.choice(STOP_CODONS))
    orfs_info.append((1, start_pos_frame1, orf1_len))

    # Adiciona ORF no frame 2 (começa na posição 2)
    # Adicionamos 2 bases
    sequence.append(random.choice(["A", "T", "G", "C"]))
    sequence.append(random.choice(["A", "T", "G", "C"]))
    start_pos_frame2 = len(sequence) - 2
    sequence.append(START_CODON)
    orf2_len = random.randint(min_orf_length, max_orf_length)
    for _ in range(orf2_len - 2):
        codon = generate_random_codon()
        while codon in STOP_CODONS:
            codon = generate_random_codon()
        sequence.append(codon)
    sequence.append(random.choice(STOP_CODONS))
    orfs_info.append((2, start_pos_frame2, orf2_len))

    # Junta tudo
    full_sequence = "".join(sequence)

    return full_sequence, orfs_info


def generate_no_orf_sequence():
    """Gera uma sequência sem ORFs completos."""
    sequence = []
    length = random.randint(30, 100)

    for _ in range(length):
        codon = generate_random_codon()
        # Evita START e STOP
        while codon == START_CODON or codon in STOP_CODONS:
            codon = generate_random_codon()
        sequence.append(codon)

    return "".join(sequence)


def generate_multiple_orfs_sequence(num_orfs=5):
    """Gera uma sequência com múltiplos ORFs no mesmo frame."""
    sequence = []
    orfs_info = []

    for i in range(num_orfs):
        # Adiciona lixo entre ORFs
        if i > 0:
            junk_len = random.randint(2, 8)
            for _ in range(junk_len):
                codon = generate_random_codon()
                while codon == START_CODON or codon in STOP_CODONS:
                    codon = generate_random_codon()
                sequence.append(codon)

        # Adiciona ORF
        orf_len = random.randint(5, 15)
        start_pos = len(sequence) * 3
        sequence.append(START_CODON)
        for _ in range(orf_len - 2):
            codon = generate_random_codon()
            while codon in STOP_CODONS:
                codon = generate_random_codon()
            sequence.append(codon)
        sequence.append(random.choice(STOP_CODONS))
        orfs_info.append((0, start_pos, orf_len))

    return "".join(sequence), orfs_info


def generate_overlapping_orfs():
    """Gera ORFs sobrepostos (situação complexa)."""
    # Frame 0: ATG AAA TAA
    # Frame 1:  TGA AAT AA
    # Frame 2:   GAA ATA A
    sequence = "ATGAAATAA"
    return sequence, [(0, 0, 3)]


def generate_test_datasets():
    """Gera todos os datasets de teste."""
    datasets = []

    # 1-10: Sequências simples com ORFs
    for i in range(10):
        seq, orfs = generate_sequence_with_orfs()
        datasets.append(
            (f"geneseeker_test_{i + 1:02d}_3orfs", seq, f"3 ORFs known: {orfs}")
        )

    # 11-15: Sem ORFs
    for i in range(10, 15):
        seq = generate_no_orf_sequence()
        datasets.append(
            (f"geneseeker_test_{i + 1:02d}_no_orfs", seq, "No complete ORFs")
        )

    # 16-25: Múltiplos ORFs
    for i in range(15, 25):
        num = random.randint(3, 8)
        seq, orfs = generate_multiple_orfs_sequence(num)
        datasets.append(
            (
                f"geneseeker_test_{i + 1:02d}_multi_{num}orfs",
                seq,
                f"{num} ORFs in frame 0",
            )
        )

    # 26-30: Muito curtas
    for i in range(25, 30):
        seq = generate_random_codon() * random.randint(5, 10)
        datasets.append((f"geneseeker_test_{i + 1:02d}_short", seq, "Short sequence"))

    # 31-35: Muito longas
    for i in range(30, 35):
        seq, orfs = generate_sequence_with_orfs(min_orf_length=20, max_orf_length=50)
        datasets.append((f"geneseeker_test_{i + 1:02d}_long", seq, "Long ORFs"))

    # 36-40: Apenas frame 0
    for i in range(35, 40):
        seq = (
            START_CODON
            + "".join([generate_random_codon() for _ in range(10)])
            + random.choice(STOP_CODONS)
        )
        datasets.append(
            (f"geneseeker_test_{i + 1:02d}_frame0_only", seq, "Only frame 0 has ORF")
        )

    # 41-45: ORFs sobrepostos
    for i in range(40, 45):
        seq, orfs = generate_overlapping_orfs()
        # Adiciona lixo antes e depois
        prefix = "".join([generate_random_codon() for _ in range(3)])
        suffix = "".join([generate_random_codon() for _ in range(3)])
        seq = prefix + seq + suffix
        datasets.append(
            (f"geneseeker_test_{i + 1:02d}_overlapping", seq, "Overlapping ORFs")
        )

    # 46-50: Aleatórios
    for i in range(45, 50):
        if random.random() < 0.5:
            seq, orfs = generate_sequence_with_orfs()
            desc = f"Random with {len(orfs)} ORFs"
        else:
            seq = generate_no_orf_sequence()
            desc = "Random without ORFs"
        datasets.append((f"geneseeker_test_{i + 1:02d}_random", seq, desc))

    # 51-55: Casos extremos
    # 51: ORF muito curto (mínimo: ATG + 1 codon + STOP)
    seq = START_CODON + generate_random_codon() + random.choice(STOP_CODONS)
    datasets.append(("geneseeker_test_51_minimal_orf", seq, "Minimal ORF (3 codons)"))

    # 52: Muitos START sem STOP
    seq = START_CODON * 5 + "".join([generate_random_codon() for _ in range(20)])
    datasets.append(("geneseeker_test_52_many_starts", seq, "Many STARTs, no STOP"))

    # 53: Muitos STOP sem START
    seq = "".join([random.choice(STOP_CODONS) for _ in range(5)]) + "".join(
        [generate_random_codon() for _ in range(10)]
    )
    datasets.append(("geneseeker_test_53_many_stops", seq, "Many STOPs, no START"))

    # 54: Apenas START
    seq = START_CODON * 10
    datasets.append(("geneseeker_test_54_only_starts", seq, "Only START codons"))

    # 55: Apenas STOP
    seq = "".join([random.choice(STOP_CODONS) for _ in range(10)])
    datasets.append(("geneseeker_test_55_only_stops", seq, "Only STOP codons"))

    return datasets


def format_fasta_sequence(sequence, line_length=60):
    """Formata sequência em linhas de 60 caracteres (padrão NCBI)."""
    lines = []
    for i in range(0, len(sequence), line_length):
        lines.append(sequence[i : i + line_length])
    return "\n".join(lines)


def save_datasets(datasets):
    """Salva os datasets em arquivos."""
    os.makedirs(TEST_DATA_DIR, exist_ok=True)

    manifest_path = os.path.join(TEST_DATA_DIR, "MANIFEST.txt")
    with open(manifest_path, "w") as manifest:
        manifest.write("=" * 70 + "\n")
        manifest.write("GeneSeeker - Dados de Teste Sintéticos\n")
        manifest.write("=" * 70 + "\n\n")
        manifest.write(f"Gerado em: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        manifest.write(f"Total de datasets: {len(datasets)}\n\n")
        manifest.write("ATENÇÃO: Estes são dados FABRICADOS com ORFs conhecidos.\n")
        manifest.write("Útil para validar se o detector está funcionando.\n\n")
        manifest.write("Lista de arquivos:\n")
        manifest.write("-" * 70 + "\n")

        for filename, sequence, description in datasets:
            filepath = os.path.join(TEST_DATA_DIR, f"{filename}.fasta")

            with open(filepath, "w") as f:
                f.write(f">{filename} {description}\n")
                f.write(format_fasta_sequence(sequence) + "\n")

            manifest.write(f"{filename}.fasta - {description} ({len(sequence)} bp)\n")
            print(f"[OK] Gerado: {filename}.fasta ({len(sequence)} bp)")

    print(f"\n[OK] Manifesto salvo em: {manifest_path}")
    print(f"[OK] Total: {len(datasets)} arquivos FASTA gerados")


def main():
    print("=" * 70)
    print("GeneSeeker - Gerador de Dados de Teste")
    print("=" * 70)
    print()

    datasets = generate_test_datasets()
    save_datasets(datasets)

    print()
    print("=" * 70)
    print("Geração concluída!")
    print("=" * 70)
    print(f"\nDados em: {TEST_DATA_DIR}/")
    print("Execute: python main.py")
    print("E edite o código para ler desses arquivos.")


if __name__ == "__main__":
    main()
