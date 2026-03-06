"""
GeneSeeker - Identificador de ORFs
"""

import argparse
import sys
from typing import List, Dict, Any
from geneseeker import sequence, reporter, utils
from geneseeker.domain import orf_finder, analysis

def parse_args() -> argparse.Namespace:
    """
    Processa os argumentos da linha de comando.
    """
    parser = argparse.ArgumentParser(description="GeneSeeker - Identificador de ORFs")
    parser.add_argument("input", help="Caminho para o arquivo FASTA de entrada")
    parser.add_argument("--output", help="Caminho para o arquivo de saída")
    parser.add_argument("--min-len", type=int, default=100, help="Tamanho mínimo de ORF em pb (padrão: 100)")
    parser.add_argument("--format", choices=["text", "csv", "json"], default="text", help="Formato do relatório")
    return parser.parse_args()

def load_sequences(input_path: str) -> Dict[str, str]:
    """
    Lê e valida sequências FASTA.
    """
    raw_sequences = sequence.read_fasta(input_path)
    valid_sequences = {}
    for seq_id, seq in raw_sequences.items():
        if sequence.validate_sequence(seq):
            valid_sequences[seq_id] = seq
        else:
            print(f"Aviso: Sequência '{seq_id}' contém caracteres inválidos e será ignorada.")
    return valid_sequences

def run_analysis(sequences: Dict[str, str], min_len: int) -> List[Dict[str, Any]]:
    """
    Identifica ORFs e realiza análises avançadas.
    """
    all_orfs = []
    for seq_id, dna in sequences.items():
        orfs = orf_finder.find_orfs(dna, min_length=min_len)
        for orf in orfs:
            orf["sequence_id"] = seq_id
            # Tradução e análises avançadas
            orf["translation"] = orf_finder.translate_orf(orf["sequence"])
            orf["domains"] = analysis.identify_protein_domains(orf["translation"])
            orf["splice_sites"] = analysis.predict_splice_sites(orf["sequence"])
            # Promotores (fita direta)
            orf["promoter"] = analysis.analyze_promoters(dna, orf["start"]) if orf["strand"] == "+" else None
            all_orfs.append(orf)
    return all_orfs

def output_results(orfs: List[Dict[str, Any]], format_type: str, output_path: str | None) -> None:
    """
    Exibe e salva os resultados.
    """
    reporter.display_summary(orfs)
    formatted = reporter.format_output(orfs, format_type)
    if output_path:
        reporter.save_results(formatted, output_path)
    else:
        if format_type != "text":
            print(formatted)

def main() -> None:
    """
    Ponto de entrada principal.
    """
    args = parse_args()
    try:
        sequences = load_sequences(args.input)
        if not sequences:
            print("Nenhuma sequência válida encontrada para análise.")
            return
            
        orfs = run_analysis(sequences, args.min_len)
        output_results(orfs, args.format, args.output)
        
    except Exception as e:
        print(f"Erro fatal: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
