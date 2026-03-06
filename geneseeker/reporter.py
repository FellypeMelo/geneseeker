"""
Módulo para exibição e salvamento de resultados da análise de ORFs.
"""

import json
import csv
from typing import List, Dict, Any

def display_summary(orfs: List[Dict[str, Any]]) -> None:
    """
    Exibe um resumo tabular dos resultados no terminal.
    """
    print("-" * 60)
    print(f"{'ID':<10} {'Início':<8} {'Fim':<8} {'Fita':<5} {'Tam (bp)':<8}")
    print("-" * 60)
    for i, orf in enumerate(orfs[:10]):
        print(f"{i+1:<10} {orf['start']:<8} {orf['end']:<8} {orf['strand']:<5} {orf['length']:<8}")
    
    if len(orfs) > 10:
        print(f"... e mais {len(orfs) - 10} ORFs.")
    print("-" * 60)
    print(f"Total de ORFs identificadas: {len(orfs)}")
    print("-" * 60)

def format_output(orfs: List[Dict[str, Any]], format_type: str) -> str:
    """
    Formata a lista de ORFs para a string final no formato escolhido.
    """
    if format_type == "json":
        return json.dumps(orfs, indent=2)
    
    if format_type == "csv":
        if not orfs:
            return ""
        header = list(orfs[0].keys())
        import io
        output = io.StringIO()
        writer = csv.DictWriter(output, fieldnames=header)
        writer.writeheader()
        writer.writerows(orfs)
        return output.getvalue()
    
    # Texto (Padrão)
    lines = ["RELATÓRIO GENESEEKER - IDENTIFICAÇÃO DE ORFs", "=" * 60, ""]
    for i, orf in enumerate(orfs):
        lines.append(f"ORF #{i+1}")
        lines.append(f"  Posição: {orf['start']} - {orf['end']} (Fita {orf['strand']}, Frame {orf['frame']})")
        lines.append(f"  Comprimento: {orf['length']} bp")
        lines.append(f"  DNA: {orf['sequence'][:50]}...")
        lines.append("-" * 20)
    return "\n".join(lines)

def save_results(formatted_output: str, output_path: str) -> None:
    """
    Grava a string formatada no caminho especificado.
    """
    try:
        with open(output_path, "w", encoding="utf-8") as f:
            f.write(formatted_output)
        print(f"Resultados salvos em: {output_path}")
    except Exception as e:
        print(f"Erro ao salvar arquivo: {e}")
