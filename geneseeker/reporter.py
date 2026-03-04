"""
Módulo responsável pela formatação e exibição dos resultados encontrados.
"""

def format_output(orfs: list, format_type: str = 'text') -> str:
    """
    Formata a lista de ORFs encontradas em um formato específico para exibição ou salvamento.
    
    Args:
        orfs (list): Lista de ORFs encontradas (resultado do orf_finder).
        format_type (str): Formato desejado para a saída ('text', 'csv', 'json', etc.).
        
    Returns:
        str: String formatada contendo os resultados.
        
    TODO:
    - Implementar a lógica para formatar os resultados de forma legível em texto (tabelas, listas).
    - Adicionar suporte a outros formatos estruturados, se necessário.
    """
    pass

def save_results(formatted_output: str, output_path: str) -> None:
    """
    Salva os resultados formatados em um arquivo no disco.
    
    Args:
        formatted_output (str): A string contendo os resultados devidamente formatados.
        output_path (str): Caminho do arquivo onde os resultados serão salvos.
        
    TODO:
    - Abrir o arquivo especificado no modo de escrita ('w').
    - Escrever o conteúdo da string `formatted_output`.
    - Lidar com possíveis exceções (ex: permissão negada, diretório inexistente).
    """
    pass

def display_summary(orfs: list) -> None:
    """
    Exibe um resumo tabular ou numérico no terminal sobre as ORFs processadas.
    
    Args:
        orfs (list): Lista de ORFs encontradas.
        
    TODO:
    - Calcular métricas a partir da lista de ORFs (quantidade total, tamanhos médios, maior e menor ORF).
    - Imprimir as informações de forma limpa e clara utilizando a saída padrão.
    """
    pass
