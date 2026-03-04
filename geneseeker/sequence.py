"""
Módulo para leitura e validação de sequências de DNA.
"""

def read_fasta(file_path: str) -> dict:
    """
    Lê um arquivo FASTA e retorna um dicionário com os cabeçalhos como chaves 
    e as sequências como valores.
    
    Args:
        file_path (str): Caminho para o arquivo FASTA.
        
    Returns:
        dict: Dicionário onde a chave é o ID da sequência e o valor é a sequência de DNA.
        
    TODO:
    - Implementar a abertura e leitura do arquivo iterando linha por linha.
    - Tratar linhas de cabeçalho (que começam com '>').
    - Concatenar linhas de sequência associadas a cada cabeçalho.
    """
    pass

def validate_sequence(sequence: str) -> bool:
    """
    Verifica se a sequência fornecida contém apenas caracteres válidos de DNA.
    
    Args:
        sequence (str): A sequência de DNA a ser validada.
        
    Returns:
        bool: True se a sequência for válida, False caso contrário.
        
    TODO:
    - Implementar a checagem de caracteres válidos (A, C, G, T, e possivelmente N).
    - Retornar True apenas se todos os caracteres da sequência estiverem no conjunto permitido.
    """
    pass
