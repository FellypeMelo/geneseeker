"""
Módulo de funções auxiliares e utilidades gerais do projeto.
"""

def clean_sequence(sequence: str) -> str:
    """
    Limpa uma sequência de DNA, removendo espaços, quebras de linha e números,
    além de padronizar para letras maiúsculas.
    
    Args:
        sequence (str): Sequência de DNA original.
        
    Returns:
        str: Sequência tratada e padronizada.
        
    TODO:
    - Remover todas as ocorrências de espaços nulos ('\n', '\t', ' ').
    - Garantir que todas as bases fiquem em maiúsculo (upper).
    """
    pass

def get_reverse_complement(sequence: str) -> str:
    """
    Gera a sequência complementar reversa de uma fita de DNA.
    
    Args:
        sequence (str): Sequência original de DNA (5' -> 3').
        
    Returns:
        str: Sequência complementar reversa (5' -> 3' da outra fita).
        
    TODO:
    - Criar o dicionário de complementariedade (A<->T, C<->G).
    - Substituir as bases da fita original pelos seus respectivos complementares.
    - Inverter a string final gerada.
    """
    pass
