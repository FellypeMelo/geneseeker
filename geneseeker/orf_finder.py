"""
Módulo contendo a lógica central de busca de Open Reading Frames (ORFs).
"""

def find_orfs(sequence: str, min_length: int = 100) -> list:
    """
    Encontra todas as ORFs em uma sequência de DNA, considerando os três frames de leitura
    na fita direta e os três na fita reversa.
    
    Args:
        sequence (str): A sequência de DNA onde as ORFs serão buscadas.
        min_length (int): Tamanho mínimo (em pares de bases) para que uma ORF seja considerada válida.
        
    Returns:
        list: Lista de dicionários (ou objetos) contendo as informações das ORFs encontradas, 
              como posição de início, fim, frame, fita e a própria sequência da ORF.
              
    TODO:
    - Identificar códons de início (ex: ATG) e de parada (ex: TAA, TAG, TGA).
    - Iterar pelos 3 frames de leitura da fita direta (forward).
    - Obter o complemento reverso para iterar pelos 3 frames de leitura da fita reversa (reverse).
    - Calcular o tamanho de cada ORF e filtrar aquelas menores que `min_length`.
    """
    pass

def translate_orf(orf_sequence: str) -> str:
    """
    Traduz a sequência de DNA de uma ORF para uma sequência de aminoácidos (proteína).
    
    Args:
        orf_sequence (str): Sequência de DNA da ORF.
        
    Returns:
        str: Sequência de aminoácidos resultante da tradução.
        
    TODO:
    - Mapear trincas de nucleotídeos (códons) para seus respectivos aminoácidos (tabela genética).
    - Iterar a sequência de 3 em 3 pares de bases.
    - Lidar adequadamente com os códons de parada.
    """
    pass
