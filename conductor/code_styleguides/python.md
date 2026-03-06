# Guia de Estilo Python (GeneSeeker / AI-XP)

Este guia estende o Google Python Style Guide com restrições específicas para engenharia científica agêntica.

## 1. Tipagem e Segurança
- **Type Hints**: Obrigatórios para todas as APIs públicas e internas (`mypy --strict`).
- **Value Objects**: Sequências biológicas devem ser instâncias de `DnaSequence`, `ProteinSequence`, etc., nunca `str` brutas.
- **Enums**: Usar para bases nitrogenadas, aminoácidos e quadros de leitura.

## 2. Estrutura e Complexidade
- **SRP (Single Responsibility)**: Cada módulo/função deve ter estritamente uma única responsabilidade.
- **Complexidade Ciclomática**: Limite máximo de **15** por função (monitorado via `radon`).
- **Nesting Depth**: Máximo de **2** níveis de aninhamento.
- **Tamanho de Funções**: Máximo de **15 linhas lógicas** por função.

## 3. Documentação (Docstrings)
- **Formato**: Google Style Docstrings obrigatórios.
- **Conteúdo**: Documentar *por que* decisões algorítmicas foram tomadas, citando referências biológicas se necessário.
- **Exemplo**:
  ```python
  def find_orfs(sequence: DnaSequence, min_len: int) -> List[Orf]:
      """Procura ORFs em 6 frames com complexidade O(n).
      
      Args:
          sequence: Objeto DnaSequence validado.
          min_len: Comprimento mínimo em nucleotídeos.
          
      Returns:
          Lista de objetos Orf encontrados.
      """
  ```

## 4. Testes e Invariantes
- **Pytest**: Framework principal.
- **Hypothesis (PBT)**: Usar para validar invariantes (ex: `len(translation(seq)) == len(seq) // 3`).
- **Mutation Testing**: Obrigatório via `mutmut` para módulos de domínio.

## 5. Clean Architecture (Mandatos)
- **Domínio Puro**: Módulos em `geneseeker.domain` não podem importar `infrastructure` ou `interface`.
- **I/O Isolado**: Proibido usar `open()`, `print()` ou `logging` dentro de funções de domínio puro.
