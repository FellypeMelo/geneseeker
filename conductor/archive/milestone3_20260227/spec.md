# Specification: Milestone 3 - Filtros e Análises Avançadas

## Objetivo
Implementar funcionalidades avançadas de filtragem e análise biológica no GeneSeeker para aumentar a precisão na identificação de genes funcionais.

## Escopo

### 1. Filtragem de ORFs
- **Requisito**: Permitir que o usuário defina um tamanho mínimo (em nucleotídeos ou aminoácidos) para que um ORF seja incluído no relatório.
- **Entrada**: Parâmetro `min_len` (padrão: 100 bp).
- **Saída**: Lista de ORFs filtrada.

### 2. Análise de Região Upstream (Promotores)
- **Requisito**: Extrair e analisar a sequência imediatamente anterior ao códon de início (ex: 50-100 bp upstream).
- **Funcionalidade**: Busca por motivos de consenso (ex: TATA box, Pribnow box) para validar o potencial de expressão do gene.

### 3. Predição de Splice Sites
- **Requisito**: Identificar potenciais junções exon-intron (GT-AG rule) para suporte básico a genomas eucariotos.
- **Funcionalidade**: Marcar regiões com alta probabilidade de splicing dentro de ORFs longos.

### 4. Identificação de Domínios Proteicos
- **Requisito**: Traduzir ORFs para aminoácidos e buscar por padrões de domínios conhecidos.
- **Funcionalidade**: Integração com bibliotecas de busca de motivos (ex: Prosite patterns ou integração simplificada com HMMER).

## Considerações Técnicas
- Manter a modularidade no `main.py` ou criar um módulo `analysis.py`.
- Utilizar funções do Biopython (`Bio.Motif`, `Bio.Seq`) para otimizar a implementação.
- Garantir que os testes unitários cubram os novos filtros e preditores.
