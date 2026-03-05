# GeneSeeker - Pilha Tecnológica (Tech Stack)

## Linguagem e Runtime
- **Python 3.9+**: Linguagem principal com foco em tipagem estática opcional via `mypy` para garantir segurança e robustez científica.

## Bibliotecas de Bioinformática
- **Biopython (v1.81)**: Manipulação avançada de sequências biológicas, leitura/escrita de arquivos FASTA e cálculos de complementaridade reversa.

## Qualidade de Código e Engenharia
- **pytest**: Framework de testes unitários e de integração.
- **Hypothesis**: Utilizado para **Property-Based Testing (PBT)**, validando invariantes biológicos (ex: toda ORF começa com ATG).
- **mypy**: Verificação estática de tipos para garantir integridade de dados entre camadas.
- **radon**: Monitoramento de complexidade ciclomática (limite máx: 15).
- **mutmut**: **Mutation Testing** para validar a eficácia da suíte de testes contra bugs lógicos.
- **black & flake8**: Padronização estética e linting rigoroso.

## Infraestrutura e Ferramentas
- **SQLite** (opcional): Cache de anotações funcionais para grandes genomas.
- **GitHub Actions**: Pipeline de CI/CD para execução automática de testes e auditoria de qualidade.
- **import-linter**: Garantia de que as fronteiras da Clean Architecture não sejam violadas.
