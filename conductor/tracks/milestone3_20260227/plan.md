# Implementation Plan: Milestone 3 - Filtros e Análises Avançadas

## Phase 1: Filtragem e Infraestrutura de Análise [checkpoint: fcbdd66]
Implementar o filtro de tamanho mínimo e reorganizar o código para suportar análises adicionais.

- [x] Task: Preparar ambiente de testes para Milestone 3 [checkpoint: 4e19c96]
    - [x] Criar arquivo de testes `tests/test_analysis.py`
- [x] Task: Implementar filtro de tamanho mínimo de ORF [checkpoint: de035a3]
    - [x] Write Tests: Testar filtragem por diferentes comprimentos (bp e aa)
    - [x] Implement Feature: Adicionar parâmetro `min_len` na busca de ORFs
- [x] Task: Conductor - User Manual Verification 'Phase 1: Filtragem' (Protocol in workflow.md)

## Phase 2: Análise de Promotores e Splicing
Adicionar inteligência biológica para identificar regiões regulatórias e estruturais.

- [x] Task: Implementar análise de região upstream (promotores) [checkpoint: 53a134a]
    - [x] Write Tests: Validar extração de sequência upstream e busca de TATA box
    - [x] Implement Feature: Criar função para análise de motivos upstream
- [x] Task: Implementar predição de Splice Sites (GT-AG) [checkpoint: 55072be]
    - [x] Write Tests: Testar identificação de sítios doadores e aceitadores
    - [x] Implement Feature: Adicionar lógica de predição de splicing básico
- [ ] Task: Conductor - User Manual Verification 'Phase 2: Bio-Inteligência' (Protocol in workflow.md)

## Phase 3: Domínios Proteicos e Relatório Final
Traduzir sequências e identificar domínios funcionais.

- [ ] Task: Implementar identificação de domínios proteicos
    - [ ] Write Tests: Testar tradução e busca de motivos proteicos conhecidos
    - [ ] Implement Feature: Integrar busca de padrões Prosite
- [ ] Task: Atualizar relatório para incluir novas análises
    - [ ] Write Tests: Validar novo formato do arquivo `orf_report.txt`
    - [ ] Implement Feature: Adicionar filtros e metadados ao relatório gerado
- [ ] Task: Conductor - User Manual Verification 'Phase 3: Finalização' (Protocol in workflow.md)
