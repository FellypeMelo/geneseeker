# Implementation Plan: Core ORF Finding and Translation (Milestone 3)

## Phase 1: Lógica de Domínio - Busca de ORFs (6 Frames)
Implementar o núcleo biológico para identificação de ORFs em todas as frentes de leitura.

- [ ] Task: Definir Value Objects `DnaSequence`, `Orf` e `Codon` (Domain)
    - [ ] Write Tests: Validar imutabilidade e integridade biológica (ATG start, STOP ends)
    - [ ] Implement: Criar `geneseeker/domain/models.py`
- [ ] Task: Implementar busca de ORFs em fita direta (Forward)
    - [ ] Write Tests: PBT para encontrar ORFs conhecidas em sequências sintéticas
    - [ ] Implement: `find_orfs` (O(n) por frame) em `geneseeker/domain/orf_finder.py`
- [ ] Task: Implementar busca de ORFs em fita reversa (Reverse)
    - [ ] Write Tests: Validar complemento reverso e varredura de frames negativos
    - [ ] Implement: Integrar `get_reverse_complement` no pipeline de busca
- [ ] Task: Adicionar filtragem por `min_len` e metadados (GC, Strand, Frame)
    - [ ] Write Tests: Testar limites de filtragem e precisão dos metadados
    - [ ] Implement: Adicionar lógica de filtragem e cálculo de GC
- [ ] Task: Conductor - User Manual Verification 'Phase 1: Busca de ORFs' (Protocol in workflow.md)

## Phase 2: Lógica de Domínio - Tradução Proteica (Standard Code)
Traduzir as ORFs encontradas usando a tabela genética padrão da NCBI.

- [ ] Task: Implementar mapeamento de códons (NCBI Table 1)
    - [ ] Write Tests: Validar tradução de todas as 64 trincas possíveis
    - [ ] Implement: Criar `geneseeker/domain/translation.py`
- [ ] Task: Integrar tradução no objeto `Orf`
    - [ ] Write Tests: PBT para garantir `len(prot) == len(dna)//3`
    - [ ] Implement: Atualizar `Orf` e `translate_orf`
- [ ] Task: Calcular composição de aminoácidos
    - [ ] Write Tests: Testar contagem de resíduos proteicos
    - [ ] Implement: Adicionar metadado `aa_composition` ao objeto `Orf`
- [ ] Task: Conductor - User Manual Verification 'Phase 2: Tradução' (Protocol in workflow.md)

## Phase 3: Infraestrutura - Relatórios (JSON, CSV, TXT)
Formatar e salvar os resultados das análises em múltiplos formatos.

- [ ] Task: Implementar formatadores JSON e CSV
    - [ ] Write Tests: Validar esquema JSON e colunas CSV
    - [ ] Implement: Criar `geneseeker/infrastructure/reporter.py`
- [ ] Task: Implementar formatador de texto (TXT) humanizado
    - [ ] Write Tests: Testar legibilidade do relatório amigável
    - [ ] Implement: Adicionar suporte a TXT no `reporter.py`
- [ ] Task: Integrar pipeline de saída (Application)
    - [ ] Write Tests: Testar fluxo completo da busca ao salvamento do arquivo
    - [ ] Implement: Criar `geneseeker/application/pipeline.py`
- [ ] Task: Conductor - User Manual Verification 'Phase 3: Relatórios' (Protocol in workflow.md)

## Phase 4: Interface e Validação Final
Expor a funcionalidade via CLI e garantir qualidade via Mutation Testing.

- [ ] Task: Atualizar CLI para suportar novos argumentos (format, min-len)
    - [ ] Write Tests: Testar parsing de argumentos via `main.py`
    - [ ] Implement: Atualizar `main.py`
- [ ] Task: Executar Validação de Cobertura e Mutmut
    - [ ] Write Tests: Garantir cobertura de testes ≥ 90%
    - [ ] Implement: Rodar `mutmut` e `radon` para auditoria final
- [ ] Task: Conductor - User Manual Verification 'Phase 4: Finalização' (Protocol in workflow.md)
