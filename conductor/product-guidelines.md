# GeneSeeker - Diretrizes de Engenharia e Produto

## Leis Invioláveis (Iron Laws)

Este projeto opera sob rigor técnico absoluto, regido pelas seguintes leis:

1.  **TDD Mandatório**: Nunca modificar código sem um teste falhando primeiro. PBT (Property-Based Testing) é obrigatório para invariantes biológicos.
2.  **Clean Architecture**: Camada de Domínio nunca importa infraestrutura. Módulos biológicos devem ser puros e testáveis.
3.  **Zero Vibe Coding**: Decisões algorítmicas são justificadas por biologia e computação; nada é aceito sem compreensão profunda.
4.  **KISS + YAGNI**: Proibido antecipar recursos não solicitados; funções com máx. 15 linhas lógicas.
5.  **Rigor Matemático**: Algoritmos de busca devem ter complexidade $O(n)$ por quadro de leitura.
6.  **Value Objects**: Usar classes tipadas (`Sequence`, `Orf`) em vez de strings brutas.

## Idioma e Comunicação
- **Idioma Principal**: Português (CLI, Erros, Documentação).
- **Tom de Voz**: Profissional, científico e técnico.

## Princípios de Experiência do Usuário (UX)
- **CLI-First**: Eficiência absoluta na linha de comando via `argparse`.
- **Feedback Científico**: Relatórios detalhados com metadados (start, end, frame, strand).
- **Tratamento de Erros**: Mensagens amigáveis para falhas de entrada (ex: "Bases inválidas encontradas em 'seq_1'").

## Estilo de Documentação
- **Acadêmico/Formal**: Incluir referências biológicas e justificativas de design.
- **AuditTrail**: Cada mudança deve ser rastreável a um requisito ou decisão algorítmica.
