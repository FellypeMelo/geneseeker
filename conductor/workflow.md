# Fluxo de Trabalho (AI-XP / Akita-Driven)

Este projeto utiliza o framework **AI-XP (Artificially Intelligent eXtreme Programming)** com foco em rigor científico e governança de IA.

## 1. Ciclo TDD Agêntico (Red-Green-Refactor)

O ciclo de desenvolvimento é estritamente controlado para evitar dívida técnica e alucinações de IA:

### 🔴 Fase 1: RED (Teste Falho)
- **Agente**: Test Analyst Agent.
- **Ação**: Criar um teste que descreve o comportamento biológico desejado (ex: `test_orf_finder.py`).
- **Restrição**: Proibido modificar código de produção. O teste DEVE falhar comprovadamente no CI/Pytest.
- **PBT**: Incluir testes de propriedade (Hypothesis) para garantir invariantes (ex: toda ORF começa com ATG).

### 🟢 Fase 2: GREEN (Implementação Mínima)
- **Agente**: Implementation Agent.
- **Ação**: Escrever o código estritamente necessário para passar no teste.
- **Restrição**: Apenas o contexto do teste falho é fornecido. Proibido sumarizar código.

### 🔵 Fase 3: REFACTOR (Melhoria de Design)
- **Agente**: Refactoring Agent.
- **Ação**: Otimizar legibilidade, reduzir complexidade ciclomática e extrair funções auxiliares (máx. 15 linhas).
- **Garantia**: Todos os testes anteriores devem continuar passando.

## 2. Protocolo AI-Pilot / Human-Navigator

- **Humano (Navigator)**: Define intenções científicas, aprova algoritmos críticos e gerencia riscos biológicos.
- **IA (Pilot)**: Gera rascunhos, implementa mudanças atômicas, escreve testes e refatora sob supervisão.
- **Regra**: A IA nunca toma decisões arquiteturais sozinha.

## 3. Gestão de Contexto e Memória

- **Context Engineering**: Fornecer apenas os grafos AST dos módulos afetados para evitar "amnésia estrutural".
- **Hard Stop-Loss**: Se a IA entrar em loop de correção recursiva, reverter (`git checkout .`), expurgar histórico e reiniciar com prompt mais restrito.

## 4. Portões de Qualidade (Hard Blockers)

Um commit ou tarefa só é considerado concluído se:
- [ ] Todos os testes passarem (Pytest).
- [ ] Cobertura de testes ≥ 90%.
- [ ] Complexidade ciclomática ≤ 15 (radon).
- [ ] Tipagem estática validada (mypy).
- [ ] Nenhuma violação de Clean Architecture (import-linter).
- [ ] Audit trail gerado com justificativa algorítmica.

## 5. Procedimentos de Emergência (Hotfixes)

1. Criar branch `hotfix/`.
2. Escrever teste de regressão que capture o bug biológico.
3. Aplicar correção mínima.
4. Executar suite completa de CI antes do merge.
