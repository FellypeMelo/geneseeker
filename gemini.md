# Gemini Guidelines para o Projeto GeneSeeker

Estas diretrizes definem o comportamento esperado da IA (Gemini/Claude Code/Antigravity) ao atuar neste projeto, baseando-se nos princípios estabelecidos na pasta `conductor/` (especialmente `workflow.md`, `tech-stack.md`, `product.md` e `product-guidelines.md`).

## 1. O Plano é a Verdade Absoluta
- Todas as tarefas e implementações devem estar de acordo com o `task.md` (o equivalente a `plan.md` no nosso fluxo de IA).
- Nenhuma nova dependência ou tecnologia deve ser introduzida sem antes planejar e validar com o usuário. A stack principal é **Python 3.7+** e **Biopython**.

## 2. Princípios de Desenvolvimento Orientado a Testes (TDD)
- **Fase Vermelha (Failing Test):** Antes de implementar qualquer nova funcionalidade (ex: ler FASTA, traduzir aminoácidos), crie ou modifique um arquivo de teste e execute-o para provar que ele falha intencionalmente.
- **Fase Verde (Passing Test):** Escreva o menor código de aplicação possível para que os testes passem.
- **Isolamento de Commits:** Cada passo de desenvolvimento (testes falhando, testes passando, refatoração) deve idealmente ser acompanhado de commits granulares detalhando as mudanças, usando Conventional Commits (ex: `feat(core): implementa leitura de arquivos FASTA usando Biopython`).
- O suporte inicial aos testes deve ser montado antes da implementação para garantir a filosofia do projeto.

## 3. Padrão de Commits
O formato obrigatório para commits no repositório é:
```
<type>(<scope>): <description>

[optional body]
```
No qual `<type>` pode ser `feat`, `fix`, `docs`, `style`, `refactor`, `test`, ou `chore`.
Todo o histórico do Git deve refletir de forma clara cada passo do progresso dos milestones.

## 4. Experiência do Usuário e CLI-First
- A ferramenta está sendo moldada para ser usada via interface de linha de comando (`argparse`). O feedback deve ser amigável, no idioma **Português**, fornecendo detalhes sobre erros (como FASTA malformado), status ou opções não reconhecidas.
- O rigor científico em relação aos termos biológicos ("Reading Frame", "ORF", "Codon") deve ser mantido nas documentações e nos logs.

## 5. Passos para Cada Nova Tarefa
1. **Ativar Task/Milestone:** Marcar no plano (`task.md`) a tarefa atual como "Em Progresso" `[/]`.
2. **Escrever Testes Falsos:** Verificar o arquivo de testes; criar testes para os novos requisitos e observar falhar.
3. **Desenvolver o Código:** Implementar os Milestones (atualmente Milestone 2, leitura FASTA, fitas complementares e conversão de aminoácidos).
4. **Validar:** Certificar que os testes passam e todas as lógicas biológicas (6 quadros de leitura, STOP codons, START codons) atuam conforme o esperado.
5. **Realizar Commit Sistemático:** Concluir com um commit bem formatado documentando a fase terminada.
6. **Atualizar o Plano:** Marcar como concluído `[x]` no arquivo de acompanhamento.

Seguir rigorosamente o fluxo de "TDD -> Commit Sistemático -> Milestone Update" é essencial para manter o código seguro e alinhado aos princípios definidos em `conductor/workflow.md`.
