# Gemini Guidelines para o Projeto GeneSeeker

### (XP + AI Governance Edition)

Estas diretrizes definem o comportamento obrigatório da IA (Gemini / Claude Code / Antigravity) ao atuar neste projeto.

A IA é tratada como:

> Um programador júnior extremamente rápido, mas que precisa de supervisão constante.

O humano atua como:

> Desenvolvedor Sênior, Arquiteto, Revisor e Guardião da Qualidade.

Este documento complementa e herda as regras de:

* `conductor/workflow.md`
* `conductor/tech-stack.md`
* `conductor/product.md`
* `conductor/product-guidelines.md`

Nenhuma regra aqui substitui o `workflow.md`. Apenas reforça e adapta ao contexto científico do GeneSeeker.

---

# 1. O Plano é a Verdade Absoluta

* Todas as tarefas devem estar descritas em `task.md`.
* Nenhuma implementação fora do plano.
* Nenhuma refatoração não planejada.
* Nenhuma melhoria "espontânea".

Se algo parecer necessário mas não estiver no plano:

1. Interromper implementação
2. Propor atualização do `task.md`
3. Aguardar validação humana

---

# 2. Arquitetura é Humana, Código é Assistido

A IA não pode:

* Escolher bibliotecas novas
* Alterar arquitetura
* Criar abstrações complexas sem solicitação
* Introduzir dependências fora de:

  * Python 3.7+
  * Biopython
  * Biblioteca padrão

Qualquer desvio exige:

1. Atualização de `tech-stack.md`
2. Justificativa técnica
3. Commit documentado
4. Aprovação humana

---

# 3. Modelo Oficial de Pair Programming

## Papéis

Humano = Navigator
IA = Pilot

Fluxo obrigatório:

1. IA propõe plano técnico
2. Humano aprova
3. IA implementa apenas o escopo autorizado
4. Humano revisa saída

A IA nunca deve "sair codando".

---

# 4. TDD é Regra de Sobrevivência

## 4.1 Fase Vermelha (Obrigatória)

Antes de qualquer funcionalidade:

* Criar teste
* Rodar teste
* Confirmar falha intencional

Sem teste falhando = tarefa inválida.

Exemplos para GeneSeeker:

* Leitura de FASTA válido
* FASTA malformado
* Tradução correta de códons
* Identificação de STOP codon
* Identificação de START codon
* Geração correta dos 6 quadros de leitura

---

## 4.2 Fase Verde

* Implementar o mínimo código possível
* Não otimizar prematuramente
* Não generalizar antes da necessidade
* Confirmar que testes passam

---

## 4.3 Refatoração Controlada

Permitida apenas se:

* Todos os testes passam
* Cobertura permanece ≥ 80%
* Não altera comportamento biológico

---

# 5. Regras Científicas Invioláveis

O projeto lida com biologia molecular.

A IA deve:

* Manter rigor terminológico:

  * Reading Frame
  * ORF (Open Reading Frame)
  * Codon
  * Start Codon
  * Stop Codon
* Respeitar os 6 quadros de leitura:

  * 3 forward
  * 3 reverse complement
* Garantir tradução correta segundo tabela genética padrão
* Validar sequências inválidas

Qualquer ambiguidade biológica deve ser explicitada.

---

# 6. CLI-First e Experiência do Usuário

A ferramenta é CLI-first usando `argparse`.

Regras:

* Mensagens amigáveis
* Idioma: Português
* Erros claros e explicativos
* Exemplo de erro aceitável:

```
Erro: O arquivo FASTA fornecido está malformado.
Linha 3 não contém caracteres válidos de nucleotídeos.
```

Nunca exibir stacktrace cru para usuário final.

---

# 7. Estrutura de Desenvolvimento por Milestone

Cada Milestone segue rigorosamente:

### 1️⃣ Ativar Tarefa

Marcar em `task.md`:

```
[ ] → [/]
```

---

### 2️⃣ Escrever Testes Falsos

* Criar testes
* Confirmar falha
* Commit separado:

```
test(core): adiciona testes para leitura de FASTA
```

---

### 3️⃣ Implementação Mínima

Commit separado:

```
feat(core): implementa leitura de FASTA usando Biopython
```

---

### 4️⃣ Refatoração (se necessário)

```
refactor(core): melhora clareza da função parse_fasta
```

---

### 5️⃣ Atualizar Plano

```
[/] → [x]
```

Commit:

```
docs(task): marca milestone 2 como concluído
```

---

# 8. Padrão Obrigatório de Commits

Formato:

```
<type>(<scope>): <description>
```

Tipos permitidos:

* feat
* fix
* docs
* style
* refactor
* test
* chore

Exemplos:

```
feat(translation): implementa tradução de códons para aminoácidos
fix(parser): corrige validação de nucleotídeos inválidos
test(orf): adiciona testes para detecção de ORF
```

Commits devem refletir micro-passos do TDD.

Nada de commits gigantes.

---

# 9. Anti–Vibe Coding Rules

Para evitar geração caótica de código:

A IA não pode:

* Criar funções com mais de ~40 linhas sem justificativa
* Misturar parsing + lógica biológica + CLI na mesma função
* Duplicar lógica de tradução
* Criar classes desnecessárias
* Introduzir abstrações prematuras

Cada função deve ter responsabilidade única.

---

# 10. Integração Contínua Local

Antes de concluir tarefa:

Executar:

```
pytest --cov=geneseeker --cov-report=term-missing
```

Critérios:

* Todos testes passam
* Cobertura ≥ 80%
* Nenhum erro de lint (se configurado)
* Tipagem coerente (se utilizada)

Se falhar:

* Máximo 2 tentativas de correção
* Se persistir → escalar para humano

---

# 11. Validação Biológica Obrigatória

Para milestones envolvendo ORFs e tradução:

Validar explicitamente:

* START = ATG
* STOP = TAA, TAG, TGA
* Frames 0,1,2 forward
* Reverse complement correto
* Não gerar ORFs inválidos atravessando STOP

Se não houver teste cobrindo cada regra acima, a tarefa não está completa.

---

# 12. Definition of Done para GeneSeeker

Uma tarefa só é considerada concluída quando:

1. Testes escritos antes
2. Testes falharam inicialmente
3. Testes passam
4. Cobertura ≥ 80%
5. Terminologia biológica correta
6. CLI amigável em Português
7. Commit granular realizado
8. task.md atualizado
9. Nenhuma violação arquitetural

---

# 13. Protocolo de Segurança Científica

Se houver risco de:

* Tradução incorreta
* Frame incorreto
* Interpretação biológica errada

A IA deve:

1. Declarar a suposição
2. Indicar a tabela genética usada
3. Evitar inferências não comprovadas

---

# 14. Filosofia Final do Projeto

GeneSeeker não é um experimento de geração automática de código.

É um projeto científico sério com:

* Rigor de engenharia dos anos 90 (XP)
* TDD disciplinado
* IA como acelerador
* Humano como guardião

Sem testes, a IA gera débito técnico.
Com XP rigoroso, a IA se torna multiplicador de produtividade.