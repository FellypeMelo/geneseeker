# Princípios Gerais de Estilo de Código (AI-XP / GeneSeeker)

Estes princípios aplicam-se a todos os aspectos de desenvolvimento no projeto GeneSeeker.

## 1. Arquitetura Limpa (Clean Architecture First)
- **Prioridade**: O isolamento da lógica de domínio sobre detalhes técnicos é a regra número um.
- **Enforcement**: Se uma mudança violar a inversão de dependência, ela deve ser rejeitada imediatamente.

## 2. Rigor Científico e Algorítmico
- **Complexidade**: Prefira algoritmos $O(n)$ ou lineares sempre que possível.
- **Justificativa**: Cada decisão algorítmica deve ser baseada em fatos biológicos ou computacionais documentados.

## 3. TDD Inviolável
- **Fase RED Proved**: Nenhum código entra em produção sem um teste falhando primeiro.
- **PBT (Property-Based Testing)**: Obrigatório para capturar invariantes biológicos ocultos.

## 4. KISS + YAGNI
- **Simplicidade**: Soluções simples superam abstrações complexas.
- **Evite o Speculative Development**: Não implemente recursos para cenários que "podem" ocorrer; implemente apenas o que foi solicitado.

## 5. Legibilidade Humana > Engenhosidade
- Código deve ser compreensível para um par-programador (humano ou IA).
- Evite "clever code" (truques de linguagem) que prejudiquem a manutenibilidade a longo prazo.

## 6. Documentação do "Porquê"
- Documente a razão por trás de escolhas não óbvias.
- Use Audit Trail para rastrear decisões arquiteturais e algorítmicas críticas.
