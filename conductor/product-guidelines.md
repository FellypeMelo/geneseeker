# GeneSeeker - Diretrizes do Produto

## Idioma e Comunicação
- **Idioma Principal**: Português. Toda a interface de linha de comando (CLI), mensagens de erro e documentação técnica devem ser redigidas em português de forma clara e profissional.
- **Tom de Voz**: Profissional, técnico e objetivo, adequado para um ambiente acadêmico e de pesquisa.

## Princípios de Experiência do Usuário (UX)
- **CLI-First**: O foco principal é a eficiência na linha de comando. O uso de argumentos (`argparse`) deve ser priorizado sobre a edição manual de scripts.
- **Feedback Informativo**: O usuário deve ser mantido informado sobre o progresso da análise (ex: barras de progresso para genomas grandes) e o status atual da execução.
- **Convenção sobre Configuração**: O GeneSeeker deve funcionar com configurações padrão sensatas, minimizando a necessidade de ajustes para tarefas comuns de identificação de ORFs.

## Tratamento de Erros
- **Foco no Usuário**: Mensagens de erro devem ser amigáveis e acionáveis, explicando o que deu errado em termos biológicos ou de sistema (ex: "Arquivo FASTA mal formatado" em vez de apenas um `ValueError`).
- **Robustez**: O programa deve validar as sequências de entrada antes de iniciar o processamento pesado.

## Estilo de Documentação
- **Acadêmico/Formal**: A documentação deve incluir explicações detalhadas sobre os algoritmos de busca de ORFs, referências biológicas e justificativas teóricas para as escolhas de implementação.
- **Rigor Científico**: Definições claras de termos como "Reading Frame", "Codon" e "ORF" devem estar presentes.
