# GeneSeeker - Identificador de ORFs

## Descrição

O **GeneSeeker** é uma ferramenta de bioinformática de linha de comando para identificação de **Open Reading Frames (ORFs)** em sequências de DNA. Análisa todos os **6 quadros de leitura** (3 na fita direta + 3 na fita reversa), filtra por tamanho mínimo, traduz as ORFs encontradas para sequências de aminoácidos e gera relatórios em múltiplos formatos.

### O que são ORFs?

Um **Open Reading Frame (ORF)** é uma sequência contínua de DNA que:

1. Começa com um **códon de início** (ATG - Metionina)
2. Continua com múltiplos códons de aminoácidos
3. Termina com um **códon de parada** (TAA, TAG ou TGA)

A identificação de ORFs é fundamental para:

- Predição de genes
- Anotação de genomas
- Descoberta de novas proteínas
- Estudos funcionais

## Arquitetura do Sistema

```
main.py (CLI / Orquestrador do Pipeline)
    │
    ├── geneseeker/sequence.py   — Leitura e validação de arquivos FASTA
    ├── geneseeker/utils.py      — Utilitários: limpeza e complemento reverso
    ├── geneseeker/orf_finder.py — Identificação de ORFs nos 6 quadros + tradução
    └── geneseeker/reporter.py   — Formatação (text/csv/json), exibição e salvamento
```

O pipeline segue 4 etapas sequenciais:

```
[Entrada FASTA] → load_sequences() → run_analysis() → output_results() → [Saída]
```

| Etapa                            | Função             | Responsável      |
| -------------------------------- | ------------------ | ---------------- |
| Leitura e validação              | `load_sequences()` | sequence + utils |
| Identificação de ORFs + tradução | `run_analysis()`   | orf_finder       |
| Formatação e saída               | `output_results()` | reporter         |

## Funcionalidades

- **Análise de 6 Quadros de Leitura**: Verifica todos os frames da fita direta (0, 1, 2) e da fita reversa (0, 1, 2 no complemento reverso)
- **Filtragem por Tamanho**: Descarta ORFs menores que o limite mínimo configurável (padrão: 100 bp)
- **Detecção de Códons**: Identifica códons START (ATG) e STOP (TAA, TAG, TGA)
- **Tradução Automática**: Converte a sequência de DNA de cada ORF em sequência de aminoácidos
- **Múltiplos Formatos de Saída**: Relatórios em texto legível, CSV ou JSON
- **Suporte a Múltiplas Sequências**: Processa arquivos FASTA com múltiplas entradas de uma vez
- **Análise Regulatória**: Identificação de motivos de promotores upstream (TATA box, Pribnow box) _(em desenvolvimento)_
- **Predição de Sítios de Splicing**: Detecção de junções canônicas GT-AG _(em desenvolvimento)_
- **Anotação Funcional**: Identificação básica de domínios proteicos (ex: Zinc Finger) _(em desenvolvimento)_

## Estrutura de Módulos

### `geneseeker/sequence.py`

| Função                        | Descrição                                                         |
| ----------------------------- | ----------------------------------------------------------------- |
| `read_fasta(file_path)`       | Lê um arquivo FASTA e retorna `{id: sequencia}`                   |
| `validate_sequence(sequence)` | Valida se a sequência contém apenas bases válidas (A, C, G, T, N) |

### `geneseeker/utils.py`

| Função                             | Descrição                                                    |
| ---------------------------------- | ------------------------------------------------------------ |
| `clean_sequence(sequence)`         | Remove espaços, quebras de linha e normaliza para maiúsculas |
| `get_reverse_complement(sequence)` | Gera o complemento reverso de uma fita de DNA (5'→3')        |

### `geneseeker/orf_finder.py`

| Função                            | Descrição                                                        |
| --------------------------------- | ---------------------------------------------------------------- |
| `find_orfs(sequence, min_length)` | Encontra todas as ORFs nos 6 quadros, filtrando por `min_length` |
| `translate_orf(orf_sequence)`     | Traduz uma sequência de DNA de ORF para aminoácidos              |

### `geneseeker/reporter.py`

| Função                                        | Descrição                                                |
| --------------------------------------------- | -------------------------------------------------------- |
| `display_summary(orfs)`                       | Imprime resumo (total, maior, menor ORF) no terminal     |
| `format_output(orfs, format_type)`            | Formata a lista de ORFs em `'text'`, `'csv'` ou `'json'` |
| `save_results(formatted_output, output_path)` | Salva o resultado formatado em arquivo                   |

## Estrutura de Dados

### 📁 `test_data/` — Dados Sintéticos (Commitados)

Contém **64 arquivos FASTA fabricados** com ORFs conhecidos:

- ✅ Commitados no GitHub
- 🧪 ORFs controlados (quantidade e posição conhecidas)
- 📊 Casos de borda: sem ORFs, múltiplas ORFs, sobreposições, sequências curtas e longas
- 🎯 Cobertura por frame: arquivos específicos para frame 0, 1 e 2

**Regenerar os dados de teste:**

```bash
python generate_test_data.py
```

**Nomenclatura dos arquivos:**

```
geneseeker_test_<N>_<descricao>.fasta
  N: número do caso de teste (01-55+)
  descricao: tipo do caso (ex: 3orfs, no_orfs, multi_4orfs, short, long, frame0_only)
```

### 📁 `data/` — Dados Reais (Gitignored)

Para dados reais do NCBI, genomas, etc.:

- 🚫 Ignorado pelo Git (não vai para o GitHub)
- 🧬 Dados reais de pesquisa
- 💾 Sem limite de tamanho

**Formatos recomendados:**

- **Nucleotide FASTA** — Genomas completos ou segmentos
- **Coding Region (CDS)** — Apenas regiões codificantes
- **mRNA** — Transcritos processados

## Instalação

### Pré-requisitos

- Python 3.7 ou superior
- pip

### Passos

```bash
git clone https://github.com/FellypeMelo/geneseeker.git
cd geneseeker
python -m venv venv
source venv/bin/activate      # Linux/macOS
# venv\Scripts\activate       # Windows
pip install -r requirements.txt
```

## Como Usar

### Execução Básica (CLI)

```bash
python main.py <arquivo.fasta>
```

### Argumentos disponíveis

| Argumento   | Tipo       | Padrão | Descrição                                     |
| ----------- | ---------- | ------ | --------------------------------------------- |
| `input`     | posicional | —      | Caminho para o arquivo FASTA de entrada       |
| `--output`  | opcional   | stdout | Caminho para o arquivo de saída               |
| `--min-len` | opcional   | `100`  | Tamanho mínimo de ORF em pares de bases       |
| `--format`  | opcional   | `text` | Formato do relatório: `text`, `csv` ou `json` |

### Exemplos

```bash
# Análise simples com saída no terminal
python main.py test_data/geneseeker_test_01_3orfs.fasta

# Filtro de tamanho mínimo de 150 bp com saída em CSV
python main.py genome.fasta --min-len 150 --format csv --output resultado.csv

# Saída em JSON
python main.py test_data/geneseeker_test_16_multi_4orfs.fasta --format json --output orfs.json
```

### Exemplo de Saída (formato `text`)

```
============================================================
GeneSeeker - Identificador de ORFs
============================================================
Total de ORFs encontradas: 3
Maior ORF: 120 bp | Menor ORF: 102 bp

Quadro de Leitura 0 (fita direta):
  ORF #1:
    Posição:     0 – 120
    Comprimento: 120 bp
    Tradução:    MXXX...
```

## Testes

O projeto segue **TDD (Test-Driven Development)**. Os testes ficam em `tests/`.

```bash
pytest tests/
```

**Cobertura dos testes (`tests/test_analysis.py`):**

- Filtragem de ORFs por tamanho mínimo (`test_find_orfs_with_min_len`)
- Análise de motivos de promotores upstream (`test_analyze_promoters`)
- Predição de sítios de splicing GT-AG (`test_predict_splice_sites`)
- Identificação de domínios proteicos Zinc Finger (`test_identify_protein_domains`)

## Dependências

| Pacote      | Versão | Uso                                                   |
| ----------- | ------ | ----------------------------------------------------- |
| `biopython` | 1.81   | Manipulação de sequências biológicas e arquivos FASTA |

## Conceitos Importantes

### Quadros de Leitura (Reading Frames)

O DNA pode ser lido em 3 quadros diferentes, dependendo de onde começamos:

```
Sequência: ATGCGATACTGA

Quadro 0: ATG CGA TAC TGA  (encontra START e STOP)
Quadro 1:  TGC GAT ACT GA   (não começa com ATG)
Quadro 2:   GCG ATA CTG A    (não começa com ATG)
```

### Codons Importantes

| Tipo      | Códons | Significado                   |
| --------- | ------ | ----------------------------- |
| **START** | ATG    | Metionina - Inicia a tradução |
| **STOP**  | TAA    | Ocre - Para a tradução        |
| **STOP**  | TAG    | Âmbar - Para a tradução       |
| **STOP**  | TGA    | Ôpal - Para a tradução        |

## Estrutura do Projeto

```
geneseeker/                        (raiz do repositório)
├── main.py                        # CLI / Orquestrador do Pipeline
├── requirements.txt               # Dependências (biopython==1.81)
├── README.md                      # Documentação
├── LICENSE                        # Licença MIT
├── .gitignore                     # Regras de exclusão do Git
├── generate_test_data.py          # Script para gerar dados de teste sintéticos
├── gemini.md                      # Configuração do agente Gemini
├── geneseeker/                    # Pacote principal
│   ├── __init__.py                # Inicializador do pacote
│   ├── sequence.py                # Leitura e validação de FASTA
│   ├── utils.py                   # Utilitários (limpeza, complemento reverso)
│   ├── orf_finder.py              # Identificação de ORFs + tradução
│   └── reporter.py                # Formatação e salvamento de resultados
├── conductor/                     # Documentação de arquitetura e guias de desenvolvimento
│   ├── architecture.md            # Modelagem do sistema (diagramas Mermaid)
│   ├── product.md                 # Visão do produto
│   ├── product-guidelines.md      # Diretrizes do produto
│   ├── tech-stack.md              # Stack tecnológica
│   ├── workflow.md                # Fluxo de trabalho
│   ├── tracks/                    # Tracks de desenvolvimento
│   └── code_styleguides/          # Guias de estilo de código
├── tests/                         # Testes automatizados (pytest)
│   └── test_analysis.py           # Testes de análise de ORFs
├── test_data/                     # 64 arquivos FASTA sintéticos (commitados)
└── data/                          # Dados reais do NCBI (gitignored)
```

## Guia de Desenvolvimento

### Milestones do Projeto

#### Milestone 1: Detecção Básica ✅

- [x] Identificar ORFs em 3 quadros de leitura
- [x] Detectar códons START e STOP
- [x] Gerar relatório simples
- [x] Documentação inicial

#### Milestone 2: Melhorias de Funcionalidade ✅

- [x] Ler sequências de arquivos FASTA
- [x] Tradução de ORFs para sequências de aminoácidos
- [x] Análise de ambas as fitas (forward e reverse)
- [x] Seis quadros de leitura (3 forward + 3 reverse)

#### Milestone 3: Filtros e Análises ✅

- [x] Filtrar ORFs por tamanho mínimo
- [x] Análise de região upstream (promotores)
- [x] Predição de splice sites (para eucariotos)
- [x] Identificação de domínios proteicos

#### Milestone 4: Integração e Automação 🔄

- [ ] Pipeline automatizado com FastaFlow
- [ ] Comparação com bancos de dados (BLAST)
- [ ] Anotação funcional automática
- [ ] Interface gráfica

### Tarefas para Contribuidores

**Nível Iniciante:**

1. Adicionar argumentos de linha de comando para input/output
2. Implementar filtro de tamanho mínimo de ORF
3. Criar testes unitários

**Nível Intermediário:**

1. Implementar tradução para aminoácidos
2. Adicionar análise da fita reversa complementar
3. Criar visualização dos ORFs

**Nível Avançado:**

1. Integrar com banco de dados de proteínas
2. Implementar algoritmos de predição de genes (HMM)
3. Criar pipeline completo de anotação

## Algoritmo

O algoritmo percorre a sequência em passos de 3 nucleotídeos:

```python
for frame in [0, 1, 2]:
    i = frame
    while i < len(sequence) - 2:
        codon = sequence[i:i+3]

        if codon == "ATG":  # START encontrado
            procura_codon_stop()  # A partir da posição i+3

        i += 3
```

## Exemplos de Aplicação

### 1. Anotação de Genomas

Identificar genes em sequências genomicas recém-sequenciadas.

### 2. Descoberta de Novos Genes

Encontrar ORFs não anotados em genomas conhecidos.

### 3. Estudos Comparativos

Comparar ORFs entre espécies para estudos evolutivos.

### 4. Engenharia de Proteínas

Identificar regiões codificantes para modificação.

## Limitações Atuais

- Identificação de domínios proteicos é limitada a padrões básicos (ex: Zinc Finger)
- Não integra com bancos de dados externos (BLAST)
- Não possui interface gráfica
- Sequências codificadas são apenas indicativas

## Próximos Passos Recomendados

1. **Integração FastaFlow**: Pipeline automatizado
2. **Integração DB (BLAST)**: Comparar com bancos de dados
3. **Anotação Funcional**: Relatório automatizado
4. **Interface Gráfica**: GUI/Web dashboard para ORFs

## Conceitos Relacionados

### Código Genético

O código genético é degenerado (mais de um códon pode codificar o mesmo aminoácido).

### Fita Reversa Complementar

O DNA é dupla-fita. A fita complementar é lida no sentido 5'→3' reverso.

### Eucariotos vs Procariotos

- **Procariotos**: ORFs geralmente contínuos
- **Eucariotos**: ORFs podem ser interrompidos por introns

## Referências

- [Biopython Tutorial](https://biopython.org/wiki/Documentation)
- [Reading Frames](https://en.wikipedia.org/wiki/Reading_frame)
- [Genetic Code](https://en.wikipedia.org/wiki/Genetic_code)
- [Gene Prediction](https://en.wikipedia.org/wiki/Gene_prediction)
- [ORF Finder - NCBI](https://www.ncbi.nlm.nih.gov/orffinder/)

## Licença

MIT License - veja arquivo LICENSE

## Contato

Abra uma issue para dúvidas ou sugestões.

---

**Status**: 🟢 Funcional - Pronto para uso e expansão
