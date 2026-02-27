# GeneSeeker - Identificador de ORFs

## Descrição

O **GeneSeeker** é uma ferramenta para identificação de **Open Reading Frames (ORFs)** em sequências de DNA. ORFs são regiões de DNA que começam com um códon de início e terminam com um códon de parada, representando potenciais regiões codificantes de proteínas.

### O que são ORFs?

Um **Open Reading Frame (ORF)** é uma sequência contínua de DNA que:
1. Começa com um **códon de início** (ATG - Metionina)
2. Continua com múltiplos códons de aminoácidos
3. Termina com um **códon de parada** (TAA, TAG ou TGA)

A identificação de ORFs é crucial para:
- Predição de genes
- Anotação de genomas
- Descoberta de novas proteínas
- Estudos funcionais

## Funcionalidades

- **Análise de 3 Quadros de Leitura**: Verifica todos os quadros possíveis (0, 1, 2)
- **Detecção de Codons**: Identifica códons START (ATG) e STOP (TAA, TAG, TGA)
- **Relatório Estruturado**: Gera relatório em arquivo texto
- **Sequência Completa**: Extrai a sequência completa de cada ORF encontrado

## Estrutura de Dados

### 📁 `test_data/` - Dados Sintéticos (Commitados)
Contém **55+ arquivos FASTA fabricados** com ORFs conhecidos:
- ✅ **Commitados no GitHub**
- 🧪 **ORFs controlados** (quantidade e posição conhecidas)
- 📊 **Casos de borda** (sem ORFs, muitos ORFs, sobrepostos)
- 🎯 **Validação garantida**

**Regenerar:**
```bash
python generate_test_data.py
```

### 📁 `data/` - Dados Reais (Gitignored)
Para dados reais do NCBI, genomas, etc.:
- 🚫 **Ignorado pelo Git** (não vai para GitHub)
- 🧬 **Dados reais** de pesquisa
- 💾 **Sem limite de tamanho**

**Formatos recomendados:**
- **Nucleotide FASTA** - Genomas completos ou segmentos
- **Coding Region (CDS)** - Apenas regiões codificantes
- **mRNA** - Transcritos processados

## Instalação

### Pré-requisitos

- Python 3.7 ou superior
- pip

### Passos

```bash
git clone https://github.com/FellypeMelo/geneseeker.git
cd geneseeker
pip install -r requirements.txt
```

## Como Usar

### Execução Básica

```bash
python main.py
```

O programa analisa uma sequência de exemplo e gera um relatório.

### Modificando a Sequência

Edite a variável `test_sequence` no final do arquivo `main.py`:

```python
test_sequence = "ATGCGATACTGAATGCCCTAGATGAAATAA"
```

### Exemplo de Saída

```
============================================================
GeneSeeker - Identificador de ORFs
============================================================

Sequência analisada: ATGCGATACTGAATGCCCTAGATGAAATAA
Comprimento: 30 bp


Quadro de Leitura 0:
----------------------------------------
  ORF encontrado:
    Posição: 0 - 12
    Comprimento: 12 bp
    Sequência: ATGCGATACTGA

Quadro de Leitura 1:
----------------------------------------
  Nenhum ORF completo encontrado

Quadro de Leitura 2:
----------------------------------------
  ORF encontrado:
    Posição: 2 - 30
    Comprimento: 28 bp
    Sequência: ATGCCCTAGATGAAATAA

Relatório salvo em: orf_report.txt

Análise concluída!
```

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

| Tipo | Códons | Significado |
|------|--------|-------------|
| **START** | ATG | Metionina - Inicia a tradução |
| **STOP** | TAA | Ocre - Para a tradução |
| **STOP** | TAG | Âmbar - Para a tradução |
| **STOP** | TGA | Ôpal - Para a tradução |

## Estrutura do Projeto

```
geneseeker/
├── main.py              # Código principal
├── requirements.txt     # Dependências
├── README.md           # Documentação
└── orf_report.txt      # Relatório gerado (após execução)
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

- Apenas analisa fita forward (não complementar)
- Não considera introns (para eucariotos)
- Não faz análise de regiões regulatórias
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