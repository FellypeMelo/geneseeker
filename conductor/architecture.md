# GeneSeeker — Arquitetura do Sistema (Clean Architecture)

## 1. Visão Geral das Camadas

O GeneSeeker segue os princípios da **Clean Architecture**, garantindo que a lógica de domínio biológico seja isolada de detalhes de infraestrutura e interfaces.

```mermaid
graph TD
    subgraph Interface["Interface de Usuário (Interface)"]
        CLI[cli.py]
        MAIN[main.py]
    end

    subgraph Application["Casos de Uso (Application)"]
        PIPE[pipeline.py]
    end

    subgraph Domain["Núcleo Biológico (Domain)"]
        ORF[orf_finder.py]
        MOD[models.py: Sequence, Orf, Codon]
        TRN[translation.py]
    end

    subgraph Infrastructure["Infraestrutura (Infrastructure)"]
        SIO[sequence_io.py]
        REP[reporter.py]
    end

    %% Fluxo de Dependência
    Interface --> Application
    Application --> Domain
    Application --> Infrastructure
    Infrastructure -.-> Domain
```

## 2. Diagrama de Componentes Detalhado

```mermaid
graph TD
    subgraph INT["Interface / Interface"]
        M[main.py]
        C[cli.py]
    end

    subgraph APP["Aplicação / Application"]
        P[pipeline.py]
    end

    subgraph DOM["Domínio / Domain"]
        OF[orf_finder.py]
        MD[models.py]
        TL[translation.py]
    end

    subgraph INF["Infraestrutura / Infrastructure"]
        SI[sequence_io.py]
        RP[reporter.py]
    end

    FASTA[(Arquivo FASTA)]
    OUTPUT[(Relatórios)]

    M --> C
    C --> P
    P --> SI
    P --> OF
    P --> RP
    OF --> MD
    OF --> TL
    SI --> MD
    SI --> FASTA
    RP --> OUTPUT
```

## 3. Diagrama de Sequência — Pipeline de Análise

```mermaid
sequenceDiagram
    participant UI as Interface (cli.py)
    participant APP as Aplicação (pipeline.py)
    participant INF as Infra (sequence_io.py)
    participant DOM as Domínio (orf_finder.py)
    participant REP as Infra (reporter.py)

    UI->>APP: run_analysis(input_path, min_len)
    APP->>INF: read_fasta(input_path)
    INF-->>APP: Dict[id, DnaSequence]
    
    loop Para cada DnaSequence
        APP->>DOM: find_orfs(dna_sequence, min_len)
        DOM->>DOM: Analisa 6 frames
        DOM-->>APP: List[Orf]
    end

    APP->>REP: format_output(orfs, format)
    REP-->>APP: content
    APP->>REP: save_results(content, output_path)
    APP-->>UI: Sucesso
```

## 4. Modelagem de Domínio (Value Objects)

As entidades de domínio são implementadas como **Value Objects** imutáveis, garantindo integridade biológica desde a construção.

```mermaid
classDiagram
    class DnaSequence {
        +id: str
        +sequence: str
        +validate() bool
    }
    class Orf {
        +start: int
        +end: int
        +frame: int
        +strand: str
        +dna_sequence: str
        +protein_sequence: str
    }
    class Codon {
        +triplet: str
        +is_start() bool
        +is_stop() bool
    }
```

## 5. Leis de Dependência

1.  **Independência de Domínio**: O módulo `geneseeker.domain` não importa nada de `infrastructure`, `application` ou `interface`.
2.  **Inversão de Dependência**: A camada de aplicação orquestra os serviços de infraestrutura para alimentar o domínio.
3.  **Pureza Algorítmica**: `orf_finder.py` é uma função pura (ou conjunto de funções) que recebe dados e retorna dados, sem efeitos colaterais de I/O.
