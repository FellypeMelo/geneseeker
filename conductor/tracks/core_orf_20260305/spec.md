# Specification: Core ORF Finding and Translation (Milestone 3 Core)

## Overview
This track implements the primary biological logic of GeneSeeker, enabling the identification of Open Reading Frames (ORFs) across all six frames (3 forward + 3 reverse) and their translation into amino acid sequences.

## Functional Requirements
- **ORF Identification**:
    - Identify ORFs starting with `ATG` and ending with `TAA`, `TAG`, or `TGA`.
    - Support for 6 reading frames (3 forward strand, 3 reverse strand using reverse complement).
    - Filter ORFs by a configurable `min_len` parameter (default: 100 bp).
- **Protein Translation**:
    - Use the **Standard Genetic Code (NCBI Table 1)** for DNA to amino acid conversion.
- **Metadata Extraction**:
    - For each ORF, capture:
        - `sequence_id`: FASTA ID.
        - `start`, `end`: Genomic coordinates (0-based).
        - `strand`: `+` (forward) or `-` (reverse).
        - `frame`: `0`, `1`, or `2`.
        - `gc_content`: Percentage of G+C.
        - `aa_composition`: Brief breakdown of amino acid counts.
        - `dna_sequence`: DNA sequence.
        - `translation`: Resulting protein sequence.
- **Reporting**:
    - Generate output in **JSON**, **CSV**, and **Text (TXT)** formats.

## Non-Functional Requirements
- **Algorithm Performance**: Implement search with $O(n)$ complexity per frame.
- **Architecture**: Follow the 4-layer Clean Architecture (Domain, Application, Infrastructure).
- **Testing**:
    - TDD approach: Write failing tests before implementation.
    - PBT (Property-Based Testing) using Hypothesis for biological invariants (e.g., all ORFs must start with ATG and end with STOP).

## Acceptance Criteria
1. `find_orfs` identifies 100% of canonical ORFs in the `test_data/` suite.
2. `translate_orf` correctly converts sequences using NCBI Table 1.
3. Reports are generated in all three requested formats.
4. Mutation testing (mutmut) score > 90% for domain logic.

## Out of Scope
- Advanced regulatory analysis (promoter motifs, splicing).
- Integration with external databases (BLAST).
