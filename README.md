# Excluscreen v3

A desktop application for screening DNA primer sequences against the NCBI nucleotide database. It identifies off-target binding sites and flags cases where two primers bind to the same genomic region, which may produce unintended amplicons.

## Requirements

- Python 3
- Internet connection (queries NCBI BLAST)
- `requests` library (`pip install requests`)

## Setup

```
python -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
python excluscreen.py
```

## Usage

1. Paste primer sequences into the input box. Accepts FASTA format or plain `name sequence` pairs, one per line. You can also import a FASTA file.
2. Enter your email address (required by NCBI).
3. Adjust parameters if needed:
   - **Min Coverage %** (default 50): minimum fraction of the primer that must align.
   - **Min Identity %** (default 70): minimum sequence similarity in the alignment.
   - **Max Distance** (default 5000): maximum distance in bp between two hits on the same subject to be reported as a pair.
   - **Database**: NCBI database to search (`nt`, `nr`, `refseq_genomic`, `refseq_rna`).
   - **Exclude TaxIDs**: comma-separated NCBI taxonomy IDs to exclude from results.
   - **RefSeq Filter**: restrict results to curated reference sequences.
4. Click **Run Analysis** and wait. Each primer is submitted to NCBI BLAST individually. The status log tab shows progress.
5. Results appear in three tabs:
   - **Hit Pairs**: deduplicated cases where two primers bind the same subject sequence within the max distance. Shows orientation (convergent, divergent, same-strand), distance, identity, coverage, and mismatch counts. Convergent pairs are the most relevant for unintended amplification.
   - **Individual Hits**: all binding sites found per primer (top 30 patterns per primer).
   - **Status Log**: timestamped progress messages.
6. Click **Export to CSV** to save results.

## How it works

Each primer sequence is submitted to NCBI's BLASTN API using the `blastn-short` task (word size 7, e-value 1000, hit list size 5000, low-complexity filter off). The program polls NCBI at 30-second intervals until results are ready, then parses the returned XML.

Hits are filtered by the coverage and identity thresholds. The program then groups all hits by subject sequence ID and finds every pair of hits from different primers (or the same primer with itself) that land on the same subject within the specified max distance. For each pair, it calculates the distance between the hits and classifies the orientation:

- **Convergent**: primers face each other (could produce an amplicon)
- **Divergent**: primers face away from each other
- **Same-strand**: both on the same strand

Identical hit-pair patterns across different accessions are deduplicated and counted. The GUI runs BLAST operations on a background thread to keep the interface responsive.
