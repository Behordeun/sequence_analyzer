# Sequence Analyzer

Bioinformatics toolkit for DNA/RNA sequence analysis, alignment, and phylogenetics.

Usable as a **Python library**, a **CLI tool**, or a **Streamlit web application**.

## Installation

Requires Python 3.11+.

```bash
# Install with uv (recommended)
uv sync                    # Core library only
uv sync --extra app        # Include Streamlit web UI
uv sync --all-extras       # Everything including dev tools

# Or install with pip
pip install .              # Core library
pip install ".[app]"       # With web UI
pip install ".[app,dev]"   # With dev tools
```

## Usage

### As a CLI tool

```bash
# Run sequence analysis
sequence-analyzer analyze --input sequences.fasta --output results.csv

# Align sequences (MSA, CLUSTAL output)
sequence-analyzer align --input sequences.fasta --method msa --format clustal

# Build a phylogenetic tree
sequence-analyzer tree --input sequences.fasta --method nj --output tree.nwk

# Launch the web UI
sequence-analyzer serve
```

### As a Python library

```python
from sequence_analyzer.io.parsers import parse_sequence_file
from sequence_analyzer.core.validation import validate_sequences
from sequence_analyzer.core.analysis import analyze_sequences

records = parse_sequence_file(open("sequences.fasta").read())
valid = validate_sequences(records, seq_type="DNA")
df = analyze_sequences(valid, is_rna=False)
print(df)
```

```python
# Contamination screening
from sequence_analyzer.core.contamination import detect_contamination

results = detect_contamination(valid, organism="escherichia_coli")
print(results[["ID", "Contamination_Risk", "Risk_Reason"]])
```

```python
# Variant calling against a reference
from sequence_analyzer.core.variants import call_variants, summarize_variants

reference = valid[0]
samples = valid[1:]
result = call_variants(reference, samples)
print(summarize_variants(result))
```

### Web application

```bash
sequence-analyzer serve
# Opens at http://localhost:8501
```

## Development

```bash
# Clone and setup
git clone https://github.com/Behordeun/sequence_analyzer.git
cd sequence_analyzer
uv sync --all-extras

# Run tests
uv run pytest

# Lint and format
uv run ruff check src/ tests/
uv run ruff format src/ tests/

# Type check
uv run mypy src/
```

## Project Structure

```text
src/sequence_analyzer/
├── cli.py              # CLI entry point (analyze, align, tree, serve)
├── core/               # Pure computation (no UI dependencies)
│   ├── alignment.py    # Pairwise + MSA
│   ├── analysis.py     # GC content, skew, composition
│   ├── contamination.py # Cross-species contamination detection
│   ├── genbank.py      # NCBI Entrez fetching
│   ├── motifs.py       # Regex motif scanning
│   ├── phylogenetics.py # Tree construction + bootstrapping
│   ├── qc.py          # Quality control assessment
│   ├── validation.py   # Sequence cleaning + type detection
│   └── variants.py    # Reference-based variant calling
├── models/             # Typed dataclasses for results
├── io/                 # File parsing (FASTA, PHYLIP, NEXUS)
└── app/                # Streamlit web interface
    ├── main.py         # Streamlit launcher
    └── pages/          # Multipage app
```

## Authors

- Muhammad Abiodun Sulaiman ([GitHub](https://github.com/Behordeun))
- Bolaji Fatai Oyeyemi ([GitHub](https://github.com/bollergene))

## License

MIT
