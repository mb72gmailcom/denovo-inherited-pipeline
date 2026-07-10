# Inherited

Analyze rare variants in family trios from a VCF, precomputed gnomAD allele-frequency JSON, and a family relations file.

## Run

**Without installing** (clone and run from the project directory):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22
```

**Or install** so the `inherited` command is on your PATH:

```bash
pip install -e ".[dev]"
inherited analyze --vcf chr22.vcf.gz --af-json gnomad_chr22.json --family-file families.tsv -o results/chr22
```

## Family file format

Tab-separated file (required via `--family-file`):

```text
spid    family_id    father_id    mother_id
child1  fam1         fa1          ma1
fa1     fam1         0            0
ma1     fam1         0            0
```

Complete trios (`father_id` and `mother_id` both not `0`) are used for analysis.

## gnomAD AF JSON

```json
{
  "var_key": 0.0001,
  "22:12345:A:G": {"AF": 0.001, "AF_EUR": 0.0005}
}
```

Variants with AF above the threshold are skipped. Missing keys default to `0` and are kept.

## Usage

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22
```

Disable multiallelic per-ALT processing:

```bash
python run.py analyze --vcf chr22.vcf.gz --af-json gnomad.json --family-file families.tsv -o results/chr22 --no-multiallelic
```

Enable debug memory logging every 50,000 variants (custom interval with `--memory-block`):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad.json \
  --family-file families.tsv \
  -o results/chr22 \
  --debug \
  --memory-block 50000
```

Skip variants inside repetitive regions (one file per chromosome):

```bash
python run.py analyze \
  --vcf chr22.vcf.gz \
  --af-json gnomad_chr22.json \
  --family-file families.tsv \
  -o results/chr22 \
  --remove-repeats chr22_repeats.bed
```

Repeat files are whitespace-separated ``chrom start end`` rows with half-open intervals ``[start, end)``.

## Output

The output directory contains:

- `inherited_XXXXX.tsv` / `mendelian_bad_XXXXX.tsv` — segmented result files (when `--segment-size > 0`)
- `inherited.tsv` / `mendelian_bad.tsv` — single files when `--segment-size 0`
- `checkpoint.json` — resume point (updated after each segment)
- `stats_cumulative.json` — running totals (updated after each segment)
- `inherited_per_variant_segXXXXX.json` — per-segment variant counts (merged at end)
- `inherited_per_variant.json` — merged variant counts
- `inherited_per_person.json` — person id → number of inherited variants
- `mendelian_bad_per_gt.json` — `mother_gt:father_gt:child_gt` → count
- `stats.json` — final summary counts
- `params.json` — parameters used for this run

Each chromosome output directory (e.g. `results/chr22/`) contains its own `params.json` alongside the result files.

Result TSV columns (default `--short-format`):

```text
#CHROM  POS  ID  REF  ALT  PATIENTS
22      3000 .   A   G    child1;child2
```

Use `--no-short-format` for full genotype output:

```text
#CHROM  POS  ID  REF  ALT  TRIO_CALLS
22      3000 .   A   G    child1=0/1|0/0|0/1|30
```

Use `--block-size` (default `10000`) for in-memory buffer flushes within a segment.

Use `--segment-size` (default `1000000`) to split output into segment files. Set `--segment-size 0` to disable segmentation.

Resume after a crash:

```bash
python run.py analyze ... -o results/chr2 --resume
```

Requires an existing incomplete `checkpoint.json` and `--segment-size > 0`.

## Tests

```bash
pytest
```
