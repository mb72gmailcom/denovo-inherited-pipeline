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

Tab-separated file with a header row (required via `--family-file`):

```text
spid    sfid    father    mother    sex
child1  fam1    fa1       ma1       Female
boy1    fam1    fa1       ma1       1
fa1     fam1    0         0         Male
ma1     fam1    0         0         2
```

Required columns: `spid`, `sfid`, `father`, `mother`, `sex`. Extra columns are ignored.

Complete trios (`father` and `mother` both not `0`) with a recognized sex are used for analysis. Sex values: males `1` / `Male` / `male`; females `2` / `Female` / `female`. Children with missing or unrecognized sex are skipped on all chromosomes.

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

## Chromosome X

VCFs are expected to be split per chromosome. Contigs `X` and `chrX` use sex-aware logic:

- Female children: same diploid QC and inheritance rules as autosomes; one output bucket
- Male PAR1 / PAR2 (hg38 closed intervals `[10001, 2781479]` and `[155701383, 156030895]`): same diploid trio logic as autosomes
- Male nonPAR: mother–son pairs only (father ignored on chrX); child uses haploid QC (`DP≥5`, `AB≥0.85`, `GQ≥20`); mother uses diploid QC; inherited if `ac>0` and `mac>0`; mendelian_bad if `ac>0` and `mac==0`

## Chromosome Y

Contigs `Y` and `chrY` are treated as nonPAR only. Analysis loops over male children and father–son pairs:

- Child and father both use haploid QC (`DP≥5`, `AB≥0.85`, `GQ≥20`)
- Inherited if `ac>0` and `fac>0`; mendelian_bad if `ac>0` and `fac==0`
- Father haploid QC failure skips the pair
- Output uses only the `male_nonPar` bucket (no empty PAR/female files)
- Full format: `child_id=father_gt|child_gt|child_gq`

## Output

Autosomal output directory contains:

- `inherited_XXXXX.tsv` / `mendelian_bad_XXXXX.tsv` — segmented result files (when `--segment-size > 0`)
- `inherited.tsv` / `mendelian_bad.tsv` — single files when `--segment-size 0`
- `inherited_per_variant.json`, `inherited_per_person.json`, `mendelian_bad_per_gt.json`, `stats.json`

chrX output uses sex/region buckets:

- `inherited_females_XXXXX.tsv`, `inherited_male_par1_XXXXX.tsv`, `inherited_male_nonPar_XXXXX.tsv`, `inherited_male_par2_XXXXX.tsv`
- matching `mendelian_bad_*` files
- matching per-bucket `inherited_per_person_*.json`, `mendelian_bad_per_gt_*.json`, `inherited_per_variant_*.json`, `stats_*.json`

Shared run metadata:

- `checkpoint.json` — resume point (updated after each segment)
- `stats_cumulative.json` — running totals (updated after each segment)
- `cumulative_detail_deltas.jsonl` — append-only per-segment summary details used for resume
- `stats.json` — final overall summary counts
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

Male nonPAR chrX full format omits the father genotype: `child_id=mother_gt|child_gt|child_gq`.

Male nonPAR chrY full format omits the mother genotype: `child_id=father_gt|child_gt|child_gq`.

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
