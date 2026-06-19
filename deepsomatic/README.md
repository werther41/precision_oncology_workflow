# DeepSomatic Somatic Calling Path

## What was run

GPU-accelerated DeepSomatic via Parabricks 4.7.0-1 (`pbrun deepsomatic`), reusing
pre-aligned BAMs from the Mutect2 pipeline — no fq2bam re-run.

| Parameter | Value |
|-----------|-------|
| Container | `nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1` |
| Model | WGS shortread (default; WES via `MODEL_TYPE=WES`; PacBio/ONT also supported) |
| Input BAMs | `HCC1395_tumor.bam` / `HCC1395_normal.bam` (no re-alignment) |
| Post-processing | PASS-only filter via `bcftools view -f PASS` — no Mutect2-era cascade |
| FILTER semantics | `PASS` = somatic call; `GERMLINE`/`RefCall`/`LowQual`/`NoCall` = suppressed |

## ⚠ Training-data caveat

**HCC1395 and COLO829 are both in DeepSomatic's training corpus.**

Benchmark numbers produced against the SEQC2 HCC1395 truth set (`eval/comparison.md`)
measure memorization on a training sample, not generalization to new patients.
**Do not cite HCC1395 F1 numbers as clinical validation evidence.**

For a generalization estimate, COLO829 held-out results will be placed in
`eval/eval_colo829/` once COLO829 BAMs are aligned and truth VCFs prepared
(`scripts/download_colo829_truth.sh`).

## Key differences from Mutect2

| Aspect | Mutect2 path | DeepSomatic path |
|--------|-------------|------------------|
| FILTER mechanism | FilterMutectCalls + 3 post-filters (POPAF, MMQ, TLOD sweep) | Single learned FILTER column |
| Germline resource | gnomAD AF VCF required | Not needed (learned) |
| Panel of Normals | 1000G PoN required | Not needed (learned) |
| Contamination model | GetPileupSummaries + CalculateContamination | Not needed |
| Orientation model | LearnReadOrientationModel | Not needed |
| Output sample columns | 2 (tumor + normal) | 1 (tumor only) |
| FFPE support | Partial (orientation model) | Via `--pb-model-file` custom FFPE model |
| vcfeval score field | `INFO/TLOD` | `FORMAT/GQ` |

## Scripts

- `scripts/run_deepsomatic_pipeline.sh` — Phases 0-2: preflight, DeepSomatic call, PASS extract
- `scripts/eval_deepsomatic.sh` — Phase 3: RTG vcfeval + comparison table
- `scripts/download_colo829_truth.sh` — Phase 4: COLO829 held-out truth prep
