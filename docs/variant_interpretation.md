# Variant Interpretation Reference

## AMP/ASCO/CAP 2017 Tier System

The gold standard for somatic variant clinical classification.

### Tier I — Strong clinical significance
- **IA**: Variant has FDA-approved therapy in this tumor type, OR is included in professional guidelines (NCCN) as predictive of response
- **IB**: Well-powered studies with expert consensus

### Tier II — Potential clinical significance
- **IIC**: FDA-approved therapy for a *different* tumor type, OR investigational therapies with clinical evidence in this tumor type
- **IID**: Preclinical trials or limited case reports

### Tier III — Variants of unknown significance (VUS)
Not seen in cancer databases, no functional data, no clear consequence prediction

### Tier IV — Benign / likely benign
Common population variants, synonymous changes, deep intronic without splice impact

## OncoKB → AMP Tier mapping (used in this pipeline)

| OncoKB Level | Meaning | AMP Tier |
|---|---|---|
| LEVEL_1 | FDA-recognized biomarker, this indication | IA |
| LEVEL_2 | Standard care biomarker, this indication | IB |
| LEVEL_3A | Clinical evidence, this indication | IIC |
| LEVEL_3B | Clinical evidence, other indication | IIC |
| LEVEL_4 | Compelling biological evidence | IID |
| LEVEL_R1 | Standard-care resistance | IA (resistance) |
| LEVEL_R2 | Investigational resistance | IID (resistance) |

## Key driver genes by tumor type (worth understanding before reviewing reports)

| Tumor | Drivers | Actionable |
|---|---|---|
| **LUAD** (lung adeno) | EGFR, KRAS, ALK, ROS1, BRAF, MET, RET, ERBB2 | Most have FDA-approved TKIs |
| **CRC** (colorectal) | KRAS, NRAS, BRAF, PIK3CA, APC, TP53 | KRAS/NRAS wt → cetuximab; BRAF V600E → encorafenib |
| **Breast** | PIK3CA, ESR1, ERBB2, BRCA1/2 | PIK3CA → alpelisib; BRCA → olaparib |
| **Melanoma** | BRAF, NRAS, KIT, NF1 | BRAF V600 → vemurafenib + cobimetinib |
| **AML** | FLT3, NPM1, IDH1/2, TP53 | FLT3 → midostaurin; IDH → ivosidenib/enasidenib |

## Variant nomenclature cheat sheet

- `p.L858R` — protein-level: leucine 858 → arginine (HGVSp short)
- `c.2573T>G` — coding-DNA-level
- `chr7:55191822 T>G` — genomic, GRCh38
- `NM_005228.5:c.2573T>G` — full HGVS with transcript

All four can refer to the same variant: the canonical EGFR L858R activating mutation.

## What makes a variant "actionable"

A variant must meet several criteria to drive treatment decisions:

1. **Detected reliably** — VAF ≥ 5%, alt depth ≥ 5 reads, tumor purity sufficient
2. **Functionally consequential** — protein-altering, splice-disrupting, or known regulatory
3. **Clinically validated** — published evidence linking it to a therapy response
4. **Patient context fits** — diagnosis, prior treatments, performance status, comorbidities

The pipeline can identify (1) and (2) computationally; (3) and (4) require expert review by a molecular tumor board.
