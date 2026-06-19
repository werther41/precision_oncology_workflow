// Parabricks-based precision oncology pipeline
// Reproducible execution via Nextflow + Docker/Singularity
//
// Run:
//   nextflow run main.nf -profile gpu \
//     --sample_id PATIENT001 \
//     --tumor_r1 data/tumor_R1.fq.gz --tumor_r2 data/tumor_R2.fq.gz \
//     --normal_r1 data/normal_R1.fq.gz --normal_r2 data/normal_R2.fq.gz \
//     --ref /refs/Homo_sapiens_assembly38.fasta

nextflow.enable.dsl = 2

params.sample_id   = null
params.tumor_r1    = null
params.tumor_r2    = null
params.normal_r1   = null
params.normal_r2   = null
params.ref         = "/data/reference/Homo_sapiens_assembly38.fasta"
params.known_sites = "/data/reference/Homo_sapiens_assembly38.known_indels.vcf.gz"
params.pon         = "/data/reference/1000g_pon.hg38.vcf.gz"
params.germline    = "/data/reference/af-only-gnomad.hg38.vcf.gz"
params.outdir      = "results"
params.num_gpus    = 4

process FQ2BAM {
    tag        "${sample}_${type}"
    container  'nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1'
    accelerator params.num_gpus, type: 'nvidia.com/gpu'
    publishDir "${params.outdir}/${sample}/bam", mode: 'copy'

    input:
    tuple val(sample), val(type), path(r1), path(r2)
    path ref
    path known_sites

    output:
    tuple val(sample), val(type), path("${sample}_${type}.bam"), path("${sample}_${type}.recal.txt")

    script:
    """
    pbrun fq2bam \\
        --ref ${ref} \\
        --in-fq ${r1} ${r2} "@RG\\tID:${sample}_${type}\\tLB:lib1\\tPL:ILLUMINA\\tSM:${sample}_${type}\\tPU:unit1" \\
        --knownSites ${known_sites} \\
        --out-bam ${sample}_${type}.bam \\
        --out-recal-file ${sample}_${type}.recal.txt \\
        --num-gpus ${params.num_gpus}
    """
}

process MUTECT2 {
    tag        sample
    container  'nvcr.io/nvidia/clara/clara-parabricks:4.3.1-1'
    accelerator params.num_gpus, type: 'nvidia.com/gpu'
    publishDir "${params.outdir}/${sample}/vcf", mode: 'copy'

    input:
    tuple val(sample), path(tumor_bam), path(tumor_recal), path(normal_bam), path(normal_recal)
    path ref
    path pon
    path germline

    output:
    tuple val(sample), path("${sample}.somatic.unfiltered.vcf.gz")

    script:
    """
    pbrun mutectcaller \\
        --ref ${ref} \\
        --tumor-name ${sample}_tumor \\
        --in-tumor-bam ${tumor_bam} \\
        --in-tumor-recal-file ${tumor_recal} \\
        --normal-name ${sample}_normal \\
        --in-normal-bam ${normal_bam} \\
        --in-normal-recal-file ${normal_recal} \\
        --pon ${pon} \\
        --germline-resource ${germline} \\
        --out-vcf ${sample}.somatic.unfiltered.vcf.gz \\
        --num-gpus ${params.num_gpus}
    """
}

process FILTER_MUTECT {
    tag        sample
    container  'broadinstitute/gatk:4.5.0.0'
    publishDir "${params.outdir}/${sample}/vcf", mode: 'copy'

    input:
    tuple val(sample), path(vcf)
    path ref

    output:
    tuple val(sample), path("${sample}.somatic.filtered.vcf.gz")

    script:
    """
    gatk FilterMutectCalls \\
        -R ${ref} \\
        -V ${vcf} \\
        -O ${sample}.somatic.filtered.vcf.gz
    """
}

workflow {
    // Create paired channel: (sample, type, r1, r2)
    samples_ch = Channel.of(
        tuple(params.sample_id, 'tumor',  file(params.tumor_r1),  file(params.tumor_r2)),
        tuple(params.sample_id, 'normal', file(params.normal_r1), file(params.normal_r2))
    )

    bams = FQ2BAM(samples_ch, file(params.ref), file(params.known_sites))

    // Pair tumor and normal BAMs by sample id
    paired = bams
        .branch {
            tumor:  it[1] == 'tumor'
            normal: it[1] == 'normal'
        }

    pair_ch = paired.tumor
        .map  { s, t, bam, recal -> tuple(s, bam, recal) }
        .join(paired.normal.map { s, t, bam, recal -> tuple(s, bam, recal) })

    raw_vcf = MUTECT2(pair_ch, file(params.ref), file(params.pon), file(params.germline))
    FILTER_MUTECT(raw_vcf, file(params.ref))
}
