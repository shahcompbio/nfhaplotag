/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { MULTIQC } from '../modules/nf-core/multiqc'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define inputs
println "Running with the following parameters:"
println "Input samplesheet: ${params.input}"
println "Output directory: ${params.outdir}"
println "skip_longphase: ${params.skip_longphase}"
println "skip_whatshap: ${params.skip_whatshap}"
println "skip_qc: ${params.skip_qc}"

no_file_name = file(params.no_file).name
haplotag_qc_script = file("${workflow.projectDir}/subworkflows/local/haplotag_qc/haplotagqc.py")
any_haplotagging = !params.skip_whatshap || !params.skip_longphase

workflow NFPHAPLOTAG {

    take:
    samplesheet // channel: samplesheet read in from --input
    main:

    inputs = Channel
         .fromPath(samplesheet)
         .splitCsv(header: true)
         .map {
            row -> [
                file(row.ref),
                file(row.ref_index),
                row.sample,
                file(row.bam),
                file(row.bam_index),
                file(row.phased_snv_vcf),
                file(row.phased_snv_vcf_index),
                row.phased_sv_vcf ? file(row.phased_sv_vcf) : file(params.no_file) // use no_file as sentinel if sv is not provided
            ]
            }

    // whatshap
    if (params.skip_whatshap) {
        whatshap_haplotag_results = Channel.empty()
    } else {
        whatshap_inputs = inputs
            .filter {
                ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index, phased_sv_vcf ->
                phased_sv_vcf.name == no_file_name
            }
            .map { ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index, phased_sv_vcf ->
                    [ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index]
                }

        whatshap_haplotag_results = whatshap_haplotag(whatshap_inputs)
    }

    // longphase
    if (params.skip_longphase) {
        longphase_haplotag_results = Channel.empty()
    } else {
    longphase_inputs = inputs.map {
        ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index, phased_sv_vcf ->
            [ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_sv_vcf]
        }

    longphase_haplotag_results = longphase_haplotag(longphase_inputs)
    }

    if (any_haplotagging) {
        indexed_bams = index_bam(longphase_haplotag_results.concat( whatshap_haplotag_results))
    }
    if (!params.skip_qc && any_haplotagging) {
        haplotag_qc_ch = indexed_bams.map {
            sample, tool, bam, bam_index ->
            [sample, tool, bam, bam_index, haplotag_qc_script]
        }
        haplotag_qc_results = haplotag_qc(haplotag_qc_ch)

        multiqc_results = MULTIQC(
            haplotag_qc_results.flatten().map{it -> file(it)}.collect(),
            [], // no multiqc_config provided
            [], // no extra_multiqc_config provided
            [], // no multiqc_logo provided
            [], // no replace_names provided
            [] // no sample_names provided
        )
    }

}

process index_bam {
    publishDir params.outdir, mode: 'copy'
    container "quay.io/biocontainers/samtools:1.22--h96c455f_0"

    input:
    tuple val(sample), val(tool), path(bam)

    output:
    tuple val(sample), val(tool), path(bam), path(indexed_bam)

    script:
    indexed_bam = "${bam}.bai"
    """
    samtools index "$bam"
    """
}

process whatshap_haplotag {
    publishDir params.outdir, mode: 'copy'
    container "quay.io/biocontainers/whatshap:2.8--py312hf731ba3_0"

    input:
    tuple path(ref), path(ref_index), val(sample), path(bam), path(bam_index), path(phased_vcf), path(phased_vcf_index)

    output:
    tuple val(sample), val("whatshap"), path (haplotagged_bam)

    script:
    haplotagged_bam = "${sample}_whatshap_haplotag.bam"
    """
    whatshap haplotag \\
        "$phased_vcf" \\
        "$bam" \\
        --reference "$ref" \\
        --ignore-read-groups \\
        -o "$haplotagged_bam" \\
        --tag-supplementary
    """
}

process longphase_haplotag {
    publishDir params.outdir, mode: 'copy'
    container "quay.io/biocontainers/longphase:1.7.3--hf5e1c6e_0"

    input:
    tuple path(ref), path(ref_index), val(sample), path(bam), path(bam_index), path(phased_snv_vcf), path(phased_sv_vcf)

    output:
    tuple val(sample), val("longphase"), path (haplotagged_bam)

    script:
    haplotagged_prefix = "${sample}_longphase_haplotag"
    haplotagged_bam = "${haplotagged_prefix}.bam"
    phased_sv_vcf = phased_sv_vcf.name != no_file_name ? phased_sv_vcf : null

    """
    longphase haplotag \\
        -r "$ref" \\
        -s "$phased_snv_vcf" \\
        ${phased_sv_vcf ? "--sv-file $phased_sv_vcf" : ""} \\
        -b "$bam" \\
        -t ${task.cpus} \\
        -o "$haplotagged_prefix"
    """
}

process haplotag_qc {
    publishDir params.outdir, mode: 'copy'
    container "quay.io/shahlab_singularity/haplotagqc:latest"

    input:
    tuple val(sample), val(tool), path(haplotagged_bam), path(haplotagged_bam_index), path(script)

    output:
    tuple path("*.txt"), path("*.png"), path("*.csv")


    script:
    """
    python haplotagqc.py \\
        --bam "$haplotagged_bam" \\
        --sample "$sample" \\
        --tool "$tool"
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
