/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Define inputs
println "Running with the following parameters:"
println "Input samplesheet: ${params.input}"

no_file_name = file(params.no_file).name
haplotag_qc_script = file("${workflow.projectDir}/subworkflows/local/haplotag_qc/haplotagqc.py")

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
    // whatshap_inputs = inputs
    // .filter { 
    //     ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index, phased_sv_vcf ->
    //     phased_sv_vcf.name == no_file_name
    // }
    // .map { ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index, phased_sv_vcf ->
    //         [ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index]
    //     }
    
    longphase_inputs = inputs.map { ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_snv_vcf_index, phased_sv_vcf ->
            [ref, ref_index, sample, bam, bam_index, phased_snv_vcf, phased_sv_vcf]
        }


    // whatshap_haplotag_results = whatshap_haplotag(whatshap_inputs)
    longphase_haplotag_results = longphase_haplotag(longphase_inputs)
    // indexed_bams = index_bam(longphase_haplotag_results.concat( whatshap_haplotag_results))
    
    // haplotag_qc_ch = indexed_bams.map { 
    //     sample, bam, bam_index ->  
    //     [sample, bam, bam_index, haplotag_qc_script]
    // }
    // haplotag_qc_results = haplotag_qc(haplotag_qc_ch)

}

process index_bam {
    publishDir params.outdir, mode: 'copy'
    container "quay.io/biocontainers/samtools:1.22--h96c455f_0"

    input:
    tuple val(sample), path(bam)

    output:
    tuple val(sample), path(bam), path(indexed_bam)

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
    tuple val(sample), path (haplotagged_bam)

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
    tuple val(sample), path (haplotagged_bam)

    script:
    haplotagged_prefix = "${sample}_longphase_haplotag"
    haplotagged_bam = "${haplotagged_prefix}.bam"
    phased_sv_vcf = phased_sv_vcf.name != no_file_name ? phased_sv_vcf : null

    """
    echo "PATH: $PATH"
    ls /opt/conda/bin

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
    tuple val(sample), path(haplotagged_bam), path(haplotagged_bam_index), path(script)

    output:
    path "*.txt", emit: phasing_qc_txt
    path "*.png", emit: phasing_qc_png


    script:
    """
    python haplotagqc.py \\
        --bam "$haplotagged_bam" \\
        --sample "${sample}"
    """
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/