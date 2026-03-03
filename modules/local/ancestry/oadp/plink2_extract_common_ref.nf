process PLINK2_EXTRACT_COMMON_REF {
    // labels are defined in conf/modules.config
    label 'process_low'
    label "${ params.copy_genomes ? 'copy_genomes' : '' }"
    label "plink2" // controls conda, docker, + singularity options

    tag "${meta.target_id ?: meta.id}"

    cachedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir
    storeDir cachedir / "ancestry" / "bed"

    conda "${task.ext.conda}"

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
    tuple val(meta), path(ref_geno), path(ref_pheno), path(ref_variants), path(target_variants)

    output:
    tuple val(meta), path("${output}.bed"), emit: geno
    tuple val(meta), path("${output}.bim"), emit: variants
    tuple val(meta), path("${output}.fam"), emit: pheno
    path "versions.yml", emit: versions

    script:
    def mem_mb = task.memory.toMega()
    output = "${ref_geno.baseName}_common"

    """
    awk '{print \$2}' $target_variants > target_variant_ids.txt

    plink2 \
        --threads $task.cpus \
        --memory $mem_mb \
        --seed 31 \
        --bed $ref_geno \
        --fam $ref_pheno \
        --bim $ref_variants \
        --extract target_variant_ids.txt \
        --make-bed \
        --out $output

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        plink2: \$(plink2 --version 2>&1 | sed 's/^PLINK v//; s/ 64.*\$//' )
    END_VERSIONS
    """
}
