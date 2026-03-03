//
// Do a PCA on reference data and project target genomes into the PCA space
//

include { EXTRACT_DATABASE } from '../../../modules/local/ancestry/extract_database'
include { INTERSECT_VARIANTS } from '../../../modules/local/ancestry/intersect_variants'
include { FILTER_VARIANTS } from '../../../modules/local/ancestry/filter_variants'
include { PLINK2_MAKEBED as PLINK2_MAKEBED_TARGET; PLINK2_MAKEBED as PLINK2_MAKEBED_REF } from '../../../modules/local/ancestry/oadp/plink2_makebed'
include { INTERSECT_THINNED } from '../../../modules/local/ancestry/oadp/intersect_thinned'
include { RELABEL_IDS } from '../../../modules/local/ancestry/relabel_ids'
include { PLINK2_ORIENT } from '../../../modules/local/ancestry/oadp/plink2_orient'
include { PLINK2_EXTRACT_COMMON_REF } from '../../../modules/local/ancestry/oadp/plink2_extract_common_ref'
include { FRAPOSA_PCA } from '../../../modules/local/ancestry/oadp/fraposa_pca'
include { FRAPOSA_PROJECT } from '../../../modules/local/ancestry/oadp/fraposa_project'

workflow ANCESTRY_PROJECT {
    take:
    geno
    pheno
    variants
    vmiss
    afreq
    reference
    target_build

    main:
    ch_versions = Channel.empty()

    // sort order is _very important_
    // input order to modules must always be: geno, pheno, variants, e.g.:
    // .pgen, .psam, .pvar.zst in plink2
    // .bed, .fam, .bim.zst in plink1
    // it's assumed variants are zstd compressed at the start of the workflow
    geno.concat(pheno, variants)
        .groupTuple(size: 3, sort: { it.toString().split("\\.")[-1] } )
        .set { ch_genomes }

    ch_genomes
        .map { it.first().subMap('id') }
        .unique()
        .set { ch_sampleset }

    //
    // STEP 0: extract the reference data once (don't do it inside separate processes)
    //

    EXTRACT_DATABASE( reference )

    EXTRACT_DATABASE.out.grch38
        .concat(EXTRACT_DATABASE.out.grch37)
        .filter { it.first().build == target_build }
        .set { ch_db }

    ch_versions = ch_versions.mix(EXTRACT_DATABASE.out.versions.first())

    ch_db.map {
        def meta = [:].plus(it.first())
        meta.is_pfile = true
        meta.id = 'reference'
        meta.chrom = 'ALL'
        return tuple(meta, it.tail())
    }
        .transpose()
        .branch {
            geno: it.last().getExtension() == 'pgen'
            pheno: it.last().getExtension() == 'psam'
            var: it.last().getExtension() == 'zst'
        }
        .set{ ch_ref }

    //
    // STEP 1: get overlapping variants across reference and target ------------
    //

    ch_genomes
        .join(vmiss, failOnMismatch: true)
        .join(afreq, failOnMismatch: true)
        .combine( ch_db.map{ it.tail() } ) // (drop hashmap)
        .flatten()
        .buffer(size: 9)
        .set { ch_ref_combined }

    INTERSECT_VARIANTS ( ch_ref_combined )
    ch_versions = ch_versions.mix(INTERSECT_VARIANTS.out.versions.first())

    //
    // STEP 2: filter variants in reference and target datasets ----------------
    //

    EXTRACT_DATABASE.out.grch37_king
        .concat(EXTRACT_DATABASE.out.grch38_king)
        .filter { it.first().build == params.target_build }
        .map { it.last() }
        .set { ch_king }

    Channel.of(
        [['build': 'GRCh37'], file(params.ld_grch37, checkIfExists: true)],
        [['build': 'GRCh38'], file(params.ld_grch38, checkIfExists: true)]
    )
        .filter{ it.first().build == params.target_build }
        .map{ it.last() }
        .set{ ch_ld }

    ch_db
        .map { it.tail() }
        .first()
        .set { ch_db_files }

    INTERSECT_VARIANTS.out.intersection
        .map { tuple(it.first().subMap('id'), it.last()) }
        .groupTuple()
        .set { ch_match_reports }

    ch_match_reports
        .combine( ch_db_files )
        .map { row ->
            def id_meta  = row[0]
            def matched  = row[1]
            def ref_geno = row[2]
            def ref_pheno = row[3]
            def ref_var  = row[4]
            tuple(id_meta.plus(['build': params.target_build]), matched, ref_geno, ref_pheno, ref_var)
        }
        .combine( ch_ld.first() )
        .map { row ->
            def meta = row[0]
            def matched = row[1]
            def ref_geno = row[2]
            def ref_pheno = row[3]
            def ref_var = row[4]
            def ld = row[5]
            tuple(meta, matched, ref_geno, ref_pheno, ref_var, ld)
        }
        .combine( ch_king.first() )
        .map { row ->
            def meta = row[0]
            def matched = row[1]
            def ref_geno = row[2]
            def ref_pheno = row[3]
            def ref_var = row[4]
            def ld = row[5]
            def king = row[6]
            tuple(meta, matched, ref_geno, ref_pheno, ref_var, ld, king)
        }
        .set{ ch_filter_input }

    FILTER_VARIANTS ( ch_filter_input )
    ch_versions = ch_versions.mix(FILTER_VARIANTS.out.versions.first())

    FILTER_VARIANTS.out.prune_in
        .map { meta, pruned -> tuple(meta.subMap('id'), pruned) }
        .set { ch_prune_by_id }

    // -------------------------------------------------------------------------
    // ref -> thinned bfile for fraposa
    //

    FILTER_VARIANTS.out.ref
        .map { build_meta, pgen, psam, pvar -> tuple(build_meta.subMap('id'), build_meta, pgen, psam, pvar) }
        .join(ch_prune_by_id, by: 0)
        .map { id_meta, build_meta, pgen, psam, pvar, pruned ->
            def m = [:].plus(build_meta)
            m.target_id = build_meta.id
            m.id = 'reference'
            m.chrom = 'ALL'
            m.is_pfile = true
            tuple(m, pgen, psam, pvar, pruned)
        }
        .set { ch_makebed_ref }

    PLINK2_MAKEBED_REF ( ch_makebed_ref )
    ch_versions = ch_versions.mix(PLINK2_MAKEBED_REF.out.versions.first())

    // -------------------------------------------------------------------------
    // 0. targets -> intersect with thinned reference variants
    // 1. combine split targets into one file
    // 2. relabel
    // 3. then convert to bim

    ch_genomes
        .map { tuple(it.first().subMap('id'), it.first()) }
        .groupTuple()
        .map { id_meta, metas -> tuple(id_meta, metas.first()) }
        .set { ch_geno_meta }

    ch_genomes
        .map { tuple(it.first().subMap('id'), it.last()) }
        .groupTuple()
        .map { id_meta, genomes -> tuple(id_meta, genomes.flatten().unique()) }
        .set { ch_genome_list }

    // [meta, [matches], pruned, geno_meta, [gigantic list of all pfiles]]
    ch_match_reports
        .join( ch_geno_meta, by: 0 )
        .join( ch_genome_list, by: 0 )
        .join( ch_prune_by_id, by: 0 )
        .map { id_meta, matched, geno_meta, genomes, pruned ->
            tuple(id_meta.plus(['build': params.target_build]), matched, pruned, geno_meta, genomes)
        }
        .set{ ch_intersect_thinned_input }

    INTERSECT_THINNED ( ch_intersect_thinned_input )
    //ch_versions = ch_versions.mix(INTERSECT_THINNED.out.versions)

    // [meta, variants, thinned] keyed by sampleset id
    INTERSECT_THINNED.out.match_thinned
        .map { meta, match_path -> tuple(meta.subMap('id'), match_path) }
        .set { ch_match_thinned_by_id }

    INTERSECT_THINNED.out.variants
        .map { meta, variant_path -> tuple(meta.subMap('id'), meta, variant_path) }
        .join(ch_match_thinned_by_id, by: 0)
        .map { id_meta, meta, variant_path, match_path -> tuple(meta, variant_path, match_path) }
        .set { ch_thinned_target }

    RELABEL_IDS( ch_thinned_target )
    ch_versions = ch_versions.mix(RELABEL_IDS.out.versions.first())

    RELABEL_IDS.out.relabelled
        .flatMap { meta, relabelled_paths ->
            relabelled_paths instanceof List ?
                relabelled_paths.collect { p -> tuple(meta, p) } :
                [tuple(meta, relabelled_paths)]
        }
        .filter { meta, relabelled_path -> relabelled_path.getName().contains('ALL') }
        .map { meta, relabelled_path -> tuple(meta.subMap('id'), relabelled_path) }
        .set { ch_ref_relabelled_variants }

    target_extract = file(projectDir / "assets" / "NO_FILE") // optional input for PLINK2_MAKEBED

    // [meta, pgen, psam, relabelled pvar, optional_input]
    INTERSECT_THINNED.out.geno
        .join(INTERSECT_THINNED.out.pheno, by: 0)
        .map { meta, pgen, psam -> tuple(meta.subMap('id'), meta, pgen, psam) }
        .join(ch_ref_relabelled_variants, by: 0)
        .map { id_meta, meta, pgen, psam, relabelled_pvar -> tuple(meta, pgen, psam, relabelled_pvar, target_extract) }
        .set{ ch_target_makebed_input }

    PLINK2_MAKEBED_TARGET ( ch_target_makebed_input )
    ch_versions = ch_versions.mix(PLINK2_MAKEBED_TARGET.out.versions.first())

    // make sure allele order matches across ref / target ----------------------
    // (because plink1 format is very annoying about this sort of thing)

    // [meta, target_bed, target_fam, target_bim, ref_bim]
    PLINK2_MAKEBED_REF.out.variants
        .map { ref_meta, ref_variants -> tuple([id: ref_meta.target_id], ref_variants) }
        .set { ch_ref_bim_by_id }

    PLINK2_MAKEBED_TARGET.out.geno
        .join(PLINK2_MAKEBED_TARGET.out.pheno, by: 0)
        .join(PLINK2_MAKEBED_TARGET.out.variants, by: 0)
        .map { meta, geno, pheno, variants -> tuple(meta.subMap('id'), meta, geno, pheno, variants) }
        .join(ch_ref_bim_by_id, by: 0)
        .map { id_meta, meta, geno, pheno, variants, ref_variants -> tuple(meta, geno, pheno, variants, ref_variants) }
        .set{ ch_orient_input }

    PLINK2_ORIENT( ch_orient_input )
    ch_versions = ch_versions.mix(PLINK2_ORIENT.out.versions.first())

    // fraposa -----------------------------------------------------------------

    PLINK2_MAKEBED_REF.out.geno
        .concat(PLINK2_MAKEBED_REF.out.pheno, PLINK2_MAKEBED_REF.out.variants)
        .groupTuple(size: 3)
        .map { meta, files -> tuple([id: meta.target_id], meta, files[0], files[1], files[2]) }
        .join(PLINK2_ORIENT.out.variants.map { meta, variants -> tuple(meta.subMap('id'), variants) }, by: 0)
        .map { id_meta, ref_meta, ref_geno, ref_pheno, ref_variants, target_variants -> tuple(ref_meta, ref_geno, ref_pheno, ref_variants, target_variants) }
        .set { ch_common_ref_input }

    PLINK2_EXTRACT_COMMON_REF( ch_common_ref_input )
    ch_versions = ch_versions.mix(PLINK2_EXTRACT_COMMON_REF.out.versions.first())

    PLINK2_EXTRACT_COMMON_REF.out.geno
        .concat(PLINK2_EXTRACT_COMMON_REF.out.pheno, PLINK2_EXTRACT_COMMON_REF.out.variants)
        .groupTuple(size: 3)
        .map { meta, files -> tuple(meta, files[0], files[1], files[2]) }
        .set { ch_fraposa_ref }

    PLINK2_ORIENT.out.geno
        .join(PLINK2_ORIENT.out.pheno, by: 0)
        .join(PLINK2_ORIENT.out.variants, by: 0)
        // use full target .fam as split_fam fallback to guarantee FRAPOSA input
        .map { meta, geno, pheno, variants -> tuple(meta.subMap('id'), meta, geno, pheno, variants, pheno) }
        .set { ch_fraposa_target }

    FRAPOSA_PCA( ch_fraposa_ref.map { meta, ref_geno, ref_pheno, ref_variants -> tuple(meta, ref_geno, ref_pheno, ref_variants) }, geno )
    ch_versions = ch_versions.mix(FRAPOSA_PCA.out.versions.first())

    // ... and project split target genomes
    ch_fraposa_ref
        .combine( ch_fraposa_target )
        .filter { row ->
            def ref_geno = row[1]
            def target_meta = row[5]
            def target_id = target_meta.id.toString()
            def ref_base = ref_geno.getBaseName()
            return ref_base.contains("_reference_${target_id}_extracted") || ref_base.endsWith("_reference_extracted")
        }
        .map { row ->
            def ref_meta = row[0]
            def ref_geno = row[1]
            def ref_pheno = row[2]
            def ref_variants = row[3]
            def target_geno = row[6]
            def target_pheno = row[7]
            def target_variants = row[8]
            def split_fam = row[9]
            tuple(ref_meta, ref_geno, ref_pheno, ref_variants, target_geno, target_pheno, target_variants, split_fam)
        }
        .combine(FRAPOSA_PCA.out.pca.map { [it] })
        .map { ref_meta, ref_geno, ref_pheno, ref_variants, target_geno, target_pheno, target_variants, split_fam, pca ->
            tuple(ref_meta, ref_geno, ref_pheno, ref_variants, target_geno, target_pheno, target_variants, split_fam, pca)
        }
        .set { ch_fraposa_input }

    // project targets into reference PCA space
    FRAPOSA_PROJECT( ch_fraposa_input )
    ch_versions = ch_versions.mix(FRAPOSA_PROJECT.out.versions.first())

    // group together ancestry projections for each sampleset
    // different samplesets will have different ancestry projections after intersection
    FRAPOSA_PCA.out.pca
        .flatten()
        .filter { it.getExtension() == 'pcs' }
        .map { [it] }
        .set { ch_ref_projections }

    FRAPOSA_PROJECT.out.pca
        .groupTuple()
        .set { ch_projections }

    emit:
    intersection = INTERSECT_VARIANTS.out.intersection
    intersect_count = INTERSECT_VARIANTS.out.intersect_count.collect()
    projections = ch_projections.combine(ch_ref_projections)
    ref_geno = ch_ref.geno
    ref_pheno = ch_ref.pheno
    ref_var = ch_ref.var
    relatedness = ch_king
    ref_afreq = FILTER_VARIANTS.out.afreq.map { it.last() }.first()
    versions = ch_versions

}
