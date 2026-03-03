include { MATCH_VARIANTS } from '../../modules/local/match_variants'
include { MATCH_COMBINE  } from '../../modules/local/match_combine'

workflow MATCH {
    take:
    geno
    pheno
    variants
    scorefile
    ch_intersection

    main:
    ch_versions = Channel.empty()

    // convert to a nested list of length 1 [[scorefile1_path, scorefile2_path]]
    scorefiles = scorefile.collect { [it] }

    variants
        .combine(scorefiles)
        .dump(tag: 'match_variants_input', pretty: true)
        .set { ch_variants }

    MATCH_VARIANTS ( ch_variants )
    ch_versions = ch_versions.mix(MATCH_VARIANTS.out.versions.first())

    // build MATCH_COMBINE inputs per sampleset ID (supports multi-sampleset runs)
    MATCH_VARIANTS.out.matches
        .map { meta, match_path -> tuple(meta.subMap('id'), meta, match_path) }
        .groupTuple()
        .map { id_meta, metas, matches ->
            def split = metas.first().chrom != 'ALL'
            tuple(id_meta, [id: id_meta.id, split: split], matches)
        }
        .set { ch_matches_grouped }

    ch_intersection
        .map { meta, shared_path -> tuple(meta.subMap('id'), shared_path) }
        .groupTuple()
        .set { ch_intersection_grouped }

    ch_matches_grouped
        .join(ch_intersection_grouped, by: 0)
        .combine(scorefiles)
        .map { id_meta, combine_meta, matches, intersections, all_scorefiles ->
            tuple(combine_meta, matches, all_scorefiles.flatten(), intersections)
        }
        .dump(tag: 'match_combine_input', pretty: true)
        .set { ch_match_combine_input }

     MATCH_COMBINE ( ch_match_combine_input )
     ch_versions = ch_versions.mix(MATCH_COMBINE.out.versions)

    emit:
    scorefiles = MATCH_COMBINE.out.scorefile
    db         = MATCH_COMBINE.out.summary
    versions   = ch_versions
}
