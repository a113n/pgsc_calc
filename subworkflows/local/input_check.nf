//
// Check input samplesheet and get read channels
//

include { FORMAT_SCOREFILES } from '../../modules/local/format_scorefiles'


workflow INPUT_CHECK {
    take:
    input_path // file: /path/to/samplesheet.csv (can be a path string or channel)
    format // csv or JSON
    scorefile // flat list of paths
    chain_files

    main:
    /* all genomic data should be represented as a list of : [[meta], file]

       meta hashmap structure:
        id: experiment label, possibly shared across split genomic files
        is_vcf: boolean, is in variant call format
        is_bfile: boolean, is in PLINK1 fileset format
        is_pfile: boolean, is in PLINK2 fileset format
        chrom: The chromosome associated with the file. If multiple chroms, null.
        n_chrom: Total separate chromosome files per experiment ID
     */

    ch_versions = Channel.empty()
    parsed_input = Channel.empty()

    input = input_path

    if (format.equals("csv")) {
        // Parse and validate each samplesheet independently
        input_path
            .map { csv_file ->
                tuple([samplesheet_path: csv_file, n_chrom: csv_file.countLines() - 1], csv_file)
            }
            .splitCsv(header: true, elem: 1)
            .groupTuple()
            .map { meta, rows ->
                def p = new SamplesheetParser(meta.samplesheet_path, meta.n_chrom as Integer, params.target_build)
                p.verifySamplesheet(rows).collect { row ->
                    def parsed = p.parseCSVRow(row)
                    return parsed
                }
            }
            .flatMap { it }
            .set { parsed_input }
    } else if (format.equals("json")) {
        // Parse and validate each JSON samplesheet independently
        input_path
            .map { json_file ->
                tuple([samplesheet_path: json_file, n_chrom: json_file.countJson()], json_file)
            }
            .splitJson(elem: 1)
            .groupTuple()
            .map { meta, jsonarray ->
                def p = new SamplesheetParser(meta.samplesheet_path, meta.n_chrom as Integer, params.target_build)
                p.verifySamplesheet(jsonarray).collect { json ->
                    def parsed = p.parseJSON(json)
                    return parsed
                }
            }
            .flatMap { it }
            .set { parsed_input }
    }

    parsed_input.branch {
                vcf: it[0].is_vcf
                bfile: it[0].is_bfile
                pfile: it[0].is_pfile
        }
        .set { ch_branched }

    // branch is like a switch statement, so only one bed / bim was being
    // returned
    ch_branched.bfile.multiMap { it ->
        bed: [it[0], it[1][0]]
        bim: [it[0], it[1][1]]
        fam: [it[0], it[1][2]]
    }
        .set { ch_bfiles }

    ch_branched.pfile.multiMap { it ->
        pgen: [it[0], it[1][0]]
        pvar: [it[0], it[1][1]]
        psam: [it[0], it[1][2]]
    }
        .set { ch_pfiles }

    FORMAT_SCOREFILES ( scorefile, chain_files )

    versions = ch_versions.mix(FORMAT_SCOREFILES.out.versions)

    ch_bfiles.bed.mix(ch_pfiles.pgen).dump(tag: 'geno').set { geno }
    ch_bfiles.bim.mix(ch_pfiles.pvar).dump(tag: 'variants').set { variants }
    ch_bfiles.fam.mix(ch_pfiles.psam).dump(tag: 'pheno').set { pheno }
    ch_branched.vcf.dump(tag: 'input').set{vcf}
    FORMAT_SCOREFILES.out.scorefiles.dump(tag: 'input').set{ scorefiles }
    FORMAT_SCOREFILES.out.log_scorefiles.dump(tag: 'input').set{ log_scorefiles }

    emit:
    geno
    variants
    pheno
    vcf
    scorefiles
    log_scorefiles
    versions
}

