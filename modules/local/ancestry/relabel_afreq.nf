process RELABEL_AFREQ {
    // labels are defined in conf/modules.config
    label 'process_medium'
    label 'pgscatalog_utils' // controls conda, docker, + singularity options

    tag "$meta.id $meta.effect_type $target_format"

    basedir = params.genotypes_cache ? file(params.genotypes_cache) : workDir
    storeDir basedir / "ancestry" / "relabel" / "afreq"

    conda "${task.ext.conda}"

    container "${ workflow.containerEngine == 'singularity' &&
        !task.ext.singularity_pull_docker_container ?
        "${task.ext.singularity}${task.ext.singularity_version}" :
        "${task.ext.docker}${task.ext.docker_version}" }"

    input:
    tuple val(meta), path(target), path(matched)

    output:
    tuple val(relabel_meta), path("${output}"), emit: relabelled
    path "versions.yml", emit: versions

    script:
    target_format = target.getName().tokenize('.')[1] // test.tar.gz -> tar, test.var -> var
    relabel_meta = meta.plus(['target_format': target_format]) // .plus() returns a new map
    output_mode = "--split --combined" // always output split and combined data to make life easier
    col_from = "ID_TARGET"
    col_to = "ID_REF"
    output = "${meta.id}.${target_format}*"
    target_for_relabel = target

    if (target_format == "afreq") {
        col_from = "ID_REF"
        col_to = "ID_TARGET"
        output_mode = "--combined"
        target_for_relabel = "${meta.id}.filtered.afreq"
    }
    """
    if [[ "$target_format" == "afreq" ]]; then
        # Keep only reference AFREQ rows that are actually present in the map(s)
        # so pgscatalog-relabel doesn't fail on unmapped IDs.
        python - <<'PY'
import gzip
from pathlib import Path

matched = Path("$matched")
target = Path("$target")
out = Path("$target_for_relabel")

ids = set()
with gzip.open(matched, "rt") as f:
    header = next(f, None)
    for line in f:
        cols = line.rstrip("\\n").split("\t")
        if len(cols) > 1:
            ids.add(cols[1])  # ID_REF

def open_target_text(path: Path):
    if path.suffix == ".zst":
        import zstandard as zstd
        return zstd.open(str(path), mode="rt")
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return open(path, "rt")

with open_target_text(target) as fin, open(out, "wt") as fout:
    header = next(fin, None)
    if header is not None:
        fout.write(header)
    for line in fin:
        cols = line.rstrip("\\n").split("\t")
        if len(cols) > 1 and cols[1] in ids:
            fout.write(line)
PY
    fi

    pgscatalog-relabel --maps $matched \
        --col_from $col_from \
        --col_to $col_to \
        --target_file $target_for_relabel \
        --target_col ID \
        --dataset ${meta.id}.${target_format} \
        --verbose \
        $output_mode \
        --outdir \$PWD

    cat <<-END_VERSIONS > versions.yml
    ${task.process.tokenize(':').last()}:
        pgscatalog.core: \$(echo \$(python -c 'import pgscatalog.core; print(pgscatalog.core.__version__)'))
    END_VERSIONS
    """
}
