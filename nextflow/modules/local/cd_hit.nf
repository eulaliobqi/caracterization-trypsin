// modules/local/cd_hit.nf — CD-HIT-EST para redução de redundância (95% identity)
// Aplicado ao assembly Trinity antes de TransDecoder

process CD_HIT {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/01_qc/cdhit", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/assembly_qc.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_nr.fasta"), emit: fasta
    tuple val(meta), path("${meta.id}_nr.fasta.clstr"), emit: clusters
    path "versions.yml", emit: versions

    script:
    def identity = params.cdhit_identity ?: 0.95
    """
    cd-hit-est \\
        -i ${fasta} \\
        -o ${meta.id}_nr.fasta \\
        -c ${identity} \\
        -n 10 \\
        -T ${task.cpus} \\
        -M ${(task.memory.toMega() * 0.9).toInteger()} \\
        -d 0

    # Log estatísticas
    echo "Before CD-HIT: \$(grep -c '>' ${fasta}) sequences" > cdhit_stats.txt
    echo "After CD-HIT:  \$(grep -c '>' ${meta.id}_nr.fasta) sequences" >> cdhit_stats.txt
    cat cdhit_stats.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cd-hit: \$(cd-hit-est 2>&1 | head -1 | awk '{print \$4}')
    END_VERSIONS
    """

    stub:
    """
    cp ${fasta} ${meta.id}_nr.fasta
    touch ${meta.id}_nr.fasta.clstr
    touch versions.yml
    """
}
