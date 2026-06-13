// modules/local/diamond.nf — DIAMOND BLASTp
// Processo 1: DIAMOND para hints de ORF (Fase 2)
// Processo 2: DIAMOND_TRYPSIN para identificação de candidatos (Fase 3 — somente IDs)
// Processo 3: INTERSECT_TRYPSINS — interseção DIAMOND ∩ HMMER → FASTAs finais (Fase 3)

process DIAMOND {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/02_orfs/diamond", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/orf_prediction.yml"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_diamond_hints.tsv"), emit: tsv
    path "versions.yml",                                    emit: versions

    script:
    """
    diamond blastp \\
        --query ${fasta} \\
        --db ${params.uniprot_db} \\
        --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs 1 \\
        --sensitive \\
        --threads ${task.cpus} \\
        --out ${meta.id}_diamond_hints.tsv

    echo "DIAMOND hits: \$(wc -l < ${meta.id}_diamond_hints.tsv)"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch ${meta.id}_diamond_hints.tsv
    touch versions.yml
    """
}

// ── FASE 3: Passo 1 — apenas DIAMOND, somente lista de IDs candidatos ─────────
process DIAMOND_TRYPSIN {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/03_trypsin_ids/diamond", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    input:
    tuple val(meta), path(pep)

    output:
    tuple val(meta), path("diamond_trypsin_ids.txt"), emit: trypsin_ids
    tuple val(meta), path("diamond_all.tsv"),         emit: diamond_tsv
    path "versions.yml",                               emit: versions

    script:
    """
    # 1. DIAMOND BLASTp contra UniProt Swiss-Prot
    diamond blastp \\
        --query ${pep} \\
        --db ${params.uniprot_db} \\
        --outfmt 6 qseqid sseqid stitle pident length qlen slen evalue bitscore \\
        --evalue ${params.diamond_evalue} \\
        --max-target-seqs 5 \\
        --sensitive \\
        --threads ${task.cpus} \\
        --out diamond_all.tsv

    # 2. Filtrar apenas sequências com "trypsin" no título do hit (case-insensitive)
    awk -F'\\t' 'tolower(\$3) ~ /trypsin/ {print \$1}' diamond_all.tsv \\
        | sort -u > diamond_trypsin_ids.txt

    N=\$(wc -l < diamond_trypsin_ids.txt)
    echo "DIAMOND candidatos com 'trypsin': \${N}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(diamond --version 2>&1 | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    echo "TRINITY_DN1_c0_g1_i1" > diamond_trypsin_ids.txt
    touch diamond_all.tsv
    touch versions.yml
    """
}

// ── FASE 3: Passo 3 — Interseção real DIAMOND ∩ HMMER (PF00089) ──────────────
// Recebe: pep + diamond_trypsin_ids.txt + trypsin_hmm_ids.txt
// Produz: trypsins_confident.fasta (ambos) + trypsins_suggestive.fasta (um só)
process INTERSECT_TRYPSINS {
    tag        "${meta.id}"
    label      'process_low'
    publishDir "${params.outdir}/03_trypsin_ids", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/annotation.yml"

    input:
    tuple val(meta), path(pep), path(diamond_ids), path(hmm_ids)

    output:
    tuple val(meta), path("trypsins_confident.fasta"),   emit: confident_fasta
    tuple val(meta), path("trypsins_suggestive.fasta"),  emit: suggestive_fasta
    tuple val(meta), path("identification_report.tsv"),  emit: report_tsv
    tuple val(meta), path("trypsin_ids_confident.txt"),  emit: trypsin_ids
    path "versions.yml",                                  emit: versions

    script:
    """
    # Interseção DIAMOND ∩ HMMER via script Python
    # confident = presente nos DOIS métodos (alta especificidade)
    # suggestive = presente em apenas UM método (requer revisão manual)
    python3 ${projectDir}/bin/filter_complete_trypsins.py \\
        --diamond ${diamond_ids} \\
        --hmm_ids ${hmm_ids} \\
        --pep ${pep} \\
        --output_confident trypsins_confident.fasta \\
        --output_suggestive trypsins_suggestive.fasta \\
        --report identification_report.tsv \\
        --ids_confident trypsin_ids_confident.txt

    N_CONF=\$(grep -c '>' trypsins_confident.fasta  2>/dev/null || echo 0)
    N_SUGG=\$(grep -c '>' trypsins_suggestive.fasta 2>/dev/null || echo 0)
    echo "=== Interseção DIAMOND ∩ HMMER ==="
    echo "  Confident (ambos):   \${N_CONF}"
    echo "  Suggestive (um só):  \${N_SUGG}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo ">TRINITY_DN1_c0_g1_i1_confident_stub" > trypsins_confident.fasta
    printf 'MAALGAVLLLCVLPALAARRGIPYSDGICSEVMPVKGRGKTFLVDDLCRSVAFLCGASITDVPDAMTQQIKAGKGLDEALITQNPYEGPVSEAQDLLQKLFGNASVNRFPSQARLSSQPLFIYVAGKLSSGNCATGKPIEVIDFRQK\\n' >> trypsins_confident.fasta
    touch trypsins_suggestive.fasta
    printf 'id\tdiamond_hit\thmm_hit\tconfident\n' > identification_report.tsv
    printf 'TRINITY_DN1_c0_g1_i1\ttrue\ttrue\ttrue\n' >> identification_report.tsv
    echo "TRINITY_DN1_c0_g1_i1" > trypsin_ids_confident.txt
    touch versions.yml
    """
}
