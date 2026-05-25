// modules/local/alphafold3.nf — AlphaFold3 para predição estrutural
// Requer GPU e modelo AF3 baixado separadamente
// Ver: https://github.com/google-deepmind/alphafold3

process ALPHAFOLD3 {
    tag        "${meta.id}"
    label      'process_gpu'
    publishDir "${params.outdir}/07_structures/alphafold3/${meta.id}", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/structure.yml"

    errorStrategy 'retry'
    maxRetries    2

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("${meta.id}_model_0.pdb"),        emit: pdb
    tuple val(meta), path("${meta.id}_confidence.json"),     emit: confidence
    tuple val(meta), path("${meta.id}.log"),                 emit: log
    path "versions.yml",                                     emit: versions

    script:
    def models_dir  = params.alphafold3_models
    def db_dir      = params.alphafold3_db
    def n_recycles  = params.af3_num_recycles ?: 10
    """
    # Converter FASTA para JSON (formato AF3)
    python3 ${projectDir}/bin/parse_alphafold3.py \\
        --fasta ${fasta} \\
        --id ${meta.id} \\
        --output ${meta.id}_input.json

    # Rodar AlphaFold3
    python3 -m alphafold3.run \\
        --json_path ${meta.id}_input.json \\
        --output_dir . \\
        --model_dir ${models_dir} \\
        --db_dir ${db_dir} \\
        --num_recycling_iters ${n_recycles} \\
        --jackhmmer_n_cpu ${task.cpus} \\
        2>&1 | tee ${meta.id}.log

    # Converter CIF → PDB se necessário
    python3 - <<'PYEOF'
import os
from pathlib import Path
for cif in Path('.').glob('*_model_0.cif'):
    try:
        from Bio.PDB import MMCIFParser, PDBIO
        p = MMCIFParser(QUIET=True).get_structure('m', str(cif))
        io = PDBIO()
        io.set_structure(p)
        io.save(str(cif).replace('.cif', '.pdb'))
        print(f'Converted {cif} to PDB')
    except Exception as e:
        print(f'WARNING: CIF→PDB conversion failed: {e}')
PYEOF

    # Verificar pLDDT
    python3 -c "
import json
conf = json.load(open('${meta.id}_confidence.json'))
if isinstance(conf, dict) and 'atom_plddts' in conf:
    plddt = sum(conf['atom_plddts']) / len(conf['atom_plddts'])
    print(f'pLDDT médio ${meta.id}: {plddt:.1f}')
    with open('${meta.id}_plddt.txt', 'w') as f:
        f.write(f'{plddt:.2f}')
"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        alphafold3: \$(python3 -c "import alphafold3; print(alphafold3.__version__)" 2>/dev/null || echo "3.0")
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    # Stub: criar arquivos vazios para teste sem GPU
    touch ${meta.id}_model_0.pdb
    echo '{"atom_plddts": [85.0, 87.0, 83.0], "ptm": 0.85, "ranking_score": 0.85}' \
        > ${meta.id}_confidence.json
    echo "Stub AlphaFold3 run for ${meta.id}" > ${meta.id}.log
    touch versions.yml
    """
}
