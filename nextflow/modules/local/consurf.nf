process CONSURF {
    tag        "${meta.id}"
    label      'process_medium'
    publishDir "${params.outdir}/08_validation/consurf", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/structure.yml"

    input:
    tuple val(meta), path(pdb)

    output:
    tuple val(meta), path("${meta.id}_consurf_grades.txt"), emit: grades
    path "versions.yml",                                     emit: versions

    script:
    """
    # ConSurf para mapear conservação evolutiva na estrutura
    # Usar alinhamento já construído na Fase 6
    ALIGNMENT="${params.outdir}/06_phylogeny/mafft/*.fasta"

    if command -v consurf &>/dev/null && ls \$ALIGNMENT 2>/dev/null; then
        consurf \\
            -pdb ${pdb} \\
            -MSA \$(ls \$ALIGNMENT | head -1) \\
            -res ${meta.id}_consurf_grades.txt \\
            -method Bayesian
    else
        # Placeholder se ConSurf não disponível
        python3 - <<'PYEOF'
from Bio.PDB import PDBParser
p = PDBParser(QUIET=True)
struct = p.get_structure("${meta.id}", "${pdb}")
model = struct[0]
with open("${meta.id}_consurf_grades.txt", "w") as f:
    f.write("# ConSurf grades placeholder — requires ConSurf installation\\n")
    f.write("POS\tSEQ\tCONSURF_GRADE\tCOLOR\\n")
    pos = 1
    for chain in model:
        for res in chain:
            if res.id[0] == ' ':
                resname = res.resname[:3]
                f.write(f"{pos}\t{resname}\t5\tblue\\n")  # conservação média
                pos += 1
PYEOF
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | awk '{print \$2}')
    END_VERSIONS
    """

    stub:
    """
    echo "# ConSurf grades stub" > ${meta.id}_consurf_grades.txt
    touch versions.yml
    """
}
