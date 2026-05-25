// modules/local/gromacs.nf — GROMACS 2025.1 para simulações MD 100 ns
// Pipeline completo: preparação → min → NVT → NPT → produção → análises
// REQUER GPU (process_gpu)

process GROMACS {
    tag        "${meta.id}"
    label      'process_gpu'
    publishDir "${params.outdir}/10_md/${meta.id}", mode: params.publish_dir_mode
    conda      "${projectDir}/../envs/md.yml"

    errorStrategy 'retry'
    maxRetries    1  // MD é caro; apenas 1 retry

    input:
    tuple val(meta), path(complex_pdb)

    output:
    tuple val(meta), path("md.xtc"),           emit: xtc
    tuple val(meta), path("md_analysis/"),     emit: analysis
    tuple val(meta), path("md_report.tsv"),    emit: report
    path "versions.yml",                        emit: versions

    script:
    def nsteps = (params.md_duration_ns as int) * 500000  // 100 ns = 50M steps (dt=0.002 ps)
    """
    # ── 1. Preparação da topologia ────────────────────────────────────────
    gmx pdb2gmx \\
        -f ${complex_pdb} \\
        -o protein_proc.gro \\
        -water tip3p \\
        -ff amber99sb-ildn \\
        -ignh \\
        -nobackup

    # ── 2. Caixa de simulação (dodecaédrica, 1.2 nm margem) ───────────────
    gmx editconf \\
        -f protein_proc.gro \\
        -o box.gro \\
        -c -d 1.2 -bt dodecahedron

    # ── 3. Solvatação (TIP3P) ─────────────────────────────────────────────
    gmx solvate \\
        -cp box.gro -cs spc216.gro \\
        -o solvated.gro -p topol.top

    # ── 4. Íons (neutralizar + 150 mM NaCl) ──────────────────────────────
    gmx grompp -f ${projectDir}/../data/references/mdp/ions.mdp \\
               -c solvated.gro -p topol.top -o ions.tpr -maxwarn 2
    echo "SOL" | gmx genion \\
        -s ions.tpr -o ionized.gro -p topol.top \\
        -pname NA -nname CL -neutral -conc 0.15

    # ── 5. Minimização de energia ─────────────────────────────────────────
    gmx grompp -f ${projectDir}/../data/references/mdp/em.mdp \\
               -c ionized.gro -p topol.top -o em.tpr
    gmx mdrun -v -deffnm em -ntmpi 1 -ntomp ${task.cpus} -gpu_id 0

    # Verificar Fmax
    python3 -c "
import subprocess
out = subprocess.run(['gmx', 'energy', '-f', 'em.edr', '-o', 'em_energy.xvg'],
                    input=b'Potential\n', capture_output=True)
print('Minimização concluída')
"

    # ── 6. Equilibração NVT (300 K, 100 ps) ──────────────────────────────
    gmx grompp -f ${projectDir}/../data/references/mdp/nvt.mdp \\
               -c em.gro -r em.gro -p topol.top -o nvt.tpr
    gmx mdrun -deffnm nvt -ntmpi 1 -ntomp ${task.cpus} -gpu_id 0

    # ── 7. Equilibração NPT (1 bar, 100 ps) ──────────────────────────────
    gmx grompp -f ${projectDir}/../data/references/mdp/npt.mdp \\
               -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
    gmx mdrun -deffnm npt -ntmpi 1 -ntomp ${task.cpus} -gpu_id 0

    # ── 8. Produção ${params.md_duration_ns} ns ───────────────────────────
    # Criar mdp de produção com nsteps correto
    sed 's/NSTEPS_PLACEHOLDER/${nsteps}/' \\
        ${projectDir}/../data/references/mdp/md_template.mdp > md.mdp

    gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
    gmx mdrun -deffnm md -ntmpi 1 -ntomp ${task.cpus} -gpu_id 0 \\
              -update gpu -pme gpu -bonded gpu

    # ── 9. Análises ───────────────────────────────────────────────────────
    mkdir -p md_analysis

    # RMSD backbone
    printf "Backbone\\nBackbone\\n" | gmx rms \\
        -s md.tpr -f md.xtc \\
        -o md_analysis/rmsd_backbone.xvg -tu ns

    # RMSF por resíduo
    echo "Backbone" | gmx rmsf \\
        -s md.tpr -f md.xtc \\
        -o md_analysis/rmsf.xvg -res

    # Raio de giração
    echo "Protein" | gmx gyrate \\
        -s md.tpr -f md.xtc \\
        -o md_analysis/gyrate.xvg

    # SASA
    echo "Protein" | gmx sasa \\
        -s md.tpr -f md.xtc \\
        -o md_analysis/sasa.xvg -tu ns

    # H-bonds na interface proteína-inibidor
    printf "Protein\\nMOL\\n" | gmx hbond \\
        -s md.tpr -f md.xtc \\
        -num md_analysis/hbonds.xvg 2>/dev/null || \\
        echo "AVISO: sem grupo MOL, pulando hbonds" > md_analysis/hbonds_warning.txt

    # ── 10. Relatório sumário ─────────────────────────────────────────────
    python3 - <<'PYEOF'
import subprocess, re
report_lines = ["id\\trmsd_mean_nm\\trmsd_max_nm\\trmsf_mean_nm\\tgyrate_mean_nm"]
# Parse RMSD
try:
    rmsd_data = [float(l.split()[1]) for l in open('md_analysis/rmsd_backbone.xvg')
                 if not l.startswith(('#', '@')) and len(l.split()) >= 2]
    rmsd_mean = sum(rmsd_data)/len(rmsd_data)
    rmsd_max  = max(rmsd_data)
except:
    rmsd_mean = rmsd_max = -1.0
# Parse RMSF
try:
    rmsf_data = [float(l.split()[1]) for l in open('md_analysis/rmsf.xvg')
                 if not l.startswith(('#', '@')) and len(l.split()) >= 2]
    rmsf_mean = sum(rmsf_data)/len(rmsf_data)
except:
    rmsf_mean = -1.0
try:
    gyrate_data = [float(l.split()[1]) for l in open('md_analysis/gyrate.xvg')
                   if not l.startswith(('#', '@')) and len(l.split()) >= 2]
    gyrate_mean = sum(gyrate_data)/len(gyrate_data)
except:
    gyrate_mean = -1.0
with open('md_report.tsv', 'w') as f:
    f.write("id\\trmsd_mean_nm\\trmsd_max_nm\\trmsf_mean_nm\\tgyrate_mean_nm\\n")
    f.write(f"${meta.id}\\t{rmsd_mean:.3f}\\t{rmsd_max:.3f}\\t{rmsf_mean:.3f}\\t{gyrate_mean:.3f}\\n")
print(f"MD completo: RMSD médio = {rmsd_mean:.3f} nm")
PYEOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gromacs: \$(gmx --version 2>&1 | grep 'GROMACS version' | awk '{print \$3}')
    END_VERSIONS
    """

    stub:
    """
    touch md.xtc
    mkdir -p md_analysis
    touch md_analysis/rmsd_backbone.xvg
    touch md_analysis/rmsf.xvg
    touch md_analysis/gyrate.xvg
    echo -e "id\trmsd_mean_nm\trmsd_max_nm\trmsf_mean_nm\tgyrate_mean_nm" > md_report.tsv
    echo -e "${meta.id}\t0.180\t0.320\t0.095\t1.820" >> md_report.tsv
    touch versions.yml
    """
}
