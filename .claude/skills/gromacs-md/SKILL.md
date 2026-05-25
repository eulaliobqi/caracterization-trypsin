---
name: gromacs-md
description: Executa simulações de dinâmica molecular (100 ns) de complexos tripsina-inibidor com GROMACS 2025. Trigger: "simulação MD", "GROMACS", "dinâmica molecular", "RMSD RMSF".
---

# GROMACS MD Simulation (Fase 10)

## Quando usar
- Fase 10: após docking selecionar top complexos tripsina-inibidor
- Requer GPU (CUDA 12.4+) e GROMACS 2025.1

## Inputs
- `results/09_docking/top_complexes/*.pdb` (complexos tripsina-inibidor do docking)
- Force field: AMBER99SB-ILDN (incluído no GROMACS)
- Topologia dos inibidores: gerada com ACPYPE/parmchk2

## Pipeline MD (por complexo)

### 1. Preparação topologia
```bash
# Separar proteína e inibidor
grep "^ATOM" complex.pdb > protein.pdb
grep "^HETATM" complex.pdb > ligand.pdb

# Topologia com pdb2gmx (AMBER99SB-ILDN)
gmx pdb2gmx -f protein.pdb -o protein_proc.gro \
    -water tip3p -ff amber99sb-ildn -ignh

# Topologia do inibidor (se peptídeo, incluir diretamente)
# Se pequena molécula: usar ACPYPE para gerar parâmetros GAFF2
acpype -i ligand.pdb -c bcc -n 0

# Combinar topologias
# Editar topol.top para incluir [ligand.itp]
```

### 2. Caixa de simulação e solvente
```bash
# Caixa dodecaédrica, mínimo 1.2 nm de proteína à borda
gmx editconf -f complex.gro -o box.gro -c -d 1.2 -bt dodecahedron

# Solvatação com TIP3P
gmx solvate -cp box.gro -cs spc216.gro -o solvated.gro -p topol.top

# Adicionar íons (neutralizar + 150 mM NaCl)
gmx grompp -f ions.mdp -c solvated.gro -p topol.top -o ions.tpr
gmx genion -s ions.tpr -o ionized.gro -p topol.top \
    -pname NA -nname CL -neutral -conc 0.15
```

### 3. Minimização de energia
```bash
gmx grompp -f em.mdp -c ionized.gro -p topol.top -o em.tpr
gmx mdrun -v -deffnm em -ntmpi 1 -ntomp ${NPROCS:-8}
# Critério: Fmax < 1000 kJ/mol/nm
```

### 4. Equilibração NVT (300 K, 100 ps)
```bash
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -deffnm nvt -ntmpi 1 -ntomp ${NPROCS:-8} -gpu_id 0
```

### 5. Equilibração NPT (1 bar, 100 ps)
```bash
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -deffnm npt -ntmpi 1 -ntomp ${NPROCS:-8} -gpu_id 0
```

### 6. Produção 100 ns
```bash
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -deffnm md -ntmpi 1 -ntomp ${NPROCS:-8} -gpu_id 0
# ~24-72h dependendo do hardware
```

### 7. Análises
```bash
# RMSD backbone (referência: estrutura inicial)
gmx rms -s md.tpr -f md.xtc -o rmsd.xvg -tu ns << EOF
Backbone
Backbone
EOF

# RMSF por resíduo
gmx rmsf -s md.tpr -f md.xtc -o rmsf.xvg -res

# Raio de giração
gmx gyrate -s md.tpr -f md.xtc -o gyrate.xvg

# Contatos H-bond interface
gmx hbond -s md.tpr -f md.xtc -num hbnum.xvg << EOF
Protein
MOL  # grupo do inibidor
EOF
```

## Parâmetros .mdp essenciais (produção)
```ini
integrator    = md
nsteps        = 50000000  ; 100 ns (dt=0.002 ps)
dt            = 0.002
nstxout-compressed = 5000  ; salvar a cada 10 ps
constraints   = all-bonds
constraint_algorithm = LINCS
tcoupl        = V-rescale
tau_t         = 0.1
ref_t         = 300
pcoupl        = Parrinello-Rahman
tau_p         = 2.0
ref_p         = 1.0
coulombtype   = PME
rcoulomb      = 1.2
rvdw          = 1.2
```

## Outputs
- `results/10_md/{complex_id}/md.xtc` (trajetória)
- `results/10_md/{complex_id}/rmsd.xvg`
- `results/10_md/{complex_id}/rmsf.xvg`
- `results/10_md/{complex_id}/analysis_summary.tsv`
