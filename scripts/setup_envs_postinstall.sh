#!/usr/bin/env bash
# =============================================================================
# setup_envs_postinstall.sh
# Pós-instalação manual após: mamba env create -f envs/docking.yml
#
# Por que este script existe:
#   meeko depende de openbabel Python bindings, mas o conda já instalou
#   openbabel=3.1.1 com bindings. O pip tenta recompilar via SWIG e falha
#   porque aponta para /usr/local/include em vez de $CONDA_PREFIX/include.
#   Solução: instalar meeko --no-deps (conda já provê todas as deps reais).
# =============================================================================

set -euo pipefail

echo "=== Pós-instalação: ambiente docking ==="

# Verificar que o ambiente existe
if ! mamba env list | grep -q '^docking '; then
    echo "ERRO: ambiente 'docking' não encontrado. Rode primeiro:"
    echo "  mamba env create -f envs/docking.yml --yes"
    exit 1
fi

# meeko: instalar sem recompilar openbabel (já presente via conda)
echo "[1/3] Instalando meeko (--no-deps)..."
mamba run -n docking pip install "meeko>=0.5.0" --no-deps --quiet
echo "      meeko OK"

# HADDOCK3: link para instalação manual em ~/gromacs/haddock3
HADDOCK_DIR="${HOME}/gromacs/haddock3"
if [ -d "$HADDOCK_DIR" ]; then
    echo "[2/3] Linkando HADDOCK3 no ambiente docking..."
    mamba run -n docking pip install -e "$HADDOCK_DIR" --quiet
    echo "      HADDOCK3 OK (instalado de $HADDOCK_DIR)"
else
    echo "[2/3] AVISO: $HADDOCK_DIR não encontrado — instale HADDOCK3 manualmente:"
    echo "      git clone https://github.com/haddocking/haddock3 ~/gromacs/haddock3"
    echo "      mamba run -n docking pip install -e ~/gromacs/haddock3"
fi

# AlphaFold3: link para instalação manual em ~/gromacs/alphafold3
AF3_DIR="${HOME}/gromacs/alphafold3"
if [ -d "$AF3_DIR" ]; then
    echo "[3/3] Linkando AlphaFold3 no ambiente structure..."
    mamba run -n structure pip install -e "$AF3_DIR" --quiet
    echo "      AlphaFold3 OK (instalado de $AF3_DIR)"
else
    echo "[3/3] AVISO: $AF3_DIR não encontrado — instale AF3 manualmente:"
    echo "      git clone https://github.com/google-deepmind/alphafold3 ~/gromacs/alphafold3"
    echo "      mamba run -n structure pip install -e ~/gromacs/alphafold3"
fi

echo ""
echo "=== Verificação rápida ==="
echo -n "docking  → openbabel: "; mamba run -n docking python -c "import openbabel; print('OK')" 2>/dev/null || echo "FALHOU"
echo -n "docking  → meeko:     "; mamba run -n docking python -c "import meeko;    print('OK')" 2>/dev/null || echo "FALHOU"
echo -n "docking  → vina:      "; mamba run -n docking python -c "from vina import Vina; print('OK')" 2>/dev/null || echo "FALHOU"
echo -n "structure → alphafold: "; mamba run -n structure python -c "import alphafold3; print('OK')" 2>/dev/null || echo "FALHOU (normal se pesos não baixados)"
echo -n "structure → jax GPU:   "; mamba run -n structure python -c "import jax; d=jax.devices(); print('OK -', d)" 2>/dev/null || echo "FALHOU"

echo ""
echo "=== Todos os passos concluídos ==="
