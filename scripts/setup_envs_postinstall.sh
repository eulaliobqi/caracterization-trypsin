#!/usr/bin/env bash
# =============================================================================
# setup_envs_postinstall.sh
# Pós-instalação manual após criar todos os envs conda do pipeline.
#
# Por que este script existe:
#   meeko depende de openbabel Python bindings, mas o conda já instalou
#   openbabel=3.1.1 com bindings. O pip tenta recompilar via SWIG e falha
#   porque aponta para /usr/local/include em vez de $CONDA_PREFIX/include.
#   Solução: instalar meeko --no-deps (conda já provê todas as deps reais).
#
#   AlphaFold3 exige Python >=3.12; structure.yml agora usa python=3.12.
#   Se o env structure existir com Python 3.11, recriar:
#     mamba env remove -n structure && mamba env create -f envs/structure.yml
# =============================================================================

# NÃO usar set -e: alguns steps são opcionais e não devem abortar o script
set -uo pipefail

ERRORS=0

echo "=== Pós-instalação: ambientes do pipeline ==="
echo ""

# ── Detectar base conda ────────────────────────────────────────────────────
CONDA_BASE=$(conda info --base 2>/dev/null || mamba info --base 2>/dev/null || echo "$HOME/miniforge3")
echo "Conda base: $CONDA_BASE"
echo ""

# ── 1. meeko no docking env ───────────────────────────────────────────────
DOCKING_ENV_PATH="${CONDA_BASE}/envs/docking"
if [ ! -d "$DOCKING_ENV_PATH" ]; then
    echo "ERRO: ambiente 'docking' não encontrado em $DOCKING_ENV_PATH"
    echo "      Rode: mamba env create -f envs/docking.yml --yes"
    ERRORS=$((ERRORS+1))
else
    echo "[1/3] Instalando meeko no docking (--no-deps, evita recompilação openbabel)..."
    if mamba run -n docking pip install "meeko>=0.5.0" --no-deps --quiet 2>&1; then
        echo "      meeko OK"
    else
        echo "      AVISO: meeko pode não ter instalado — verificar manualmente"
    fi
fi
echo ""

# ── 2. HADDOCK3 no docking env ────────────────────────────────────────────
HADDOCK_DIR="${HOME}/gromacs/haddock3"
if [ ! -d "$DOCKING_ENV_PATH" ]; then
    echo "[2/3] SKIP: ambiente docking não existe (ver erro acima)"
elif [ -d "$HADDOCK_DIR" ]; then
    echo "[2/3] Linkando HADDOCK3 no ambiente docking..."
    # pip pode emitir warnings de conflito de openbabel — não é erro fatal
    if mamba run -n docking pip install -e "$HADDOCK_DIR" --quiet 2>&1 | grep -v "^ERROR: pip" | grep -v "openbabel"; then
        true
    fi
    # Verificar se haddock3 importa
    if mamba run -n docking python -c "import haddock3; print('HADDOCK3 OK')" 2>/dev/null; then
        echo "      HADDOCK3 OK (instalado de $HADDOCK_DIR)"
    else
        echo "      AVISO: HADDOCK3 instalado mas import falhou — pode ser problema de deps openbabel"
    fi
else
    echo "[2/3] AVISO: $HADDOCK_DIR não encontrado — instale HADDOCK3 manualmente:"
    echo "      git clone https://github.com/haddocking/haddock3 ~/gromacs/haddock3"
    echo "      mamba run -n docking pip install -e ~/gromacs/haddock3"
fi
echo ""

# ── 3. AlphaFold3 no structure env ────────────────────────────────────────
STRUCTURE_ENV_PATH="${CONDA_BASE}/envs/structure"
AF3_DIR="${HOME}/gromacs/alphafold3"

if [ ! -d "$STRUCTURE_ENV_PATH" ]; then
    echo "[3/3] ERRO: ambiente 'structure' não encontrado"
    echo "      Rode: mamba env create -f envs/structure.yml --yes"
    ERRORS=$((ERRORS+1))
else
    # Verificar versão Python do structure env
    PY_VER=$(mamba run -n structure python -c "import sys; print(sys.version_info.minor)" 2>/dev/null || echo "0")
    if [ "$PY_VER" -lt 12 ]; then
        echo "[3/3] AVISO: structure env tem Python 3.${PY_VER} — AlphaFold3 exige >=3.12"
        echo "      Recriar o env com python=3.12:"
        echo "        mamba env remove -n structure --yes"
        echo "        mamba env create -f envs/structure.yml --yes"
        echo "      (envs/structure.yml já foi atualizado para python=3.12)"
        ERRORS=$((ERRORS+1))
    elif [ -d "$AF3_DIR" ]; then
        echo "[3/3] Linkando AlphaFold3 no ambiente structure (Python 3.${PY_VER})..."
        if mamba run -n structure pip install -e "$AF3_DIR" --quiet 2>&1; then
            echo "      AlphaFold3 OK (instalado de $AF3_DIR)"
        else
            echo "      ERRO ao instalar AlphaFold3"
            ERRORS=$((ERRORS+1))
        fi
    else
        echo "[3/3] AVISO: $AF3_DIR não encontrado — instale AlphaFold3 manualmente:"
        echo "      git clone https://github.com/google-deepmind/alphafold3 ~/gromacs/alphafold3"
        echo "      mamba run -n structure pip install -e ~/gromacs/alphafold3"
    fi
fi
echo ""

# ── Verificação final ─────────────────────────────────────────────────────
echo "=== Verificação de imports críticos ==="
check_import() {
    local env=$1; local code=$2; local label=$3
    printf "%-35s" "$label"
    result=$(mamba run -n "$env" python -c "$code" 2>&1)
    if echo "$result" | grep -q "OK\|CUDA\|cuda\|gpu\|GPU\|CudaDevice\|[0-9] device"; then
        echo "✅ $result"
    else
        echo "❌ FALHOU: $result" | head -1
    fi
}

check_import docking   "import openbabel; print('OK')"                      "docking  → openbabel:"
check_import docking   "import meeko; print('OK')"                          "docking  → meeko:"
check_import docking   "from vina import Vina; print('OK')"                 "docking  → vina:"
check_import structure "import foldseek; print('OK')" 2>/dev/null || \
  printf "%-35s%s\n"  "docking  → foldseek (CLI):"  "✅ CLI tool (sem import Python)"
check_import structure "import jax; d=str(jax.devices()); print('OK',d)"   "structure → jax+GPU:"
check_import structure "import alphafold3; print('OK')" 2>/dev/null || true

echo ""
if [ "$ERRORS" -gt 0 ]; then
    echo "=== $ERRORS erro(s) encontrado(s) — ver mensagens acima ==="
    exit 1
else
    echo "=== Todos os passos concluídos com sucesso ==="
fi
