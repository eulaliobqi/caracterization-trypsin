#!/usr/bin/env bash
# pre-commit.sh — Validações obrigatórias antes de qualquer commit
set -euo pipefail

echo "🔍 Pre-commit checks — trypsin-agemmatalis-structural"

# 1. Lint arquivos Nextflow (.nf)
NF_FILES=$(find . -name "*.nf" -not -path "./.git/*" | head -1)
if [[ -n "$NF_FILES" ]]; then
    if command -v nf-core &>/dev/null; then
        echo "  → nf-core lint..."
        nf-core lint --release 2>&1 || { echo "❌ nf-core lint falhou"; exit 1; }
    else
        echo "  ℹ️  nf-core não disponível — pulando lint Nextflow"
    fi
fi

# 2. Lint Python
PY_FILES=$(find . -name "*.py" -not -path "./.git/*" | head -1)
if [[ -n "$PY_FILES" ]]; then
    if command -v ruff &>/dev/null; then
        echo "  → ruff check..."
        ruff check nextflow/bin/ tests/ 2>&1 || { echo "❌ ruff falhou"; exit 1; }
    elif command -v pylint &>/dev/null; then
        pylint nextflow/bin/*.py --disable=C,R 2>&1 || { echo "❌ pylint falhou"; exit 1; }
    else
        echo "  ℹ️  ruff/pylint não disponível — pulando lint Python"
    fi
fi

# 3. Validar YML de environments
for env in envs/*.yml; do
    [[ -f "$env" ]] || continue
    python3 -c "import yaml; yaml.safe_load(open('$env'))" \
        || { echo "❌ YAML inválido: $env"; exit 1; }
done

# 4. Verificar secrets hardcoded
if grep -rE "(API_KEY|PASSWORD|TOKEN|SECRET)=['\"][^'\"]+['\"]" \
     --exclude-dir=.git --include="*.py" --include="*.sh" \
     --include="*.nf" --include="*.config" . > /dev/null 2>&1; then
    echo "❌ Possível secret hardcoded detectado! Revise antes de commitar."
    exit 1
fi

# 5. Verificar paths absolutos pessoais (exceto projectDir/baseDir)
if grep -rE "/home/[a-z]|/Users/[A-Z]" \
     --exclude-dir=.git --include="*.nf" --include="*.config" . > /dev/null 2>&1; then
    echo "❌ Path absoluto pessoal detectado em .nf/.config — use \${projectDir} ou params"
    exit 1
fi

echo "✅ Pre-commit OK — $(date '+%Y-%m-%d %H:%M')"
