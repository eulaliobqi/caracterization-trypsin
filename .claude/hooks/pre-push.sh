#!/usr/bin/env bash
# pre-push.sh — Valida pipeline com stub-run antes de push
set -euo pipefail

echo "🚀 Pre-push: rodando stub-run do pipeline..."

if ! command -v nextflow &>/dev/null; then
    echo "  ⚠️  Nextflow não disponível localmente — pulando stub-run"
    echo "  AVISO: Execute 'nextflow run nextflow/main.nf -stub-run -profile test' no servidor antes de usar"
    exit 0
fi

cd "$(dirname "$0")/../../nextflow" 2>/dev/null || cd nextflow

echo "  → nextflow run main.nf -stub-run -profile test..."
nextflow run main.nf \
    -stub-run \
    -profile test \
    --outdir ../results_stub \
    --input_fasta ../tests/data/test_assembly.fasta \
    2>&1

echo "  → Limpando resultados de teste..."
rm -rf ../results_stub ../work

echo "✅ Stub-run passou — OK para push"
