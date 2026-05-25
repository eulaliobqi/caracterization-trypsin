#!/usr/bin/env bash
# post-execute.sh — Acionado após cada fase concluída
set -euo pipefail

PHASE="${1:-unknown}"
TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
LEARNINGS="$(dirname "$0")/../../LEARNINGS.md"

echo "📊 Post-execute: fase ${PHASE} concluída em ${TIMESTAMP}"

# Registra a conclusão da fase no LEARNINGS.md
if [[ -f "$LEARNINGS" ]]; then
    echo "" >> "$LEARNINGS"
    echo "## EXEC-$(date '+%Y%m%d%H%M') [${TIMESTAMP}] Fase ${PHASE} concluída" >> "$LEARNINGS"
    echo "- Status: completed" >> "$LEARNINGS"
    echo "- Outputs em: results/$(printf '%02d' $PHASE)_*/" >> "$LEARNINGS"
fi

# Verificar se há outputs esperados
case "$PHASE" in
    1) CHECK_DIR="results/01_qc" ;;
    2) CHECK_DIR="results/02_orfs" ;;
    3) CHECK_DIR="results/03_trypsin_ids" ;;
    4) CHECK_DIR="results/04_complete" ;;
    5) CHECK_DIR="results/05_primary" ;;
    6) CHECK_DIR="results/06_phylogeny" ;;
    7) CHECK_DIR="results/07_structures" ;;
    8) CHECK_DIR="results/08_validation" ;;
    9) CHECK_DIR="results/09_docking" ;;
    10) CHECK_DIR="results/10_md" ;;
    11) CHECK_DIR="results/11_report" ;;
    *) CHECK_DIR="results" ;;
esac

if [[ -d "$CHECK_DIR" ]] && [[ -n "$(ls -A "$CHECK_DIR")" ]]; then
    echo "✅ Outputs encontrados em $CHECK_DIR"
else
    echo "⚠️  AVISO: Diretório $CHECK_DIR vazio ou não encontrado após fase $PHASE"
fi

echo "Phase ${PHASE} completed at ${TIMESTAMP}" >> pipeline_execution.log
