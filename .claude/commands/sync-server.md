---
description: Sincroniza local → GitHub → Servidor Debian para execução
---

Sincroniza o pipeline do ambiente local para o servidor Debian de execução.

## Procedimento

1. **Verifique status local:**
   ```bash
   git status
   git log --oneline -5
   ```

2. **Se houver código novo:** invoque `code-reviewer` antes de continuar

3. **Push para GitHub:**
   ```bash
   git push origin dev
   ```

4. **Acesse o servidor:**
   ```bash
   ssh ${SERVER_USER}@${SERVER_HOST}
   ```

5. **No servidor — pull e setup:**
   ```bash
   cd ~/trypsin-agemmatalis-structural
   git pull origin dev
   
   # Se environments novos foram adicionados:
   for env in envs/*.yml; do
       mamba env create -f "$env" --no-default-packages -y || \
       mamba env update -f "$env" -y
   done
   ```

6. **Rodar pipeline no servidor:**
   ```bash
   nextflow run nextflow/main.nf \
       -profile debian \
       --input_fasta data/raw/trinity_assembly.fasta \
       --outdir results \
       -resume \
       -bg  # background
   tail -f .nextflow.log
   ```

7. **Reportar:** hash do commit deployado + status do Nextflow

## Variáveis de configuração (defina no .env local)
```bash
export SERVER_USER=seu_usuario
export SERVER_HOST=ip_ou_hostname_do_servidor
```
