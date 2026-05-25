---
description: Compacta o contexto da conversa preservando estado essencial (acione quando contexto > 60%)
---

Invoque `token-optimizer` com as seguintes instruções:

Compacte o contexto desta conversa preservando OBRIGATORIAMENTE:
1. **Fase atual** do pipeline e seu status
2. **Último entry do LEARNINGS.md** (erro ou decisão mais recente)
3. **Próximos 3 passos** concretos a executar
4. **Outputs produzidos** nesta sessão (caminhos + N resultados)
5. **Parâmetros customizados** ativos (qualquer override de params.*)

Podem ser compactados (→ referência de arquivo):
- Logs de execução de ferramentas
- Conteúdo de arquivos já lidos integralmente
- Histórico de erros já resolvidos
- Outputs longos de comandos anteriores

Meta: redução de ≥ 40% no tamanho do contexto.
