# Go — Implementation Notes

Decision log for `mlp-main/go`. Format per `code_optimization.md` step 7.
Only `main` é otimizada (ver `CLAUDE.md`).

## Contiguidade de `seq` — adotado (2026-07-15)

### Contexto

`tSolution.seq` era `[][]tSubseq` (slice-de-slices) — cada linha
(`seq[i]`) uma alocação separada no heap, fragmentado mesmo em `main`.
Era o único gap real na técnica G2 de contiguidade: Go tem arrays
multidimensionais nativos que poderiam ser achatados, ao contrário de
Fortran/Julia (contíguos nativamente, sem escolha) — ver
`TECNICAS_G1_G2.md` §1, onde isso já estava marcado como pendência.

### O que foi testado

Achatado pra `[]tSubseq` único, tamanho `(dimen+1)*(dimen+1)`, indexado
via `i*dim+j` através de um método helper (`solut.at(i, j) *tSubseq`),
substituindo as 41 linhas de acesso `solut.seq[i][j]` no arquivo
inteiro (troca mecânica, mesmo padrão usado no C/C++).

### Microbenchmark (leve + intensiva)

- **Leve** (`burma14`, N=15): 0.00301s ± 0.00064 (fragmentado) vs.
  0.00289s ± 0.00055 (flat) — inconclusivo, instância trivial.
- **Intensiva** (`pr299`, N=5): **257.95s ± 16.73** (fragmentado) vs.
  **243.22s ± 22.92** (flat) — flat 5.7% mais rápido em média, mas as
  distribuições **se sobrepõem** (mínimo do fragmentado entra na faixa
  do flat) — sinal direcionalmente consistente, não estatisticamente
  limpo no N testado (mesmo padrão do MATRIX/FLAT original em C, que
  precisou de N=15 pra sinal ficar limpo).
- **Corretude**: `COST` idêntico nas duas variantes, todas as execuções
  (burma14 e pr299).

### Decisão

**Adotado mesmo com sinal não-decisivo** — decisão do autor em
2026-07-15: aceitar o resultado direcionalmente favorável (~5.7%) sem
estender pra N=15, dado o custo de tempo desproporcional (~200s+ por
execução em `pr299`). Fecha o gap de contiguidade que faltava em Go.

**Consequência pro TCC**: os dados já coletados na campanha `main`
completa (90 linhas de Go, `mlp_testao/go.csv`) refletem a versão
*fragmentada* antiga — precisam ser recoletados pra refletir o
código atual (ver `CRONOGRAMA.md`).
