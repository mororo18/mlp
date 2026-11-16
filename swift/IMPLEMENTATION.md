# Swift — Implementation Notes

Decision log for `mlp-main/swift`. Format per `code_optimization.md` step 7.
Só `main` é otimizada (ver `CLAUDE.md`).

## Matriz de custo 2D `[[Double]]` (2026-11-11, f2045ab)

Veio da ablação da tradução C++ → Swift por LLM (relatório
`SWIFT_PERFORMANCE_ANALYSIS.md`). Com a matriz 1D (`matrix[i*dimension + j]`),
cada consulta nos loops internos faz um `imul`, porque o índice da linha vem
da permutação em tempo de execução. Com `[[Double]]`, o LLVM tira o ponteiro
da linha (`matrix[s[i]]`) de dentro do loop interno. Ganho medido no
relatório: cerca de 3% (bier127, gr202).

## `_slowPath` na atualização do melhor movimento (2026-11-16)

### Problema

Nas três buscas (`searchSwap`, `searchTwoOpt`, `searchOrOpt`), o LLVM
converte o `if costNew < costBest { ... }` em `cmov`, executado em toda
iteração do loop interno. Isso cria uma dependência entre iterações em
`costBest` e nos índices do melhor movimento. Como uma melhoria é rara, um
branch comum seria previsto quase sempre e sairia mais barato.

O f2045ab tentou resolver o caso do `searchOrOpt` com uma terceira variável
(`bestJ = j`), que no relatório passou o bloco do limiar de if-conversion e
deu ~10%. Mas o `bestJ` só era lido num `assert`, que some no build de
release. O compilador eliminou a variável e o `cmov` continuou. Por isso
quem rodou os testes não viu ganho.

### Mudança

`if _slowPath(costNew < costBest)` nos 5 pontos de atualização. `_slowPath`
é público na stdlib do Swift e vira `llvm.expect`, que impede a
if-conversion. O `bestJ` do `searchOrOpt` foi removido.

### Verificação (Swift 6.1.3, `-c release -Ounchecked`)

Número de `cmov` em cada função, contado com `objdump`:

| Versão                         | Swap | TwoOpt | OrOpt |
|--------------------------------|-----:|-------:|------:|
| b83d040 (antes)                |    4 |      2 |     4 |
| f2045ab (`bestJ` morto)        |    4 |      2 |     4 |
| f2045ab + `bestJ` usado        |    4 |      2 |     0 |
| esta mudança (`_slowPath`)     |    0 |      0 |     0 |

Custo final idêntico em berlin52, bier127 e gr202.

Os tempos ainda não foram medidos com confiança. A máquina estava com load
alto (~9 em 8 CPUs, governor `powersave`), e o tempo do bier127 variou entre
5 e 7,7 s numa mesma variante. Nessas medições, as versões sem `cmov` no
`searchOrOpt` ficaram de 15% a 25% abaixo do b83d040, mas falta remedir numa
máquina sem carga.

O resultado depende da versão do compilador: o limiar de if-conversion muda
entre versões do LLVM, então a mesma fonte pode gerar `cmov` numa versão e
branch em outra.
