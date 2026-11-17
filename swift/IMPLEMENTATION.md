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

## `bestJ` usado no `searchOrOpt` (2026-11-17)

### Problema

No `searchOrOpt`, o LLVM converte o `if costNew < costBest { ... }` em `cmov`,
executado em toda iteração do loop interno. Isso cria uma dependência entre
iterações em `costBest`, `bestI` e `bestPos`. Uma melhoria é rara, então um
branch seria previsto quase sempre e sairia mais barato. É a busca mais
quente (cerca de 2/3 do tempo no bier127). O bloco atualiza só dois índices
e fica abaixo do limiar de if-conversion; nas outras linguagens a busca
equivalente guarda três (`I`, `J`, `POS`) e já fica com branch.

O f2045ab tentou resolver com uma terceira variável (`bestJ = j`). Mas o
`bestJ` só era lido num `assert`, que some no build de release. O compilador
eliminou a variável e o `cmov` continuou. Por isso quem rodou os testes não
viu ganho.

### Mudança

O tamanho do bloco passado ao `reinsert` agora sai de `bestJ`
(`bestJ - bestI + 1`), no lugar do `assert`. Com o `bestJ` lido, o bloco
atualiza três variáveis e o LLVM mantém o branch.

`searchSwap` e `searchTwoOpt` continuam com `cmov`. Entre 2026-11-16 e
2026-11-17 elas usaram `_slowPath` (d5fd2c3), revertido para manter só a
mudança no Or-opt, sem dicas explícitas ao compilador.

### Verificação (Swift 6.1.3, `-c release -Ounchecked`)

Número de `cmov` em cada função, contado com `objdump`:

| Versão                         | Swap | TwoOpt | OrOpt |
|--------------------------------|-----:|-------:|------:|
| b83d040 (antes)                |    4 |      2 |     4 |
| f2045ab (`bestJ` morto)        |    4 |      2 |     4 |
| esta mudança (`bestJ` usado)   |    4 |      2 |     0 |

Custo final idêntico em berlin52, bier127 e gr202.

O mesmo resultado (OrOpt 4 → 0) apareceu no Swift 5.10.1, 6.0.3, 6.1.3 e
6.2.4. Mesmo assim, o limiar de if-conversion pode mudar entre versões do
LLVM.

Os tempos ainda não foram medidos com confiança. Com a máquina carregada
(load ~9), as versões sem `cmov` no `searchOrOpt` ficaram de 15% a 25% abaixo
do b83d040; falta remedir numa máquina sem carga.
