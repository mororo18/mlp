# Julia — Implementation Notes

Decision log for `mlp-main/julia`. Format per `code_optimization.md` step 7.
Só `main` é otimizada (ver `CLAUDE.md`).

## Contador `nimp` na atualização do melhor movimento (2026-11-17)

### Problema

Em `search_swap` e `search_two_opt`, o LLVM converte o
`if cost_new < cost_best ... end` em `cmov`. Como o bloco atualiza só o
custo e dois índices (`I`, `J`), ele fica abaixo do limiar de
if-conversion. O `cmov` roda em toda iteração do loop interno e cria uma
dependência entre iterações. Uma melhoria é rara, então um branch seria
previsto quase sempre e sairia mais barato. A `search_reinsertion` atualiza
três índices (`I`, `J`, `POS`) e já ficava com branch.

É o mesmo problema que o Swift e o Rust tinham (ver os `IMPLEMENTATION.md`
deles).

### Alternativas testadas

Julia 1.12 não tem diretiva equivalente ao `_slowPath` do Swift ou ao
`#[cold]` do Rust:

- `llvm.expect` via `Base.llvmcall`: o Julia descarta a dica e o `cmov`
  continua.
- Chamar uma função ou um `asm` vazio dentro do `if`: remove o `cmov`, mas o
  LLVM passa a supor que a chamada pode escrever na memória e recarrega os
  dados depois dela (+14% de instruções, mais lento).
- `asm sideeffect` vazio com o atributo `memory(none)`, via `llvmcall`:
  funciona (−4,5% de ciclos no gr202), mas depende de IR do LLVM escrita à
  mão, que muda entre versões.
- `JULIA_LLVM_ARGS=-x86-cmov-converter-force-all=true`: cerca de −3% de
  tempo, mas é uma flag global do backend.

### Mudança

Um contador de melhorias, `nimp`, incrementado dentro do `if` e testado
antes de aplicar o movimento:

```julia
if cost_new < cost_best
    cost_best = cost_new - data.EPSILON
    I = i
    J = j
    nimp += 1
end
...
if nimp > 0 && cost_best < solut.cost - data.EPSILON
```

Com três variáveis inteiras vivas no bloco, ele passa do limiar de
if-conversion e o LLVM mantém o branch. O `nimp > 0` garante que o contador
é lido, para o compilador não eliminá-lo. É Julia puro, sem IR.

### Verificação (Julia 1.12.5, `--check-bounds=no --math-mode=fast`)

Nenhum `cmov` depende mais da comparação de ponto flutuante (`code_native`):

| Função           | `cmov` após comparação FP, antes | depois |
|------------------|---------------------------------:|-------:|
| `search_swap`    |                                2 |      0 |
| `search_two_opt` |                                2 |      0 |

Custo final idêntico em berlin52, bier127 e gr202. `perf stat`,
`taskset -c 3`, melhor de 3 execuções:

| Instância | Instruções | Ciclos | Tempo         |
|-----------|-----------:|-------:|--------------:|
| bier127   |      −0,8% |  −0,5% | 6,25 → 6,22 s |
| gr202     |      −1,9% |  −4,4% | 28,9 → 27,7 s |

A máquina tinha outros processos disputando o mesmo núcleo, por isso o
melhor de 3. O gr202 bate com a medição anterior (investigação de
2026-11-16: −4,6% de ciclos, −1,9% de instruções). No bier127 o efeito
ficou dentro do ruído.
