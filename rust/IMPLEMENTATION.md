# Rust — Implementation Notes

Decision log for `mlp-main/rust`. Format per `code_optimization.md` step 7.
Só `main` é otimizada (ver `CLAUDE.md`).

## `#[cold]` na atualização do melhor movimento (2026-11-17)

### Problema

Em `search_swap` e `search_two_opt` (as duas inlined em `main`), o LLVM
converte o `if cost_new < cost_best { ... }` em `cmov`. Como o bloco
atualiza só o custo e dois índices (`I`, `J`), ele fica abaixo do limiar de
if-conversion. O `cmov` roda em toda iteração do loop interno e cria uma
dependência entre iterações em `cost_best`, `I` e `J`. Uma melhoria é rara,
então um branch seria previsto quase sempre e sairia mais barato. A
`search_reinsertion` atualiza três índices (`I`, `J`, `POS`) e já ficava
com branch.

É o mesmo problema que o Swift tinha (ver `swift/IMPLEMENTATION.md`).

### Mudança

Uma função vazia marcada como fria, chamada dentro dos 5 `if` de melhoria:

```rust
#[cold]
#[inline(never)]
fn cold() {}
```

A chamada não pode ser executada especulativamente, então o LLVM mantém o
branch, e o `#[cold]` marca o caminho como improvável. O `std::hint::cold_path`
faria o mesmo, mas ainda é instável no Rust 1.93.

### Verificação (rustc 1.93.0, `cargo build --release --features flat`)

`cmov` contados com `objdump`:

| Função                | Antes | Depois |
|-----------------------|------:|-------:|
| `main` (swap, 2-opt)  |    11 |      5 |
| `search_reinsertion`  |     0 |      0 |

Os 5 `cmov` que sobram em `main` ficam fora das buscas (leitura da entrada e
`core::str`).

Custo final idêntico em berlin52, bier127 e gr202. `perf stat`, `taskset -c 3`,
melhor de 3 execuções:

| Instância | Instruções | Ciclos | Tempo        |
|-----------|-----------:|-------:|-------------:|
| bier127   |      −5,3% |  −6,2% | 5,32 → 4,97 s |
| gr202     |      −5,3% |  −7,1% | 26,5 → 24,4 s |

A máquina tinha outros processos disputando o mesmo núcleo, por isso o
melhor de 3 em vez da mediana. A contagem de instruções é determinística e
bate com a medição anterior (investigação de 2026-11-16: −7,3% e −4,1% de
ciclos em bier127 e gr202).
