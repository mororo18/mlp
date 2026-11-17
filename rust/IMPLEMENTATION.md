# Rust — Implementation Notes

Decision log for `mlp-main/rust`. Format per `code_optimization.md` step 7.
Só `main` é otimizada (ver `CLAUDE.md`).

## `cold_path` na atualização do melhor movimento (2026-11-17)

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

`std::hint::cold_path()` no início dos 5 `if` de melhoria. Ele marca o
caminho como improvável, e o LLVM mantém o branch.

`cold_path` é estável desde o Rust 1.95, por isso o `Cargo.toml` declara
`rust-version = "1.95"`, e o toolchain foi atualizado de 1.93.0 para
1.100.0. A primeira versão desta mudança (7dfcaf1), ainda no 1.93, chamava
uma função vazia marcada com `#[cold] #[inline(never)]`, que tem o mesmo
efeito.

### Verificação (rustc 1.100.0, `cargo build --release --features flat`)

`cmov` que dependem da comparação de custo, contados com `objdump`:

| Versão                 | swap | 2-opt |
|------------------------|-----:|------:|
| sem hint               |    4 |     2 |
| `cold_path`            |    0 |     0 |

Os `cmov` que sobram em `main` ficam fora das buscas (leitura da entrada,
`core::str` e o laço do RVND).

Custo final idêntico em bier127 e gr202. `perf stat`, `taskset -c 7`,
mediana de 5 execuções, em relação ao 1.100.0 sem hint:

| Variante                  | Instância | Instruções | Ciclos | Tempo           |
|---------------------------|-----------|-----------:|-------:|----------------:|
| 1.100.0 + `cold_path`     | bier127   |      −6,6% |  −4,8% | 4,00 → 3,76 s   |
| 1.100.0 + `cold_path`     | gr202     |      −6,7% |  −6,7% | 21,8 → 20,1 s   |
| 1.93.0 + `#[cold]` (7dfcaf1) | bier127 |     +17,7% |  +2,7% | 4,00 → 4,55 s   |
| 1.93.0 + `#[cold]` (7dfcaf1) | gr202   |     +18,0% |  +3,6% | 21,8 → 23,4 s   |

O hint dá de 5% a 8% no 1.100, como dava no 1.93. A troca de compilador
também pesa: mesmo sem hint, o 1.100 executa cerca de 15% menos instruções
que o 1.93 com hint.
