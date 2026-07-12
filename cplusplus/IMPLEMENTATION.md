# C++ — Implementation Notes

Decision log for `mlp-main/cplusplus`. Format per `code_optimization.md`
step 7. Only `main` é otimizada (ver `CLAUDE.md`).

## MATRIX vs FLAT — default trocado pra FLAT (2026-07-12)

Mesmo achado e mesma correção do C (ver `c/IMPLEMENTATION.md` pra
metodologia completa): `main.cpp` também implementa os dois layouts
(`#ifdef FLAT` / `#elif defined(MATRIX)`), mas tinha `#define MATRIX`
hardcoded no topo do arquivo — o caminho `FLAT` existia como código
morto, nunca compilado pelo build default.

### Resultado do microbenchmark de desempate

Mesmo teste do C (algoritmo completo real, `burma14` leve e `pr299`
intensiva, N=10/N=15):

| instância | MATRIX (mediana) | FLAT (mediana) | vencedor |
|---|---|---|---|
| burma14 (leve) | 0.0015s | 0.0012s | FLAT (~18%) |
| pr299 (intensiva) | 81.21s (σ=1.67) | 77.75s (σ=1.00) | **FLAT, ~4.3%** |

FLAT vence nas duas escalas, mesma direção do C.

### Ação

- `main.cpp`: `#define MATRIX` (linha 21) trocado pra `#define FLAT`.
- Corretude validada: `COST` bate com o valor esperado em `burma14`
  (20315), `att48` (209320) e `pr299` (6556628) após a troca.
