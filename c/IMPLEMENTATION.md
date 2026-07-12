# C — Implementation Notes

Decision log for `mlp-main/c`. Format per `code_optimization.md` step 7.
Only `main` é otimizada (ver `CLAUDE.md`).

## MATRIX vs FLAT — default trocado pra FLAT (2026-07-12)

`types.h` já implementava dois layouts alternativos pro campo `seq` de
`tSolution`, selecionados por macro de pré-processador: `MATRIX` (array
de ponteiros, alocação por linha — fragmentado) e `FLAT` (array 1D
contíguo, indexado via `to_1D(i,j,size)`). O `makefile` usava `MATRIX`
como default — ou seja, o binário de `main` que o harness builda e mede
sempre usou o layout fragmentado, contradizendo a regra de performance
do `CLAUDE.md` ("sem fragmentação de memória, arrays contíguos como
estrutura primária"). Não achado nenhum commit/doc explicando essa
escolha — provável esquecimento de quando o layout `FLAT` foi
implementado sem trocar o default.

### Microbenchmark de desempate (2026-07-12)

Metodologia: build real do binário nas duas configurações (`make` com
`MACRO=-DMATRIX` vs `-DFLAT`), rodando o algoritmo completo (não um
microbenchmark sintético isolado) contra instâncias reais do conjunto
do TCC — **leve** (`burma14`, 14 nós) e **intensiva** (`pr299`, 299
nós, maior instância do conjunto oficial). N=10 pra `burma14`, N=15
pra `pr299` (amostra maior por pedido do autor, já que a primeira
rodada com N=5 deu resultado ambíguo/dentro do ruído).

| instância | MATRIX (mediana) | FLAT (mediana) | vencedor |
|---|---|---|---|
| burma14 (leve) | 0.0012s | 0.0012s | FLAT (~1.6%, quase ruído nessa escala) |
| pr299 (intensiva) | 80.01s (σ=5.69) | 77.41s (σ=3.36) | **FLAT, ~3.3%** |

**FLAT vence nas duas escalas** (mesmo padrão em C++, ver
`cplusplus/IMPLEMENTATION.md`). A vantagem é modesta (~3-4% na escala
intensiva) mas consistente — a primeira rodada com N=5 deu MATRIX
"vencendo" por ruído (diferença de ~6s contra desvio-padrão de ~6s,
~1 sigma, não significativo); com N=15 o sinal ficou claro.

### Ação

- `makefile`: `MACRO` default trocado de `-DMATRIX -DNDEBUG` pra
  `-DFLAT -DNDEBUG`. Alvo `flat:` mantido (agora redundante com o
  default, por clareza). Alvo `matrix:` adicionado pra referência
  futura, caso seja preciso comparar de novo.
- Corretude validada: `COST` bate com o valor esperado em `burma14`
  (20315), `att48` (209320) e `pr299` (6556628) após a troca.

### Nota metodológica

Diferente do ciclo do Lua (que usou microbenchmarks sintéticos
isolados pra decompor causas), aqui o desempate foi feito rodando o
**algoritmo completo real** nas duas configurações — já que o código
das duas variantes já existia pronto atrás da macro, não havia
necessidade de isolar a estrutura de dados num teste sintético
separado; o teste end-to-end já é o teste mais direto e realista
possível pra essa pergunta específica.
