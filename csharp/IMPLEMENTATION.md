# C# — Implementation Notes

Decision log for `mlp-main/csharp`. Format per `code_optimization.md` step 7.
Only `main` é otimizada (ver `CLAUDE.md`).

## Bounds-check via `unsafe`/`fixed` — testado e descartado (2026-07-14)

### Contexto

`tSolution.GetSeq`/`SetSeq` (`tSolution.cs`) fazem indexação em `seq`
(`double[][]`, jagged — cada linha achatada em `j*3+k`), com
bounds-checking normal do CLR nas duas indexações (`seq[i]` e
`[j*3+k]`). Outras linguagens do benchmark já desligam bounds-check
como técnica G2 (Rust via `get_unchecked`, Julia via
`--check-bounds=no`, Go via `-gcflags=-B` — ver `TECNICAS_G1_G2.md`).
C# nunca tinha essa técnica testada; nula hipótese natural era "C#
também ganha com bounds-check off".

### O que foi testado

Variante `unsafe`, condicional via `#if UNSAFE_SEQ` (mesmo padrão de
macro do C/C++), pinning a linha (`fixed (double* p = seq[i])`) e
acessando via ponteiro (`p[j*3+k]`) em vez do indexador gerenciado —
isolado só ao bounds-check da indexação interna (`j*3+k`), sem mexer
na estrutura jagged em si (não confundir com a técnica de
contiguidade, que é um eixo separado — ver `TECNICAS_G1_G2.md` §1).

```csharp
#if UNSAFE_SEQ
[MethodImpl(MethodImplOptions.AggressiveInlining)]
public unsafe double GetSeq(int i, int j, int k) {
    fixed (double* p = seq[i]) { return p[j*3 + k]; }
}
[MethodImpl(MethodImplOptions.AggressiveInlining)]
public unsafe void SetSeq(int i, int j, int k, double v) {
    fixed (double* p = seq[i]) { p[j*3 + k] = v; }
}
#else
// original, indexador gerenciado
#endif
```

`<AllowUnsafeBlocks>true</AllowUnsafeBlocks>` no `.csproj`, alternando
build via `dotnet build -p:DefineConstants=UNSAFE_SEQ`.

### Microbenchmark (mesmo workflow do C/C++/Python: leve + intensiva, N amostras)

- **Leve** (`burma14`, N=15): checked 0.0887s ± 0.0346, unsafe 0.1027s
  ± 0.0341 — diferença menor que meio desvio-padrão, inconclusivo
  (esperado — instância trivial, tempo de processo/JIT domina).
- **Intensiva** (`pr299`, N=8): checked **227.97s ± 19.05s**, unsafe
  **416.76s ± 31.88s** — **unsafe é 1.83x mais lento**, distribuições
  sem sobreposição nenhuma (checked máx. 266.34s, unsafe mín. 384.27s).
  Sinal forte, não é ruído.
- **Corretude**: `COST: 6556628` idêntico nas duas variantes, todas as
  16 execuções.
- **Reprodução independente** (2026-07-14, mesmo dia, suspeita de que a
  1ª medição tivesse CPU concorrente elevada): repetido do zero — build
  limpo, `distance_matrix` de `pr299` regerado, CPU confirmada ociosa
  antes de começar (`uptime`/`top`), load average logado antes de cada
  execução individual (ficou entre 0.8 e 4.0 durante o teste, sempre
  bem abaixo dos 8 cores da máquina). Resultado: checked **202.07s ±
  8.80s**, unsafe **396.52s ± 19.53s** — **1.96x mais lento**, ainda
  mais forte que a primeira medição, distribuições de novo sem
  sobreposição (checked máx. 213.13s, unsafe mín. 371.19s). `COST`
  idêntico nas 16 execuções. Descarta a hipótese de artefato de
  contenção de CPU — o achado é real e reproduzível.

### Por que piorou (hipótese, não medida separadamente)

`GetSeq`/`SetSeq` são chamados de dezenas de sites espalhados por
`GILS_RVND.cs`, não dentro de um laço único fechado — cada chamada
individual entra e sai do próprio bloco `fixed`, ou seja, o *pinning*
da linha (`seq[i]`) é feito e desfeito a cada chamada, não uma vez só
por fora de um laço quente. O overhead de pinning por chamada
(interação com o GC) parece dominar amplamente a economia de pular uma
única checagem de bounds — que o JIT moderno já executa de forma
barata (branch previsível) e possivelmente já elimina em alguns casos
via análise de intervalo. Não medido isoladamente (não é o foco: o
resultado agregado já é conclusivo o suficiente pra decisão).

**Possível variante melhor, não testada**: pinning uma única vez por
fora de uma chamada quente (ex.: dentro do próprio RVND, mantendo o
ponteiro pela duração de uma passada inteira) em vez de por
Get/Set individual — exigiria mudar a assinatura de chamada em todo o
código, intervenção bem maior que o testado aqui. Fica como possível
trabalho futuro, não perseguido nesta sessão.

### Decisão

**Não adotado.** Bounds-check-off via `unsafe`/`fixed` no granularidade
testada (por chamada individual de `GetSeq`/`SetSeq`) é uma regressão
clara e mensurável em C#, ao contrário de Rust/Julia/Go, onde a mesma
categoria de técnica (G2, "bounds-check desligado") ajudou. Código de
teste descartado (não commitado); `tSolution.cs` de `main` permanece
inalterado, com bounds-checking normal do CLR.

**Consequência pro TCC**: é um exemplo direto e citável de "a mesma
categoria de técnica (G2) não necessariamente se transfere entre
linguagens da mesma forma" — a técnica de bounds-check-off é boa em
linguagens onde a garantia removida tem custo real por si (Rust panics,
checagem de índice do Julia), mas em C# o *mecanismo* disponível pra
remover essa garantia (`unsafe`/`fixed`, pinning por chamada) tem custo
próprio que supera o benefício, pelo menos na forma testada aqui.
Relevante pra RQ3 (as técnicas de otimização transferem entre
linguagens?) e pro ponto (a) da tese (a melhora de técnicas não-óbvias
vale "na medida que a linguagem permite" — aqui a linguagem *permite*
mecanicamente, mas o mecanismo específico não compensa).
