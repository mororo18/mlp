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

### Validação da hipótese de pinning (2026-08-02)

A hipótese acima ("overhead de pinning por chamada domina a economia do
bounds-check") estava marcada como não medida isoladamente. Validada com
dois experimentos fora do solver real, pro artigo companion
(`doc/artigo-linguagens/main.tex`, §5.3):

**Microbenchmark isolado** (`net8.0`, Release, array de 1M `double`,
4M índices "aleatórios" reciclados 25x, N=10 processos por variante,
CPU com ruído de fundo mas sinal limpo — distribuições sem sobreposição
entre percall e as outras duas):

| variante | ns/acesso (média ± dp) | vs. checked |
|---|---|---|
| indexador gerenciado (checked) | 11,09 ± 1,63 | 1,00x |
| `unsafe`/`fixed` por chamada | 26,97 ± 2,65 | **2,43x** |
| `unsafe`, pin único fora do laço | 9,91 ± 1,94 | 0,89x |

Confirma a hipótese: pin único (fora do laço quente) é estatisticamente
igual ou mais rápido que o indexador checked; o problema é
especificamente re-pinar a cada chamada, não o acesso via ponteiro em
si. Ordem de grandeza consistente com o 1,83x–1,96x medido no solver
real (a diferença absoluta é esperada — o microbenchmark isola só o
acesso, sem o resto do trabalho do GILS-RVND diluindo o overhead).

**Assembly gerado** (Godbolt, compilador `dotnet80csharpcoreclr`,
Release/FullOpts, mesmos dois métodos `GetChecked`/`GetUnsafePerCall`
do microbenchmark acima) — mostra o mecanismo exato:

```asm
; GetChecked — indexador gerenciado (7 instruções)
cmp    esi, dword ptr [rdi+0x08]      ; index < length?
jae    SHORT RNGCHKFAIL               ; 1 branch previsível
mov    eax, esi
vmovsd xmm0, qword ptr [rdi+8*rax+0x10]
ret
RNGCHKFAIL: call CORINFO_HELP_RNGCHKFAIL; int3

; GetUnsafePerCall — fixed por chamada (~17 instruções)
push   rbp
sub    rsp, 16
lea    rbp, [rsp+0x10]
mov    gword ptr [rbp-0x08], rdi
test   rdi, rdi                       ; null-check do array
je     SHORT ...
mov    rax, gword ptr [rbp-0x08]
cmp    dword ptr [rax+0x08], 0        ; length == 0 ?
jne    SHORT ...
...                                    ; (branch de array vazio → ptr null)
mov    rax, gword ptr [rbp-0x08]
cmp    dword ptr [rax+0x08], 0
jbe    SHORT RNGCHKFAIL
mov    rax, gword ptr [rbp-0x08]
add    rax, 16                        ; pula header p/ 1o elemento
movsxd rcx, esi
vmovsd xmm0, qword ptr [rax+8*rcx]    ; load — sem bounds check aqui
add    rsp, 16
pop    rbp
ret
```

O bounds check em si é 1 `cmp`+`jae` (branch previsível, quase de
graça). O `fixed` por chamada, pra materializar o ponteiro pinado
segundo a semântica do C# (`fixed` sobre array precisa retornar
`null` se o array for `null` ou vazio, daí os dois checks extras),
monta um frame de pilha e faz 2 branches adicionais mais aritmética de
ponteiro — mais caro que o único check que substitui, e isso é só o
custo determinístico por instrução. A interação com o GC (pin bloqueia
compactação, custo de bookkeeping se uma coleta ocorrer durante a
janela pinada) é um custo probabilístico *adicional* em cima desse,
consistente com a doc oficial do .NET sobre bounds-check elimination
no JIT e sobre custo de pinning, mas não a causa dominante — a
dominante é o overhead determinístico de instruções visível acima.

Scripts/fontes do microbenchmark e do assembly não commitados (mesma
convenção do teste original acima — resultado documentado aqui, código
de teste descartado).

**Confirmação de que o assembly acima é build otimizado** (pergunta do
autor, 2026-08-02): o dump do JIT do CoreCLR rotula os dois métodos
`(FullOpts)` (ex.: `Bench:GetChecked(double[],int):double (FullOpts):`)
— tag que o próprio JIT usa quando aplica otimização completa,
oposto de `MinOpts` (o que sairia em Debug, mesmo problema descoberto
antes nesta sessão pro build do solver via `dotnet build -c Release`).
Adicionalmente, o backend de C#/.NET do Compiler Explorer roda com
`buildConfig=Release` e desliga tiered compilation
(`DOTNET_TieredCompilation=0`), forçando o JIT a compilar já na tier
otimizada, sem passar por Tier0/QuickJit. Não usa `dotnet build -c
Release` por trás (invoca `csc.dll` direto, sem MSBuild), mas o efeito
relevante pra comparação — otimização completa aplicada, sem o
throttling que `-c Debug` causaria — está confirmado pela tag
`(FullOpts)` do próprio JIT, que é a evidência mais direta possível.

**Validação cruzada local com BenchmarkDotNet** (mesma dúvida,
2026-08-02): pra ter uma garantia mais forte e citável que "li o
source do Compiler Explorer no GitHub", reproduzido o mesmo par de
métodos (`GetChecked`/`GetUnsafePerCall`) com BenchmarkDotNet 0.15.8
(`dotnet run -c Release`, `.NET 8.0.25`, `RyuJIT`) usando
`[DisassemblyDiagnoser]`. Confirmado empiricamente que o
BenchmarkDotNet recusa rodar em build não-otimizado: `dotnet run -c
Debug` no mesmo projeto falha com "Benchmark was built without
optimization enabled (most probably a DEBUG configuration). Please,
build it in RELEASE." — validação automática embutida na ferramenta,
não depende de inspecionar config de terceiro.

Resultado (`BenchmarkDotNet.Artifacts/results/Bench-asm.md`):

```asm
; Bench.GetChecked()  — 33 bytes
push   rax
mov    rax,[rdi+8]        ; this.arr
mov    ecx,[rdi+10]       ; this.idx
cmp    ecx,[rax+8]
jae    short RNGCHKFAIL
vmovsd xmm0,qword ptr [rax+rcx*8+10]
ret
RNGCHKFAIL: call CORINFO_HELP_RNGCHKFAIL; int 3

; Bench.GetUnsafePerCall()  — 79 bytes
push   rbp
sub    rsp,10
lea    rbp,[rsp+10]
mov    rax,[rdi+8]        ; this.arr
mov    [rbp-8],rax
test   rax,rax            ; null-check
je     L01
cmp    dword ptr [rax+8],0 ; zero-length check
je     L01
cmp    dword ptr [rax+8],0
jbe    RNGCHKFAIL
add    rax,10             ; pula header
L00: movsxd rcx,dword ptr [rdi+10]  ; this.idx
     vmovsd xmm0,qword ptr [rax+rcx*8]
     ret
L01: xor eax,eax; jmp L00
RNGCHKFAIL: call CORINFO_HELP_RNGCHKFAIL; int 3
```

Estruturalmente idêntico ao assembly do Godbolt (mesmo `cmp`+`jae`
único no checked; mesmo null-check + zero-length-check duplicado +
aritmética de ponteiro no `fixed`) — as únicas diferenças são 2 `mov`
extras pra carregar campos de instância (`this.arr`/`this.idx`), porque
aqui os métodos são de instância (exigência do BenchmarkDotNet),
enquanto no Godbolt eram estáticos com array/índice como parâmetro.
Cross-validação confirma que os dois resultados (Godbolt e local) são
consistentes. Artigo companion passou a citar o BenchmarkDotNet como
fonte primária (validação de Release embutida e verificável) e o
Godbolt como cross-check secundário.

Timing do próprio BenchmarkDotNet (chamada isolada, array pequeno,
não é o cenário de acesso não-sequencial em array grande do
microbenchmark isolado acima, só complementa o assembly): `GetChecked`
1.57ns, `GetUnsafePerCall` 2.13ns (1.36x) — ordem de grandeza menor
que os outros dois experimentos porque aqui não há pressão de cache
(array de 1024 doubles, cabe em L1) nem chamadas repetidas
amplificando o custo de pin/unpin; o objetivo desse experimento
específico era o assembly, não reproduzir a magnitude da regressão.

## Contiguidade de `seq`: jagged → flat — adotado (2026-07-15)

### Contexto

`tSolution.seq` era `double[][]` (jagged) — cada linha (`seq[i]`) uma
alocação separada no heap, cada linha internamente achatada
(`j*3+k`). Categoria "parcial" em `TECNICAS_G1_G2.md` §1: contíguo
*dentro* da linha, não *entre* linhas — diferente de um array 1D
verdadeiro pra toda a estrutura.

### O que foi testado

`seq` trocado pra `double[]` único, tamanho `(dimen+1)*(dimen+1)*3`,
indexado via `i*(d+1)*3 + j*3 + k` — muda só `GetSeq`/`SetSeq` (a
camada de getter/setter já existente absorve a mudança, sem precisar
tocar nos ~40+ call sites em `GILS_RVND.cs`, ao contrário do Go).

### Microbenchmark (leve + intensiva)

- **Leve** (`burma14`): tempos de processo/JIT dominam, não medido
  separadamente (mesmo padrão de todos os outros testes desta sessão).
- **Intensiva** (`pr299`, N=5): **200.83s ± 12.78** (jagged) vs.
  **180.89s ± 13.22** (flat) — flat 9.9% mais rápido em média,
  distribuições com sobreposição leve (jagged mín. 192.99s, flat máx.
  203.33s) — sinal mais forte que o do Go, mas ainda não
  estatisticamente limpo (sem sobreposição) no N testado.
- **Corretude**: `COST: 6556628` idêntico nas duas variantes, e
  `COST: 20315` confirmado em `burma14` após aplicar em `main`.

### Decisão

**Adotado mesmo com sinal não-decisivo** — mesma decisão do Go (autor,
2026-07-15): aceitar o resultado direcionalmente favorável (~9.9%) sem
estender pra N=15. `tSolution.cs` de `main` simplificado pra usar só a
versão flat (sem manter a variante jagged condicional — decisão já
tomada, não vale manter os dois caminhos no código de produção).

**Consequência pro TCC**: os dados já coletados na campanha `main`
completa (90 linhas de C#, `mlp_testao/csharp.csv`) refletem a versão
*jagged* antiga — precisam ser recoletados pra refletir o código atual
(ver `CRONOGRAMA.md`).

## Pin único de `seq` via `GCHandle` — adotado (2026-08-03)

### Contexto

Motivado pela seção de bounds-checking do artigo companion
(`doc/artigo-linguagens/main.tex`, §5.3): o teste original de
bounds-check-off (`unsafe`/`fixed` por chamada, seção acima) mostrou
regressão, com a hipótese de que o overhead de pin/unpin por chamada
dominava a economia do bounds-check. Dois microbenchmarks isolados
(fora do solver real) validaram essa hipótese: pin único fora do laço
performa igual ou melhor que o indexador checked, só o pin por chamada
perde. `seq` desde 2026-07-15 é um array 1D flat único (não mais
jagged) — diferente da época do teste original, pinar `seq` inteiro
uma única vez agora cobre a estrutura toda, não só uma linha.

### O que foi testado

`GetSeq`/`SetSeq` trocados de indexador gerenciado pra ponteiro obtido
via `GCHandle.Alloc(seq, GCHandleType.Pinned)` **uma única vez**, no
construtor de `tSolution` — sem `fixed` nenhum, sem tocar nos ~40+ call
sites em `GILS_RVND.cs`. `solve()` cria só 3 instâncias de `tSolution`,
uma vez, reaproveitadas pro programa inteiro, então "pin único" aqui
significa pin pelo tempo de vida do processo (GCHandle não é liberado
explicitamente — só 3 objetos pinados, sem custo relevante).

```csharp
private GCHandle seqHandle;
private unsafe double* seqPtr;

public unsafe tSolution(int dimen, double c) {
    seq = new double [(dimen+1)*(dimen+1)*3];
    seqHandle = GCHandle.Alloc(seq, GCHandleType.Pinned);
    seqPtr = (double*)seqHandle.AddrOfPinnedObject();
    d = dimen; cost = c;
}

[MethodImpl(MethodImplOptions.AggressiveInlining)]
public unsafe double GetSeq(int i, int j, int k) {
    return seqPtr[i*(d+1)*3 + j*3 + k];
}
[MethodImpl(MethodImplOptions.AggressiveInlining)]
public unsafe void SetSeq(int i, int j, int k, double v) { seqPtr[i*(d+1)*3 + j*3 + k] = v; }
```

`<AllowUnsafeBlocks>true</AllowUnsafeBlocks>` adicionado ao `.csproj`.

### Microbenchmark isolado (fora do solver, aproximando o cenário real)

Array flat `(d+1)*(d+1)*3` com `d=299` (mesma instância do teste
abaixo), indexação idêntica à de `GetSeq`/`SetSeq`, stream de acesso
`(i,j,k)` uniformemente aleatório com ~15% de escritas (aproximando
leitura dominante na avaliação de vizinhança + escrita ocasional na
aplicação de movimento — não é um replay real do GILS-RVND). N=10,
`sink` idêntico em todas as execuções (correção confirmada):

| variante | ns/op (média ± dp) | vs. checked |
|---|---|---|
| checked (indexador gerenciado) | 18,24 ± 2,72 | 1,00x |
| unsafe/fixed por chamada | 24,15 ± 3,95 | 1,32x mais lento |
| unsafe, fixed no escopo do bench inteiro | 14,75 ± 2,93 | 0,81x (mais rápido) |

Margem menor que o teste isolado anterior (array totalmente aleatório,
sem escritas: 2,43x/0,89x) — provavelmente por causa do array menor
(~2,16MB, cabe melhor em cache) e da mistura de escritas.

### Microbenchmark no solver real (`pr299`, N=5 cada variante)

Mesma metodologia dos testes anteriores desta seção. CPU com ruído de
fundo (load average ~4 durante o teste, mesma máquina de sempre).
`COST: 6556628` idêntico nas 10 execuções (5 pinned + 5 baseline) —
corretude confirmada. `COST: 20315` confirmado em `burma14` antes do
teste intensivo.

| variante | n | média ± dp | min–max |
|---|---|---|---|
| pinned (GCHandle único) | 5 | 214,27s ± 14,66 | 201,95–239,42 |
| baseline (checked) | 5 | 223,25s ± 7,09 | 216,81–232,36 |

Pinned ~4,2% mais rápido em média, mas sinal fraco: a faixa do pinned
contém inteiramente a faixa do baseline (sem separação nenhuma entre
as distribuições) — mais fraco que o sinal da mudança jagged→flat
acima (9,9%, sobreposição leve) e muito mais fraco que a regressão
original do bounds-check-off por chamada (1,83x–1,96x, sem
sobreposição). Direção bate com os dois microbenchmarks isolados
(pin único ajuda), mas não é estatisticamente limpo nesse N.

### Decisão

**Adotado mesmo com sinal fraco** (autor, 2026-08-03) — mesmo
precedente da mudança jagged→flat: aceitar o resultado direcionalmente
favorável sem estender N. `tSolution.cs` de `main` passa a usar
ponteiro pinado único em vez do indexador gerenciado.

**Consequência pro TCC**: mais uma mudança em `tSolution.cs` desde a
última coleta completa da campanha de C# — dados de
`mlp_testao/csharp.csv` (`main`) precisam de recoleta por esse motivo
também, além do bug de Debug/Release já registrado (ver `TODO.md` da
raiz). Não afeta `standard` (só `main` é otimizada).

**Consequência pro artigo companion**: fortalece a explicação da
seção de bounds-checking — a hipótese de que o overhead é
especificamente o *pin por chamada*, não o `unsafe` em si, agora tem
validação em três níveis (assembly gerado, microbenchmark isolado, e
o solver real), embora o último com sinal mais fraco que os dois
primeiros.
