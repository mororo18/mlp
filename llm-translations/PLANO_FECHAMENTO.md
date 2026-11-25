# Plano: fechamento dos gaps das traduções (Seção 7.7)

Registrado em 2026-11-23, a pedido do autor. Continua o `PLANO.md` (a análise de
fatores nova, Fases 0 a 5) e substitui o escopo dele: agora entram **todos os
pares**, inclusive os que só tinham a análise antiga, e as traduções para C++
feitas de versões anteriores das nossas implementações (medidas como estão:
o autor decidiu não traduzir de novo). O `PLANO.md` continua
sendo o diário da análise anterior e a referência das decisões técnicas.
Problemas novos deste plano entram no fim da seção "Problemas e soluções"
deste arquivo.

## Objetivo

Para cada par de tradução, encontrar o menor conjunto de mudanças com nome que,
aplicadas **juntas**, levam a tradução ao tempo da nossa implementação
(Tabela `tab:translation-closing`, Seção `sec:LLM-closing` do artigo). Os fatores
medidos um por vez não se somam (Fortran → C++: ε 8,9%, matriz 23,6%, os dois
juntos −2,6%; JavaScript: o layout ajuda a tradução e quase não muda a nossa),
por isso a tabela usa combinações.

## Critérios

- **Fechado:** gap residual dentro de ±2% (CV das execuções em torno de 1%).
- **Só mudanças com nome:** as escolhas das Seções 4.3 e 5 (Tabela de desvios)
  ou um detalhe identificado no código. Trocar o código inteiro até a tradução
  virar a nossa não explica nada. O que não fechar entra como "não
  explicado", com o resíduo.
- **Direção:** quando a tradução é mais lenta, a nossa escolha vai na
  tradução; quando ela é mais rápida, a escolha dela vai na nossa (como no
  Java, Julia e Go). A direção oposta só confirma.
- **Mais de um caminho:** reportar o menor conjunto encontrado; se outro
  conjunto também fechar, registrar aqui.
- **Iteração:** depois das variantes de um fator e da combinação, se o
  resíduo passar de 2%, volta à leitura do código (Fase 2) no máximo duas
  vezes; depois disso, "não explicado".

## Decisões técnicas

Valem as decisões 1 a 9 do `PLANO.md` (base das nossas sem cópia, traduções
copiadas uma vez para `translations/`, variante = base + patches, build
herdado, `wall clock (s):` em volta da busca, custo conferido, `-n 10 -j 4`
em `taskset -c 8,9,10,11`, gap = média das razões, nomes `<par>.<fator>`),
com duas mudanças:

1. **Instâncias:** berlin52, kroA100, gr202 e lin318 para todos os pares,
   menos C++ → Lua e C++ → Python interpretados, que usam berlin52, kroA100 e
   lin105 (`instances_interp.txt`; justificativa no F16 do `PLANO.md`).
   Lua → C++ e Python → C++ medem um programa C++ e usam as 4. LuaJIT e PyPy
   usam as 4 (35 e 153 s no lin318).
2. **Combinações:** cada par termina com uma variante que aplica juntos todos
   os fatores com efeito acima de 2% (várias linhas de patch no manifesto,
   como o `fortran_cpp.epsilon_row_matrix`).
3. **Código no servidor:** enquanto o autor não liberar o push, os commits vão
   por `git bundle` para o clone do servidor (`~/mlp-test/fa/mlp`), como na
   `p26`. O `commit.txt` da rodada continua citando o hash local.

## Situação dos 24 pares (2026-11-25)

| Par | Gap | Situação | Frente |
|---|---|---|---|
| C++ → Fortran | −0,1 | fechado (nada) | — |
| C++ → Go | −2,8 → −1,4 | fechado (acessor na nossa) | — |
| C++ → Swift | −1,5 | fechado (nada) | — |
| C++ → Julia | −7,4 → −0,7 | fechado (índice base 0 na nossa, aquecimento) | — |
| Fortran → C++ | +34,1 → −2,6 | fechado (ε, matriz); fonte antiga | C |
| Go → C++ | +28,7 → +1,1 | fechado (ε, matriz); fonte antiga | C |
| Rust → C++ | +8,4 → −1,1 | fechado (matriz) | — |
| Swift → C++ | +12,4 → +1,5 | fechado (matriz); fonte antiga | C |
| C# → C++ | +12,5 → −0,2 | fechado (matriz) | — |
| C++ → JavaScript | +19,0 → +1,7 | fechado (tipo de array, layout intercalado e matriz *packed* na tradução; posições constantes e cópias locais na nossa) | A |
| C++ → Java | −10,6 → −1,8 | fechado (acessores, cópias locais, SoA e `BufferedReader` na nossa; B7) | A |
| C++ → C | +0,7 (+1,1 sem berlin52) | fechado (nada) | — |
| C++ → Rust | +69,5 → +1,1 | fechado (verificação de limites, R1) | B |
| C++ → C# | +8,4 → −3,5 | fechado (S1+S2; +0,2 sem berlin52, onde a tradução é 14% mais rápida também com aquecimento nos dois; causa não achada) | B |
| C++ → LuaJIT | +73,4 → −6,6 → +11,1 → +0,3 | fechado (L1-L5 na tradução, LJ1 na nossa; B9; gr202 1,125 em nota de rodapé) | B |
| C++ → PyPy | +39,6 → −4,1 → −1,4 | fechado (P1+P2 na tradução, PP1 na nossa; B7) | B |
| C++ → Lua | +89,2 → −14,0 → −1,4 | fechado (L1-L5 na tradução, LZ1+LZ2 na nossa; B9) | B |
| C++ → Python | +102,7 → −4,7 → −0,4 | fechado (tudo na tradução, PY1 na nossa; B7) | B |
| C → C++ | +2,3 → −1,3 | fechado (matriz, K2) | B |
| Lua → C++ | +1,9 → −0,5 | fechado (matriz, M2; ε não pesa) | B |
| Python → C++ | +6,0 → +0,7 | fechado (`seq` plano e matriz, Y1+Y2) | B |
| Java → C++ | +7,2 → +1,6 | fechado (`seq` plano e matriz, V1+V2); fonte antiga | B, C |
| JavaScript → C++ | +23,2 → 0,0 | fechado (deslocamentos e matriz, J1+J2); fonte antiga | B, C |
| Julia → C++ | +31,7 → −1,9 | fechado (deslocamentos e matriz, U1+U3); fonte antiga | B, C |

\* Análise antiga (outra máquina, berlin52/bier127/gr202, 1 a 10 execuções,
versões antigas das nossas): só ponto de partida.

Frente B depois das rodadas B4 (B4-comp em 4 instâncias, B4-interp em 3):
16 pares fechados (Java → C++ na borda, +2,0), 2 que passaram do ponto (a
tradução corrigida fica mais rápida que a nossa, fora de −2%: pela regra de
direção, falta medir a escolha dela na nossa) e 6 abertos.
Depois da B5-comp (2026-11-25): 20 fechados, 2 que passaram do ponto
(C++ → LuaJIT e C++ → Python) e 2 abertos (C++ → Lua, L4 na B6; C++ → Java,
A2). Os gaps de base dos cinco pares da B5 são os da própria rodada.
Depois da B6 e da B7 (2026-11-25): **22 fechados** (entram Java, Python e
PyPy, agora com a berlin52) e 2 abertos (C++ → LuaJIT, +11,1 com LJ1 na
nossa; C++ → Lua, B8). Artigo atualizado (`459010c` no artigo).

Gaps de base da Frente B: rodada `b-comp-20261124` (B1, 4 instâncias, 10
execuções, custo ok em todas). † Estimativa com a mediana do nosso C++ no
berlin52 (duas execuções fora da curva nesta rodada; ver Problemas); vale a
repetição do berlin52.

## Frente A. Pares abertos

- **A1. C++ → JavaScript.** Rodada `p26-20261123` (commit `d947dc0`), sobre o
  `Array`: a nossa, `three_arrays_array`, `const_offsets`, `local_info`,
  `local_info_consts`, a tradução com `Array` (`cpp_javascript.array`) e a
  tradução com `Array` intercalado (`cpp_javascript.interleaved_array`).
  Hipótese (P21, medida com `Float64Array`): tipo do array + layout +
  deslocamentos constantes e cópias locais explicam tudo. Falta: a
  combinação na tradução (`array` + `interleaved_array`) e, se a nossa com
  `local_info_consts` ficar mais rápida que a tradução, registrar que o par
  fecha na direção da nossa. Saída: linha da Tabela 13 e o parágrafo da
  7.6 sobre o `Array`.
  **Resultado (p26, 2026-11-23; custos iguais, CV ≤ 1,4%, sem bimodalidade).**
  Diferenças em relação à nossa sobre o `Array` (média das razões por
  instância, berlin52/kroA100/gr202/lin318):
  - tradução com `Array` (três vetores): +5,0% (reproduz o +4,9% da jsf64);
  - tradução com `Array` intercalado: −6,6% (o layout intercalado deixa a
    tradução 11% mais rápida);
  - nossa com três `Array`s: +1,4% (o intercalado é melhor na nossa);
  - nossa com `const_offsets`: −10,4%; com `local_info`: −0,8%; com os dois
    (`local_info_consts`): −10,7%;
  - mesmo tipo, mesmo layout e mesmas constantes e cópias locais nos dois
    programas (tradução `interleaved_array` × nossa `local_info_consts`): a
    tradução fica **+4,6%** mais lenta (+7,9, +4,5, +2,5, +3,7% por
    instância). Mesmo quadro do P21 com `Float64Array` (5 a 8%).
  **Nova hipótese (2026-11-24): matriz de distâncias "holey".** A tradução
  cria o vetor externo da matriz com `new Array(sz)` (elementos *holey* no
  V8: toda leitura `cost[s[k]]` confere buraco); o nosso usa `[]`
  preenchido em ordem (*packed*). No Node 20 local (kroA100, 2 execuções),
  só essa troca (X1, `packed_matrix.patch`) leva a tradução igualada de
  11,8/11,6 s para 11,1/10,7 s, empatando com a nossa igualada (11,1/10,8 s).
  Variantes `cpp_javascript.packed_matrix` e
  `cpp_javascript.interleaved_array_packed` na próxima rodada.
  Conclusão (p26): o par **não fecha**. Igualando as três escolhas, sobra uma
  diferença de 4,6% a favor da nossa, de causa não identificada (a
  tradução ainda difere em outros pontos do código). Para a Tabela 13:
  mudanças "tipo de array, layout intercalado (na tradução); posições
  constantes e cópias locais (na nossa)", resíduo +4,6%, Aberto.
  **Fechado na B5-comp (2026-11-25):** com a matriz *packed* (X1) junto do
  `Array` intercalado, a tradução fica +1,7% contra a nossa com
  `local_info_consts` (ver B5-comp na Frente B).
- **A2. C++ → Java.** Resíduo de −3,1% com acessores, cópias locais e
  armazenamento das subsequências na nossa (`ours_java.soa`), −9,3% com
  aquecimento. Hipóteses H-J5 e H-J6 do `PLANO.md`: o caminho quente do
  nosso Java tem mais classes e chamadas (`tArray`, `tSolution`, `tInfo`),
  o que atrasa a compilação e a inlining do C2. Passos: (1) contar e
  cronometrar as compilações (`-Xlog:jit+compilation`, ou
  `-XX:+PrintCompilation`) nos dois programas no berlin52 e no lin318; (2)
  comparar o inlining dos métodos quentes (`-XX:+PrintInlining`); (3)
  variante da nossa sem as classes intermediárias, se (1) e (2) apontarem
  para elas. Sem rodada grande até ter a hipótese.

- **A3. Traduções que ficam mais rápidas que a nossa (2026-11-25, local).**
  Pedido do autor: investigar C++ → LuaJIT (−6,6 com L1), C++ → Python
  (−4,7 com Q1 a Q4) e o JIT do C++ → Java (A2). Medições locais (notebook
  i5-8265U, um núcleo, ciclos de usuário com `perf stat`, porque o relógio
  oscila 10 a 15% aqui; Lua 5.1/LuaJIT 2.1, CPython 3.12, OpenJDK 21).
  - **LuaJIT.** Reproduz: kroA100, nossa 15,9 G ciclos (37,9 G instruções),
    tradução com L1 13,8 (30,8). Causa no código de máquina (`-jdump`): o
    traço do laço interno começa no meio da função, então `T`, `C` e `W`
    (locais copiadas de `info`) entram como valores de pilha (`SLOAD`), e
    cada acesso faz soma, conversão e deslocamento próprios; na tradução,
    as posições são literais (`+1`, `+2`, `+3`) e viram deslocamento fixo do
    endereço, com um só cálculo de base para T, C e W. Variantes da nossa:
    LJ1 posições literais nos índices: 13,9 (−12%); LJ1 + LJ2
    (`cost_concat_*` e `cost_new` declaradas dentro do laço, como na
    tradução, em vez de fora): 13,3 (−16%), abaixo da tradução. Só
    `local T = 1` (sem literal) não muda nada (16,1). Custo 143721/983128.
  - **Python.** Reproduz: berlin52, nossa 48,5 G ciclos, tradução
    corrigida 46,2. A nossa `update_subseq_info_matrix` cria `seq` e `s`
    locais, mas os laços continuam lendo `solut.seq[...]` e `solut.s[...]`:
    usar as locais (PY1) dá 46,7 (−3,7%); índices literais `[0] [1] [2]`
    sozinhos não mudam nada (48,3); os dois juntos (PY1 + PY2) dão 46,0,
    empate com a tradução.
  - **Java.** Com a `ours_java.soa` (acessores, cópias locais e
    armazenamento da tradução na nossa), o gr202 e o lin318 já empatavam
    (0,998 e 1,000 na `p11p12`); o resíduo é do berlin52 e do kroA100,
    execuções curtas. O `bench.py` fixa cada execução num núcleo, então as
    threads do C1/C2 dividem a CPU com a busca medida. Contagem com
    `-XX:+PrintCompilation` e o horário do início do relógio (kroA100): depois
    do início, o C2 da nossa ainda compila 31 métodos do JDK (`Scanner`,
    `java.util.regex`, `CharBuffer`, `String`), todos da leitura da
    instância, e 23 do programa; a tradução, 6 e 15 (lê com
    `BufferedReader` e `parseInt`). Tempo de compilação dentro do trecho
    medido (`CompilationMXBean`): nossa 1,1 s (berlin52) e 1,5 a 1,8 s
    (kroA100); tradução 0,6 a 0,7 s. Com `BufferedReader` na nossa (JV1),
    as compilações do JDK depois do início caem para 2; tempos locais
    (mínimo/mediana de 6): kroA100 `soa` 3,84/4,01, `soa` + JV1 3,50/3,67,
    tradução 3,81/4,08; berlin52 0,87/1,11, 0,82/1,17 e 0,75/0,97 (ainda
    sobra um pouco no berlin52). O `TIME` impresso pelo nosso Java é o
    tempo de CPU do processo (com as threads do JIT) e o da tradução, o da
    thread; o harness usa o `wall clock` dos dois, então a tabela não é
    afetada.
  - **Rodada B7** (disparada pela sessão local depois da B6, commits até
    `e004b08` pelo `b7.bundle`): `b7-comp-20261125` (4 instâncias: nossa
    LuaJIT, `cpp_lua.luajit_subseq_array`, `ours_lua.literal_offsets`,
    `ours_lua.literal_offsets_loop_locals`; nosso Java, `cpp_java.base`,
    `ours_java.soa`, `ours_java.bufread`, `ours_java.soa_bufread`) e
    `b7-interp-20261125` (3 instâncias: nosso Python, `cpp_python.all`,
    `ours_python.update_locals`, `ours_python.update_locals_literals`).
    Cerca de 2 h.
  - **Berlin52 do C# e do PyPy (2026-11-25, local).** No servidor, a
    diferença do C# é fixa (berlin52: nossa 1,46 s, tradução S1+S2 1,25 s;
    kroA100 5,00 e 4,97 s). Sem efeito na nossa (berlin52, local, 5
    execuções): `Stopwatch` no lugar do `DateTime.Now`, `InvariantGlobalization`
    e `ServerGarbageCollection` (os dois estão no `.csproj` da tradução e
    não no nosso). A compilação em camadas domina a berlin52 nos dois
    programas: com `DOTNET_TieredCompilation=0` (ou
    `DOTNET_TC_QuickJitForLoops=0`), nossa 1,8 → 0,5 s e tradução 1,45 →
    0,45 s; em ciclos do processo, a diferença cai de ~1,0 G para ~0,2 G.
    Com uma busca não medida antes da medida nos dois (aquecimento), a
    berlin52 empata: nossa 1,41 s, tradução 1,42 s (custos ok); a tradução
    quase não muda com o aquecimento e a nossa perde 0,4 s. Leitura
    provável: no tier 0 os acessores da nossa (`GetSeq`, `GetCost`, o
    indexador da `List<int>`) são chamadas de verdade até a recompilação.
    No PyPy, o tempo local oscila demais (kroA100 de 11,4 a 17,7 s na
    mesma variante) para decidir. Variantes com aquecimento nos dois
    programas (`ours_csharp.warm`, `cpp_csharp.unchecked_seq_flat_matrix_warm`,
    `ours_python.pypy_warm`, `cpp_pypy.float_w_interleaved_warm`) na
    `b7-short-20261125` (berlin52 e kroA100, `instances_short.txt`), entre
    a `b7-comp` e a `b7-interp`.
  - **PyPy no servidor (2026-11-25, a pedido do autor).** Diagnóstico
    fora do harness, no ambiente nix das rodadas (PyPy 7.3.20), container
    preso ao núcleo 14 (livre; a d66 ocupava 8 a 11), 5 a 6 execuções
    alternadas, CV abaixo de 2%, custos ok. Scripts, saídas e
    `PYPYLOG=jit-summary` em `diag/pypy-20261125/`. Programas: nossa
    (`main_pypy.py`), tradução P1+P2 (`float_w_interleaved`), as duas com
    uma busca não medida antes (aquecimento) e a nossa com o índice
    calculado em linha (PP1: `3*(i*(dimen+1)+j)+K` com K literal no lugar
    de `to_1D(i, j, info.X, info.dimen+1)`). Razões tradução/nossa
    (medianas):

    | | berlin52 | kroA100 |
    |---|---|---|
    | sem aquecimento | 0,886 / 0,893 | 0,969 / 0,968 |
    | com aquecimento nos dois | 0,973 / 0,972 | 0,983 / 0,980 |
    | nossa com PP1 | 0,946 | 0,987 |
    | nossa com PP1, com aquecimento nos dois | 0,973 | 0,989 |

    O aquecimento explica a maior parte da berlin52 (0,713 → 0,507 s na
    nossa, 0,631 → 0,493 s na tradução). O `jit-summary` mostra traços
    quase duas vezes maiores na nossa (operações gravadas na berlin52:
    64 mil contra 35 mil; tempo de tracing 0,097 contra 0,066 s): as
    chamadas a `to_1D` e as leituras de `info.T/C/W` e `info.dimen` entram
    nos traços. Com PP1, 46 mil operações e 0,079 s; sem aquecimento, PP1
    leva a berlin52 de 0,893 para 0,946 e o kroA100 de 0,968 para 0,987.
    Com aquecimento, PP1 só muda o kroA100; sobra cerca de 2,7% na
    berlin52 (cópias locais e outras leituras de campo ainda nos traços).
    Variante `ours_python.pypy_inline` na `b7-comp` (4 instâncias, contra
    `cpp_pypy.float_w_interleaved` e a nossa no PyPy), para a Tabela 13.
  - **Parâmetros do JIT do PyPy (2026-11-25, pedido do autor: induzir o
    aquecimento como no C#).** Grade no servidor (`diag3.sh`, `out3.txt`;
    núcleo 14, 3 execuções, custos ok), medianas em s, berlin52:

    | `--jit` | nossa | nossa PP1 | tradução | trad./nossa | trad./PP1 |
    |---|---|---|---|---|---|
    | padrão | 0,712 | 0,662 | 0,630 | 0,884 | 0,952 |
    | `threshold=200,function_threshold=300` | 0,830 | 0,690 | 0,655 | 0,790 | 0,950 |
    | `threshold=50,function_threshold=80` | 0,707 | 0,677 | 0,651 | 0,920 | 0,960 |
    | `trace_eagerness=50` | 0,692 | 0,646 | 0,615 | 0,888 | 0,951 |
    | `threshold=200,...,trace_eagerness=50` | 0,792 | 0,677 | 0,647 | 0,817 | 0,956 |
    | `trace_limit=20000` | 0,720 | 0,669 | 0,637 | 0,885 | 0,952 |

    No kroA100 todas as configurações ficam entre 0,966 e 0,991. Leitura:
    não há no PyPy um equivalente ao `DOTNET_TieredCompilation=0`; o JIT
    de traços só compila o que já executou, e baixar os limiares não
    antecipa o aquecimento (o melhor caso, `trace_eagerness=50`, tira 3%
    da nossa e 2% da tradução; `threshold=200` piora a nossa em 17%, com
    traços de caminhos ainda atípicos). Nenhuma configuração chega perto da
    busca não medida antes (nossa 0,507 s, tradução 0,493 s). A razão
    tradução/nossa com PP1 fica em 0,95 a 0,96 em todas: o resíduo sem
    aquecimento é do código (índice por função), não dos limiares.

## Frente B. Pares só com a análise antiga (13)

Fases como no `PLANO.md`:

- **B0. Protocolo.** Copiar as traduções de `claude_new/claude/` para
  `translations/` (`cpp_c`, `cpp_rust`, `cpp_csharp`, `cpp_lua`, `cpp_pypy`,
  `cpp_python`, `c_cpp`, `lua_cpp`, `python_cpp`, `java_cpp`,
  `javascript_cpp`, `julia_cpp`), com `RESULTS.md`,
  `CLAUDE.md` e `prompt.txt`. Conferir em cada uma o relógio, o trecho
  medido e a impressão de `wall clock (s):`; escrever os `measure.patch`.
  Conferir como a nossa LuaJIT e a nossa Lua rodam (mesmo código? FFI só no
  LuaJIT?) e qual programa é a referência de cada par (C++ → LuaJIT: a
  tradução para Lua rodando no LuaJIT contra a nossa Lua no LuaJIT).
  Catálogo e smoke de build de todas as bases, local e no servidor.
- **B1. Gap de base.** `<par>.base` e as nossas, nas instâncias do lote.
  Lote compilado feito (`runs/b-comp-20261124`, 2026-11-24, cerca de 1 h 30);
  razões por instância (berlin52, kroA100, gr202, lin318):

  | Par | Gap | Razões |
  |---|---|---|
  | C++ → C | +0,7 | 0,995 1,007 1,018 1,007 |
  | C++ → Rust | +69,5 | 1,726 1,699 1,718 1,636 |
  | C++ → C# | +8,5 | 0,883 1,095 1,195 1,165 |
  | C++ → LuaJIT | +73,4 | 1,567 1,651 1,957 1,758 |
  | C++ → PyPy | +39,3 | 1,321 1,410 1,440 1,400 |
  | C → C++ | +0,2 | 0,933 1,027 1,028 1,019 |
  | Lua → C++ | −0,1 | 0,912 1,024 1,029 1,033 |
  | Python → C++ | +3,8 | 0,922 1,045 1,089 1,096 |
  | Java → C++ | +5,1 | 0,941 1,066 1,102 1,094 |
  | JavaScript → C++ | +20,6 | 1,078 1,236 1,264 1,246 |
  | Julia → C++ | +29,3 | 1,192 1,348 1,346 1,286 |

  O berlin52 dos pares para C++ está puxado para baixo pelas duas execuções
  lentas do nosso C++ (ver Problemas). No C++ → C# o berlin52 é real: o
  nosso C# leva 1,46 s no berlin52 (a tradução, 1,29 s); a diferença do
  tempo fixo entra na B2. Lote interpretado (`b-interp-20261124`, berlin52,
  kroA100 e lin105, cerca de 1 h 50): C++ → Lua +90,3 (1,886 a 1,920) e
  C++ → Python +102,7 (2,010 a 2,044), custo ok.
- **B2. Mapa e lista de fatores**, a partir da Tabela de desvios, da
  análise antiga e da leitura do código. **Revisão do autor antes da B3.**
  Pontos de partida:

  | Par | Fatores da análise antiga e da tabela de desvios |
  |---|---|
  | C++ → C | tipo do índice (a tradução copiou o `unsigned` de 32 bits do nosso C++; gap antigo +4 → −1) |
  | C++ → Rust | verificação de limites (+81 → +3), matriz de distâncias aninhada (→ +39); combinação |
  | C++ → C# | verificação de limites (+6 → −2; o nosso usa ponteiro fixado com `GCHandle`), matriz aninhada (→ +2); combinação |
  | C++ → LuaJIT | armazenamento das subsequências (+97 → −19); chamadas e campos relidos (o trace compiler deve eliminar; conferir) |
  | C++ → Lua | chamadas de função para o índice (+75 → +32), campos relidos nos laços (→ +39), armazenamento das subsequências (→ +67); os três juntos levavam ao tempo da nossa |
  | C++ → Python | chamadas de função (+92 → +44), armazenamento das subsequências (→ +23), tipos dos elementos (`j + k` inteiro numa lista de floats), cópias locais; resíduo antigo de 8 a 9% atribuído às chamadas que ficam |
  | C++ → PyPy | tipos dos elementos (+29 → −8); armazenamento: as três listas da tradução são mais rápidas (com o nosso layout, +35), então vai na direção da nossa |
  | C → C++ | gap antigo 0: tipo do índice (com o nosso, +12) e matriz (−4) se compensam; conferir se o gap novo fica em ±2% sem mudança |
  | Lua → C++ | sem análise: matriz em `vector`s aninhados (Tabela de desvios), o resto pela leitura |
  | Python → C++ | armazenamento das subsequências (+6 → 0, um `vector` por linha), matriz aninhada |
  | Java → C++ | armazenamento das subsequências (+6 → +2, um `vector` por linha, como os arrays irregulares do nosso Java de julho), matriz aninhada |
  | JavaScript → C++ | deslocamentos de T, C e W lidos de campos (+25 → +1), matriz aninhada |
  | Julia → C++ | deslocamentos de T, C e W lidos de campos (+27 → +2), matriz plana (coluna) |

  **Leitura do código (2026-11-24), contra as nossas atuais.** Flags de
  compilação iguais nos dois lados (o harness impõe as do nosso `build.sh`).
  Gap de base da B1 entre parênteses.

  | Par | Tradução | Nossa | Fatores a medir (B3) |
  |---|---|---|---|
  | C++ → C (+0,7) | — | — | nenhum: fechado sem mudança |
  | C++ → Rust (+69,5) | `seq` plano com índice verificado; matriz `Vec<Vec<f64>>`; índices `i32` convertidos com `as usize` | `get_unchecked` nas duas matrizes, ambas planas (`--features flat`); `usize` | R1 verificação de limites; R2 matriz aninhada; R3 tipo do índice; combinação |
  | C++ → C# (+8,5; +15,2 sem berlin52) | `seq` plano de `struct` com índice verificado; matriz `double[][]` | `seq` por ponteiro fixado (`GCHandle`), sem verificação; matriz plana | S1 verificação de limites no `seq`; S2 matriz irregular; combinação. Berlin52: a tradução ganha (0,883); o nosso C# tem 1,46 s no berlin52 (aquecimento do JIT em camadas?), olhar à parte |
  | C++ → LuaJIT (+73,4) | uma tabela `{T, C, W}` por par (i, j); `to_1D` como função; `solut.seq`, `data.cost`, `solut.s` relidos a cada acesso | array `double[]` da FFI, índice em linha, campos em locais | L1 armazenamento das subsequências; L2 chamada do `to_1D`; L3 campos relidos; combinação |
  | C++ → PyPy (+39,3) | três listas (T, C, W); `seqW[idx] = j + k` inteiro numa lista de floats | uma lista intercalada; `float(j + k)` | P1 tipo dos elementos (`float(j + k)`); P2 layout (as três listas devem ajudar a tradução; com o nosso layout a análise antiga dava +35) |
  | C → C++ (+0,2; ≈+2,6) | `vector<tSeq>` plano; matriz `vector<vector<double>>`; índice `int` | ponteiro `tSeqInfo*`; `double**`; `to_1D` devolve `unsigned` | depende da repetição do berlin52; se ficar fora de ±2%: tipo do índice e matriz |
  | Lua → C++ (−0,1; ≈+2,3) | `seq` plano base 1 (`vector<double>`, intercalado); matriz aninhada | idem C → C++ | idem C → C++ |
  | Python → C++ (+3,8; ≈+6,2) | `vector<vector<array<double,3>>>` (um `vector` por linha); matriz aninhada | `seq` plano | Y1 armazenamento das subsequências; Y2 matriz aninhada |
  | Java → C++ (+5,1; ≈+7,5) | `seq2d`: um `vector<double>` por linha (`TEST = true`); matriz `vector<vector<double>>`; `tArray` para a solução | `seq` plano | V1 armazenamento das subsequências; V2 matriz aninhada; V3 `tArray` |
  | JavaScript → C++ (+20,6; ≈+23,4) | `seq` plano com deslocamentos `info.T/C/W` lidos de campos; matriz aninhada | deslocamentos constantes (campos da `struct`) | J1 deslocamentos constantes; J2 matriz aninhada |
  | Julia → C++ (+29,3; ≈+32,4) | `seq` plano (i, j, t) com deslocamentos `data.T/C/W` lidos de campos; solução base 1 (`S(i)` = `s[i - 1]`); matriz plana com `dimen + 1` por linha via método | deslocamentos constantes; base 0 | U1 deslocamentos constantes; U2 índice base 1; combinação |

- **B3. Variantes:** um patch por fator, a combinação, custo conferido no
  berlin52 antes de subir.
  Feita (2026-11-24): 29 variantes no manifesto (25 compiladas/LuaJIT/PyPy e
  4 do Lua interpretado), todas com 143721 no berlin52 pelo harness local.
  Ficaram de fora, por ora, R3 (índice `i32` no Rust), U2 (base 1 no Julia
  → C++) e V3 (`tArray` no Java → C++): entram se a combinação deixar
  resíduo. Na leitura apareceu um fator novo nos pares Java → C++ e Python →
  C++: o ε (`1e-16` e `1e-13` contra o `DBL_EPSILON` do nosso C++; no Python
  → C++ o custo da solução também desconta o ε), como no Fortran → C++ e no
  Go → C++ (V0 e Y0).
  C++ → Python (interpretado), lido contra o nosso `main.py`: a tradução
  guarda um objeto `SeqInfo` (`__slots__`) por par numa lista plana, chama
  `to_1D` em todo acesso, relê `solut.seq`, `data.cost` e `solut.s` e grava
  `W = j + k` inteiro; o nosso usa listas aninhadas `seq[i][j][k]`, cópias
  locais e `float(j + k)`. Variantes Q1 (`to_1D` em linha), Q2 (cópias
  locais), Q3 (listas aninhadas; tira também o `to_1D` dos laços), Q4
  (`float`) e a combinação. Lote B4-interp: as 4 do Lua e as 5 do Python,
  com as nossas e as duas bases (13 programas, cerca de 6 h).
- **B4. Rodadas:**
  B4-comp feita (`runs/b4-comp-20261124`, 1640 execuções, 41 programas,
  custo ok; o nosso C++ no berlin52 voltou a 0,092 a 0,094 s). Gaps (4
  instâncias):

  | Par | Base | Variantes | Leitura |
  |---|---|---|---|
  | C++ → Rust | +69,1 | R1 sem verificação +1,1; R2 matriz plana +46,0; R1+R2 +5,6 | fechado pela verificação de limites (R1); com a matriz plana junto fica pior |
  | C++ → C# | +8,4 | S1 `seq` por ponteiro 0,0; S2 matriz plana +5,6; S1+S2 −3,5 | fechado (S1); berlin52 ainda 0,86 a 0,89 em todas (tempo fixo do nosso C#) |
  | C++ → LuaJIT | +73,7 | L1 array −6,6; L2 `to_1D` +73,5; L3 locais +75,2; todos −11,2 | L1 explica tudo; com ele a tradução fica mais rápida que a nossa |
  | C++ → PyPy | +39,6 | P1 `float` +0,5; P2 intercalada +41,1; P1+P2 −4,1 | fechado (P1) |
  | C → C++ | +2,5 | — | berlin52 normal (1,020 a 1,037): um pouco fora de ±2% |
  | Lua → C++ | +2,2 | — | idem |
  | Python → C++ | +6,0 | Y0 ε +5,8; Y1 `seq` plano +2,4; Y2 matriz +3,7; todos +0,7 | fechado pela combinação (Y1 + Y2); ε não pesa |
  | Java → C++ | +7,2 | V0 ε +7,4; V1 `seq` plano +4,4; V2 matriz +4,2; todos +2,0 | na borda de ±2% (V1 + V2); ε não pesa |
  | JavaScript → C++ | +23,2 | J1 deslocamentos +5,5; J2 matriz +17,9; J1+J2 0,0 | fechado pela combinação |
  | Julia → C++ | +32,3 | U1 deslocamentos +14,7 | aberto: falta U2 (base 1) e o resto |

  - Lote B-comp (4 instâncias): C++ → C, Rust, C#, LuaJIT, PyPy; C, Lua,
    Python, Java, JavaScript, Julia → C++. Uns 11 pares x 4 a 6 programas.
    Estimativa: 3 a 4 h (o PyPy pesa: 153 s no lin318, cerca de 1 h
    sozinho).
  - Lote B-interp (`instances_interp.txt`): C++ → Lua e C++ → Python. Uns 7
    programas por par a cerca de 1,4 vez o tempo da nossa: cerca de 3 h por
    par, 6 h os dois (uma noite).
- **Resíduos depois da B4 (2026-11-24).** Julia → C++ (+14,7% com U1): a
  diferença não aparece com o GCC 14 local (nossa 2,41 s, tradução 2,31 s
  na kroA100) e aparece com os binários do GCC 15.3 do servidor rodados aqui
  (nossa 2,10 s, U1 2,52 s): o GCC 15 acelera o nosso C++ e não a tradução.
  No laço quente da reinserção, a tradução calcula `S(k) * (dimen + 1)`
  (`imul` + `cltq`) para a matriz plana; o nosso lê o ponteiro da linha
  (`double**`). Variantes U3 (`row_matrix`) e U1+U3. C → C++ (+2,5%) e
  Lua → C++ (+2,2%): matriz `double**` (K2, M2; a análise antiga já dava −4%
  com a matriz no C → C++) e, no Lua, o ε `1e-15` (M0). Java → C++ (+2,0%):
  o `tArray::get` é expandido no lugar (sem V3); a rodada mede de novo a base
  e a combinação.
- **B4-interp feita (2026-11-25):** `runs/b4-interp-20261124`, 390
  execuções, 13 programas, 3 instâncias (`instances_interp.txt`), custo ok.
  Gaps contra as nossas (Lua 5.4 e CPython):

  | Par | Base | Variantes | Leitura |
  |---|---|---|---|
  | C++ → Lua | +87,8 | L1 array +98,1; L2 `to_1D` em linha +43,0; L3 locais +50,3; todos +5,7 | L2 e L3 pesam e se somam (juntos rendem mais que a soma); L1 sozinho piora; sobra +5,7, fora de ±2% |
  | C++ → Python | +102,3 | Q1 `to_1D` em linha +55,7; Q2 locais +86,4; Q3 listas aninhadas +17,8; Q4 `float` +105,8; todos −4,7 | fechado pela combinação (Q3 explica a maior parte); `float` não pesa no CPython |

  Q1 no lin105 com CV de 11,4% (uma execução fora da curva); não muda a
  leitura. Resíduo do Lua (+5,7): comparando `lua_all` com a nossa, a
  tradução relê `s[i]`, `s[j]`, `s[i_prev]` e `s[j_next]` dentro dos laços
  internos e recalcula `dimen + 1` em cada índice; a nossa guarda tudo em
  locais. Variante L4 (`lua_all_hoist`, custo 143721 no berlin52). O teste
  local não decide (só há Lua 5.3 aqui, e nele a nossa nem fica atrás da
  `lua_all`); precisa do Lua 5.4 do servidor.
- **B5-comp feita (2026-11-25):** `runs/b5-comp-20261125` (commit
  `bd12ee0`), 760 execuções, 19 programas, 4 instâncias, custo ok, CV
  máximo 1,2%. Gaps contra as nossas (média das razões por instância,
  berlin52/kroA100/gr202/lin318):

  | Par | Base | Variantes | Leitura |
  |---|---|---|---|
  | Julia → C++ | +31,7 | U1 deslocamentos +14,3; U3 matriz `double**` +15,3; U1+U3 −1,9 (−3,5 −1,7 −2,8 +0,4) | fechado pela combinação; cada fator sozinho deixa cerca de 15% |
  | C → C++ | +2,3 | K2 matriz `double**` −1,3 | fechado (K2) |
  | Lua → C++ | +1,9 | M0 ε +1,9; M2 matriz −0,5; M0+M2 −0,5 | base já dentro de ±2% nesta rodada; a matriz leva a −0,5; ε não pesa |
  | Java → C++ | +7,2 | V0+V1+V2 +1,6 | fechado (reproduz a base da B4 e confirma a combinação) |
  | C++ → JavaScript | +19,0 | X1 matriz *packed* +12,9 (−5,1% sobre a base); intercalado −6,7; intercalado + X1 −9,3 | contra a nossa com `local_info_consts` (−10,8% sobre a nossa): intercalado +4,5 (reproduz a `p26`), intercalado + X1 **+1,7** (+3,0 +2,2 +1,0 +0,6): fechado |

  C++ → JavaScript fecha com as cinco escolhas: tipo de array (`Array`),
  layout intercalado e vetor externo da matriz *packed* na tradução;
  deslocamentos constantes e cópias locais na nossa (direção da nossa, que
  com elas fica mais rápida). O resíduo de 4,6% da `p26` era a matriz
  *holey*.
- **Artigo atualizado (2026-11-25, antes da d66 e da B6).** O `tables.py`
  passou a gerar as linhas da Frente B e da Tabela 13 (tabelas `pitfalls`,
  `other` e `closing` do `comparisons.csv`; Lua e Python em
  `berlin52/kroA100/lin105`, com o custo do lin105 tirado da `b4-interp`).
  No artigo: Tabelas 8 a 11 sem asteriscos (linhas da análise antiga
  trocadas pelas da B4/B5; saem C++ → C e C → C++ com o tipo do índice, que
  não têm mais gap a explicar), Tabela 13 completa, texto da 7.5 e da 7.6
  com os números novos, parágrafo da matriz *holey* no JavaScript e
  rascunho do texto da 7.7. Decisão a confirmar com o autor: C++ → C# e
  C++ → PyPy entram na Tabela 13 com a combinação (S1+S2: −3,5, +0,2 sem
  berlin52; P1+P2: −4,1, −1,7 sem berlin52) e não com o fator mínimo pela
  média (S1: 0,0, mas +4,5 sem berlin52; P1: +0,5, mas +3,9 sem berlin52),
  porque o fator mínimo só fecha na média, compensando o berlin52 (a
  tradução 11 a 14% mais rápida) com as três instâncias maiores. Falta:
  linha do C++ → Lua com a B6 e a Figura 7 com a d66.
- **Artigo, texto das causas (2026-11-25, `b6aeb6a` no artigo).** Na 7.6,
  o resíduo do C++ → Java vai para o JIT (leitura com `Scanner`, um núcleo
  por execução); na 7.7, as causas do LuaJIT (posições literais) e do
  Python (cópias locais não usadas na atualização das subsequências), o
  aquecimento na berlin52 do C# e do PyPy e o índice por função no PyPy.
  Números só das rodadas do protocolo: cada ponto que depende da B7 está
  com `\todo{B7: ...}`; a Tabela 13 não muda até a B7 (e a B6).
- **B6-lua feita (2026-11-25, 12h20 a 13h50):** `runs/b6-lua-20261125`
  (commit `6aad636`), 4 programas, 3 instâncias interpretadas, 10
  execuções, custo ok, CV máximo 2,8% (`cpp_lua.lua`). Gaps contra a nossa
  Lua 5.4: base +89,2; `lua_all` (L1+L2+L3) +5,6 (o mesmo da B4-interp);
  `lua_all_hoist` (+L4) **−14,0** (0,858 a 0,861 por instância). Com L4 a
  tradução passa do ponto, como no LuaJIT e no Python. Causa achada no
  código: a tradução usa posições a partir de 0 e calcula o índice da
  subsequência como `3*(i*n1+j)+K` (4 operações da VM); a nossa usa
  posições a partir de 1 e escreve `3*((s_size)*(i - 1) + (j - 1)) + T`
  (6 operações; o `luac` dá 99 a 263 instruções nas funções quentes da
  nossa contra 91 a 222 da tradução). Variante LZ1
  (`ours_lua.zero_based_index`): os índices sem os `- 1`, com o
  deslocamento somado nas posições (`T0 = T - 3*(dimen+1) - 3`, idem
  `C0`, `W0`). Local (Lua 5.3, berlin52, núcleo 2): nossa 33,4 s, LZ1
  28,2 s, `lua_all_hoist` 28,4 s, custos ok.
  - **Diagnóstico no servidor (2026-11-25, `diag/lua-20261125/`).**
    Ambiente nix das rodadas, container no núcleo 14 (a B7 ocupava 8 a
    11), execuções alternadas, custos ok. O `lua` do ambiente é o
    **5.2.4** (ver Problemas). Medianas, berlin52 (5 execuções) e kroA100
    (1 a 2): nossa 17,36 e 214,2 s; LZ1 13,91 e 172,6 s;
    `lua_all_hoist` 14,94 e 183,5 s. Razões `lua_all_hoist`/nossa 0,861
    e 0,857 (a B6 deu 0,858 a 0,861); LZ1/nossa 0,801 e 0,806;
    `lua_all_hoist`/LZ1 1,074 e 1,063: com LZ1 a nossa passa a ser mais
    rápida. Diferenças que sobram, no código: a tradução lê `s[j]` duas
    vezes no swap e no 2-opt, relê `s[j - 1]` no 2-opt e, na atualização
    das subsequências, relê do `seq` os valores da posição anterior; a
    nossa leva tudo em locais. Variante L5 (`cpp_lua.lua_all_hoist_carry`,
    `cpp_lua/carry_locals.patch`, sobre a L4): L5/LZ1 0,982 (berlin52) e
    0,968 (kroA100), cerca de −2,5%, na borda.
  - **Rodada B8** (disparada pela sessão local depois da B7, `b8.bundle`):
    `b8-lua-20261125`, 3 instâncias interpretadas, nossa Lua,
    `cpp_lua.lua_all_hoist`, `cpp_lua.lua_all_hoist_carry` e
    `ours_lua.zero_based_index`.
- **B7 feita (2026-11-25, 13h50 a 16h20; commit `441dcd0`).** Três lotes,
  custo ok em todos. Gaps (média das razões por instância; entre
  parênteses, sem berlin52), em `tables.py`/`comparisons.csv`:
  - `b7-comp-20261125` (4 instâncias, carga máxima 18,6): Java base −9,8
    (−10,6); tradução / nossa com SoA −3,8 (−2,5), com `BufferedReader`
    −7,4 (−7,6), com os dois **−1,8 (+0,8): fecha**. LuaJIT: L1 −6,6
    (o mesmo da B5); LJ1 deixa a nossa 15,9% mais rápida e a tradução
    fica **+11,1** (+12,9) contra ela (1,058, 1,073, 1,315 no gr202 e
    0,998); LJ1+LJ2 +8,2. PyPy: P1+P2 −3,9 (−1,9); PP1 na nossa: **−1,4
    (−0,8): fecha**, com berlin52 0,971.
  - `b7-short-20261125` (berlin52 e kroA100, carga máxima 8): células
    bimodais no berlin52 da nossa PyPy (0,71 a 1,18 s) e dos dois C#
    sem aquecimento (1,46 a 2,21 s e 1,24 a 2,03 s); as com aquecimento
    não. C# com aquecimento nos dois: −6,6 (berlin52 0,857, kroA100
    1,011): o aquecimento deixa os dois 20% mais rápidos no kroA100, mas
    **não tira a diferença da berlin52** (o diagnóstico no notebook dizia
    que tirava; no servidor, não). PyPy com aquecimento nos dois: −2,4
    (0,972 e 0,979).
  - `b7-interp-20261125` (3 instâncias, carga máxima 20): Python
    `cpp_python.all` −4,8; com PY1 na nossa **−0,4: fecha**; PY1+PY2
    −0,4 (as posições literais não mudam nada no CPython).
  - Artigo: Tabela 13 (Java, LuaJIT, PyPy, Python, Lua com a B6; nota d
    só para o C#), Tabela 12 (linhas do `BufferedReader`), 7.6 e 7.7
    reescritas, contagem 22 de 24.
- **B9-interp feita (2026-11-25, 18h15 a 19h20; commit `1724ac4`).**
  `runs/b9-interp-20261125`, Lua 5.2, 3 instâncias, 10 execuções, custo
  ok, carga máxima 17. Contra a nossa: L1-L5 −21,2; LZ1 −19,5; LZ1+LZ2
  −20,1. Tradução L1-L5 / nossa com LZ1 −2,1; **com LZ1+LZ2 −1,4** (0,986,
  0,976 e 0,996): **fecha**. Com isso, os 24 pares fecham. Artigo
  (`67b1469`): linha do Lua fechada, texto da 7.7 com a LZ2, contagem
  "all 24 pairs close", pendência da seção removida.
- **B9-comp feita (2026-11-25, 17h50 a 18h15; commit `1724ac4`).**
  `runs/b9-comp-20261125`, 4 instâncias, 10 execuções, custo ok, carga
  máxima 18,8 (CV de 20% só no C# padrão da berlin52, bimodal como na
  B7-short). C#: S1+S2 / nossa −3,5 (+1,0 sem berlin52; berlin52 0,832);
  com `QuickJitForLoops=false` nos dois **+2,5** (1,018 a 1,029 em todas;
  a berlin52 fica 4,9x mais rápida na nossa e 4,0x na tradução). LuaJIT:
  L1 / nossa LJ1 +10,7 (gr202 1,294); **L1-L5 / nossa LJ1 +0,3** (+1,3 sem
  berlin52; 0,973, 1,125 no gr202, 0,977 e 0,938): fecha pela média, com
  a nota de rodapé sobre o gr202 (decisão do autor).
  Artigo (`paper: C# tiered compilation...`): LuaJIT fechado com nota e,
  C# com a causa (camadas do JIT) na nota d e na 7.7, linha do Lua com a
  B8 (−2,2, aberto até a B9-interp), contagem 23 de 24.
- **B8 feita (2026-11-25, 16h23 a 17h50; commit `228e4bd`).**
  `runs/b8-lua-20261125`, Lua 5.2, 3 instâncias, 10 execuções, custo ok,
  CV máximo 1,2%, carga máxima 21. Contra a nossa: `lua_all_hoist` −13,6,
  L1-L5 (`lua_all_hoist_carry`) −21,5, LZ1 −19,7. Tradução / nossa com LZ1:
  L1-L4 +7,6; **L1-L5 −2,2** (−2,4 sem berlin52; 0,982, 0,980 e 0,972), na
  borda. Diferença que sobra, no código: a tradução com L4 tira dos laços
  internos `s[i_prev]` (swap) e `s[i_prev]` e `s[i]` (2-opt), que a nossa
  relê a cada iteração. Variante LZ2 (`ours_lua.zero_based_index_hoist`,
  `ours_lua/hoist_reads.patch`, sobre a LZ1), para a B9.
- **LuaJIT no gr202 (2026-11-25, `diag/luajit-20261125/`).** O LuaJIT
  do ambiente nix é o 2.1.1774638290 (o do notebook, 2.1.1703358377, não
  mostra a anomalia: tradução L1 / nossa LJ1 = 1,004 em ciclos). No
  servidor (núcleo 14, médias de 2), tempo / nossa com LJ1:

  | Tradução | gr202 | kroA100 |
  |---|---|---|
  | L1 | 1,30 a 1,31 | 1,07 |
  | L1 + matriz sem buracos | 1,27 | 1,07 |
  | L1 + reinserção com `for` | 1,27 | 1,05 |
  | L1-L4 | 1,17 | 1,02 |
  | L1-L5 | 1,13 | 0,97 a 0,98 |
  | L1-L5 + matriz sem buracos | 1,11 | 0,97 |
  | L1-L5 + reinserção com `for` | 1,13 | 0,97 |

  Perfil (`-jp`): mesma distribuição por função nos dois programas (busca
  de reinserção 52%, swap 22%, recálculo das subsequências 15%, 2-opt 9%),
  100% em código compilado, sem GC. A tradução é cerca de 13% mais lenta
  por igual no gr202. O `-jv` mostra 505 traços na tradução contra 219 na
  nossa, com muitas saídas laterais no traço da reinserção e 179 traços
  no quicksort recursivo da construção. Causa não achada. A matriz da
  tradução tem a diagonal nunca escrita (buraco em cada linha), mas isso
  pesa só 1 a 3%. Variante `cpp_lua.luajit_all_hoist_carry` (L1-L5 no
  LuaJIT) para a B9.
  **Decidido (autor, 2026-11-25):** ficam os números do LuaJIT do
  servidor (o do ambiente nix de todas as rodadas), com uma nota de rodapé
  no artigo sobre a anomalia do gr202 e a outra versão do LuaJIT.
- **C# na berlin52: causa achada (2026-11-25, `diag/csharp-20261125/`).**
  Ambiente nix das rodadas (.NET SDK 8.0.423), núcleo 14, binários da
  `b7-short`, execuções alternadas (mediana; berlin52 4 a 5, kroA100 2).
  O GC não pesa (`gcServer` 0/1, `gcConcurrent=0`: igual ao padrão). O JIT
  em camadas pesa tudo:

  | berlin52 | nossa | tradução | tr./nossa |
  |---|---|---|---|
  | padrão | 1,459 | 1,243 | 0,852 |
  | padrão, aquecimento nos dois | 1,446 | 1,235 | 0,854 |
  | `TC_QuickJitForLoops=0` | 0,332 | 0,340 | 1,024 |
  | `TieredCompilation=0` | 0,332 | 0,347 | 1,045 |
  | `TC_QuickJit=0` | 0,331 | 0,340 | 1,027 |
  | `TC_CallCountingDelayMs=0` | 0,430 | 0,456 | 1,060 |
  | `TC_OnStackReplacement=0` | 2,081 | 1,955 | 0,939 |
  | `TC_BackgroundWorkerTimeoutMs=0` | 1,459 | 1,245 | 0,853 |
  | `TieredPGO=0` | 1,398 | 1,195 | 0,855 |

  kroA100: padrão 5,00 e 4,93 s; `QuickJitForLoops=0` 3,99 e 4,07 (1,019).
  Leitura: os métodos quentes (com laços) rodam a execução inteira no
  código de tier 0 com OSR, por isso a busca de aquecimento não muda nada
  (o diagnóstico do notebook, outro .NET, dizia o contrário). A vantagem
  da tradução na berlin52 é a qualidade desse código não otimizado; com os
  laços compilados otimizados desde o início, a tradução fica cerca de 2%
  mais lenta, como nas instâncias maiores. Isso também vale para a Seção 6:
  o C# nas instâncias pequenas mede em boa parte o tier 0 (nossa berlin52
  1,46 s contra 0,33 s). Variantes CS1 (`TieredCompilationQuickJitForLoops`
  `false` no `.csproj`): `ours_csharp.quickjit_loops_off` e
  `cpp_csharp.unchecked_seq_flat_matrix_quickjit_loops_off`, para a rodada
  B9 (4 instâncias) depois da B8.
- **B5. Tabela:** linhas da Tabela 13 (`tables.py` passa a gerar também a
  tabela de fechamento) e saída dos asteriscos das Tabelas 8 a 12.

## Frente C. Traduções para C++ feitas de versões anteriores

**Decidido (autor, 2026-11-23): não traduzir de novo.** As traduções para C++
de Fortran, Go, Swift, Java, JavaScript e Julia partiram das
versões de julho das nossas implementações, que mudaram depois. Ficam como
foram entregues e são medidas contra o nosso C++ atual, que não mudou desde
julho. Consequências:

- **C1. Medição:** Java, JavaScript e Julia → C++ entram na frente B (B0 a
  B5). Fortran, Go e Swift → C++ já fecham na análise nova.
- **C2. Texto:** a nota "tradução de uma versão anterior da nossa
  implementação (julho de 2026)" fica nas Tabelas 8 a 13 para os sete
  pares, sem o "to be translated again". O texto da 7.7 explica que o
  fechamento compara com o nosso C++, e que as diferenças das fontes antigas
  (por exemplo, os deslocamentos de T, C e W lidos de campos no nosso
  JavaScript e Julia de julho) entram como fatores herdados da fonte.
- **C3. Figura 7 e durações:** o lado to-C++ continua válido (mesma tradução,
  mesmo C++ de referência), salvo o F3 (método de tempo). A
  `tab:translation-durations` não muda.
- **C4. C# → C++ (autor, 2026-11-23):** não leva a nota de versão anterior;
  o par entra na Tabela 13 sem marca.

## Frente D. Figura 7, lado from-C++

O tempo relativo dos from-C++ Fortran, Go, Swift, Java e JavaScript usa as
nossas versões antigas. Conferir se os tempos das traduções nas 66
instâncias são do mesmo servidor e do mesmo método das nossas atuais
(`results/main/cleaned_data_2026-11-23.csv`); se forem, só recalcular a
razão; se não, medir as 5 traduções nas 66 instâncias (cerca de 4 h).

**Conferido (2026-11-24):** a `results/translations/mlp_result/` (16 a 18/11)
tem as traduções e as nossas, mas com mais de 4 execuções por socket e antes
das mudanças do nosso JavaScript (flag, volta ao `Array`) e da repetição do
Java; a `main/cleaned_data` (17/11) também não segue o protocolo atual (só o
JavaScript e 4 instâncias foram refeitos). Com 10 execuções, as 5 traduções e
as nossas 5 nas 66 instâncias levariam cerca de 18 h. **Decidido (autor,
2026-11-24): 3 execuções**, traduções e nossas na mesma rodada
(`d66-20261125`, cerca de 5,5 h), disparada depois da B5-comp.

**Feita (2026-11-25, 09h44 a 15h10 no relógio do servidor).** As 10 séries
(5 traduções e as 5 nossas) com 3 execuções nas 66 instâncias, custos iguais
em todas, cv acima de 5% só nas instâncias de milissegundos (Fortran) e nas
pequenas do Java (JIT). Razão = média, nas 66 instâncias, de média da
tradução / média da nossa (tempo de parede da busca):

| C++ → | Antes (julho) | d66 | mín. a máx. por instância |
|---|---|---|---|
| Java | 0,928 | 0,889 | 0,64 a 1,18 |
| Go | 0,955 | 0,968 | 0,95 a 0,98 |
| Fortran | 0,986 | 0,992 | 0,67 a 1,16 |
| Swift | 0,950 | 0,995 | 0,95 a 1,04 |
| JavaScript | 1,147 | 1,162 | 0,97 a 1,25 |

Artigo: Figura 7 com esses valores (nova ordem: Java, Julia, Go, Fortran,
Swift), nota d da Tabela 13 aponta o aquecimento do JIT, e a 7.6 diz que
Fortran e Swift são mais rápidas por menos de 1%.

## Ordem proposta

1. A1 (dados da `p26`, já rodando).
2. B0 a B3 para os 13 pares (cerca de 2 dias de trabalho).
3. Lote B-comp numa tarde; lote B-interp numa noite.
4. A2 (Java) entre as rodadas, com as medições de JIT.
5. D (Figura 7, lado from-C++), numa noite, se precisar medir.
6. Tabela 13 completa, saída dos asteriscos, texto da 7.7 e o parágrafo de
   fechamento da Seção 7.

Pedido do autor (2026-11-25): quando B5-comp, d66 e B6-lua terminarem,
atualizar as tabelas e o texto do artigo com as novas descobertas (Tabela 13,
asteriscos das Tabelas 8 a 12, 7.7, Figura 7 com a d66) e registrar tudo
aqui. Fila no q-dual: B5-comp (até cerca de 08h), d66 (até cerca de 13h30),
B6-lua (até cerca de 14h30), cada uma disparada pela anterior a partir da
sessão local.

Servidor no total: cerca de 9 a 10 h para as variantes (B) e até 4 h para a
Figura 7 (D). O servidor está mais carregado desde 2026-11-23 (outros
dois usuários nos CPUs 0 a 6, 12 a 39 e 40 a 46, e um `podman` sem pinning que
já contaminou núcleos irmãos dos nossos): conferir o log de carga de toda
rodada e refazer as células com tempo bimodal por núcleo, como no P21.

## Pendências de decisão do autor

1. ~~C0: traduzir de novo três ou sete pares.~~ Decidido: não traduzir
   de novo (2026-11-23).
2. Push do `mlp-main` e do artigo (hoje só local; as rodadas vão por bundle).
3. Figura 7: medir de novo só as traduções afetadas ou todas com o método
   unificado (F3).

## Problemas e soluções

- **B0 feita (2026-11-24).** As 12 traduções copiadas de `claude_new/claude/`
  para `translations/` (sem binários, `.git` e `.claude`). Um
  `measure.patch` por tradução: a busca devolve o custo e o `main` imprime
  `COST:` depois de parar o relógio; relógio de parede em nanossegundos onde
  era em milissegundos (`cpp_c`, `cpp_rust`, `cpp_csharp`) ou só CPU
  (`c_cpp`, `lua_cpp`, `python_cpp`); a tradução para Lua usa
  `clock_gettime` pela FFI no LuaJIT (como a nossa) e `os.clock()` no `lua`.
  As 12 bases compilam e dão 143721 no berlin52 na máquina local.
- **Catálogo do Rust (2026-11-24).** O `catalog_base.py` decide o binário do
  `cpp_rust` olhando uma pasta que não existe aqui e caía em
  `cpp_to_rust`; agora também olha `translations/cpp_rust/Cargo.toml`
  (binário `solve`).
- **LuaJIT e Lua separados (2026-11-24).** A base `cpp_lua` gera duas
  execuções (lua e luajit). O lote compilado só quer o LuaJIT e o lote
  interpretado só o lua: variantes `cpp_lua.luajit` e `cpp_lua.lua`, cada
  uma com um `run_*_only.sh`; o `make_variants.py` passou a tratar o
  `run_*.sh` de um patch como o único comando da variante
  (`runtimes_extras=[]`).
- **Nosso C++ lento no berlin52 (2026-11-24).** Na `b-comp-20261124`, duas
  das 10 execuções do nosso C++ no berlin52 levaram 0,138 e 0,147 s (as
  outras, 0,093 a 0,094, como em todas as rodadas anteriores): cv de 19%,
  média 10% acima da mediana, e as razões das seis traduções para C++ no
  berlin52 caem para 0,91 a 0,94. É a única série com cv acima de 2,4%.
  Solução: repetir o berlin52 do nosso C++ e das seis traduções para C++ na
  mesma rodada (`INST=smoke.txt`, `REPS=10`) depois do lote interpretado, e
  usar essa rodada no berlin52.
- **Conexão (2026-11-24).** O jump parou de responder ao SSH (aceitava TCP
  sem mandar o banner) e voltou sozinho. O `~/.ssh/config` local passou a ter
  `ServerAliveInterval 30` nos dois hosts e `ControlPersist 12h` no q-dual;
  o acompanhamento das rodadas usa o master aberto. Antes de toda rodada,
  `git checkout -- c/solve` (o lote compilado recompila o binário
  versionado e o `run.sh` recusa a rodada seguinte).

- **Versão do Lua (2026-11-25).** O `lua` do `flake.nix` é o Lua **5.2.4**
  (`lua -v` dentro do ambiente nix das rodadas, diagnóstico de
  `diag/lua-20261125/`), não o 5.4 que as notas c, d e f das Tabelas 8, 10
  e 13 do artigo diziam. Corrigido no artigo para "Lua 5.2". O
  `versions.txt` das rodadas não registra `lua`, `luajit` nem `pypy3`;
  acrescentar ao `run.sh` depois da B8 (não mudar o harness com rodadas na
  fila). Decisão do autor: manter o 5.2 (o que todas as rodadas usaram) ou
  trocar para `lua5_4` e refazer as rodadas interpretadas do Lua.

- **Versões do ambiente nix (2026-11-25).** Consultadas no servidor, dentro
  do container das rodadas (uma sessão SSH): GCC 15.3.0, Go 1.26.5, rustc e
  cargo 1.97.0, Swift 5.10.1, .NET SDK 8.0.423 (runtime 8.0.29), OpenJDK
  21.0.12, Node.js 24.18.0, Julia 1.12.6, CPython 3.14.6, PyPy 7.3.20
  (Python 3.11.13), Lua 5.2.4, LuaJIT 2.1.1774638290; Nixpkgs `e2587ca`
  (2026-07-23) e, para Swift, `597283a` (2026-07-24). Artigo: Apêndice C
  (`tab:versions`) e frase na seção do ambiente, com `\todo` sobre as Figs.
  4 e 6 e os mapas de calor (rodada de 17/11, sem registro de versões).
  Monografia: Tabela 2 na Seção 3.2. O `run.sh` ainda não grava lua,
  luajit, pypy3 nem dotnet no `versions.txt`.
