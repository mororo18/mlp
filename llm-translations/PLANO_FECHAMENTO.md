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

## Situação dos 24 pares (2026-11-23)

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
| C++ → JavaScript | +18,9 → +4,9 | aberto | A |
| C++ → Java | −10,6 → −3,1 | aberto | A |
| C++ → C | +0,7 (+1,1 sem berlin52) | fechado (nada) | — |
| C++ → Rust | +69,5 | aberto | B |
| C++ → C# | +8,5 (+15,2 sem berlin52) | aberto | B |
| C++ → LuaJIT | +73,4 | aberto | B |
| C++ → PyPy | +39,3 | aberto | B |
| C++ → Lua | +90,3 | aberto (interpretado) | B |
| C++ → Python | +102,7 | aberto (interpretado) | B |
| C → C++ | +0,2 (+2,6†) | a confirmar (berlin52) | B |
| Lua → C++ | −0,1 (+2,3†) | a confirmar (berlin52) | B |
| Python → C++ | +3,8 (+6,2†) | aberto | B |
| Java → C++ | +5,1 (+7,5†) | aberto; fonte antiga | B, C |
| JavaScript → C++ | +20,6 (+23,4†) | aberto; fonte antiga | B, C |
| Julia → C++ | +29,3 (+32,4†) | aberto; fonte antiga | B, C |

\* Análise antiga (outra máquina, berlin52/bier127/gr202, 1 a 10 execuções,
versões antigas das nossas): só ponto de partida.

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

## Ordem proposta

1. A1 (dados da `p26`, já rodando).
2. B0 a B3 para os 13 pares (cerca de 2 dias de trabalho).
3. Lote B-comp numa tarde; lote B-interp numa noite.
4. A2 (Java) entre as rodadas, com as medições de JIT.
5. D (Figura 7, lado from-C++), numa noite, se precisar medir.
6. Tabela 13 completa, saída dos asteriscos, texto da 7.7 e o parágrafo de
   fechamento da Seção 7.

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
