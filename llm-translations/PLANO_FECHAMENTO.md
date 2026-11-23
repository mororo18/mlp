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
| C# → C++ | +12,5 → −0,2 | fechado (matriz); fonte antiga? | C |
| C++ → JavaScript | +18,9 → +4,9 | aberto | A |
| C++ → Java | −10,6 → −3,1 | aberto | A |
| C++ → C | +4* | a medir | B |
| C++ → Rust | +81* | a medir | B |
| C++ → C# | +6* | a medir | B |
| C++ → LuaJIT | +97* | a medir | B |
| C++ → PyPy | +29* | a medir | B |
| C++ → Lua | +75* | a medir (interpretado) | B |
| C++ → Python | +92* | a medir (interpretado) | B |
| C → C++ | 0* | a medir | B |
| Lua → C++ | sem análise | a medir | B |
| Python → C++ | +6* | a medir | B |
| Java → C++ | +6* | a medir; fonte antiga | B, C |
| JavaScript → C++ | +25* | a medir; fonte antiga | B, C |
| Julia → C++ | +27* | a medir; fonte antiga | B, C |

\* Análise antiga (outra máquina, berlin52/bier127/gr202, 1 a 10 execuções,
versões antigas das nossas): só ponto de partida.

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

- **B3. Variantes:** um patch por fator, a combinação, custo conferido no
  berlin52 antes de subir.
- **B4. Rodadas:**
  - Lote B-comp (4 instâncias): C++ → C, Rust, C#, LuaJIT, PyPy; C, Lua,
    Python, Java, JavaScript, Julia → C++. Uns 11 pares x 4 a 6 programas.
    Estimativa: 3 a 4 h (o PyPy pesa: 153 s no lin318, cerca de 1 h
    sozinho).
  - Lote B-interp (`instances_interp.txt`): C++ → Lua e C++ → Python. Uns 7
    programas por par a cerca de 1,4 vez o tempo da nossa: cerca de 3 h por
    par, 6 h os dois (uma noite).
- **B5. Tabela:** linhas da Tabela 13 (`tables.py` passa a gerar também a
  tabela de fechamento) e saída dos asteriscos das Tabelas 8 a 12.

## Frente C. Traduções para C++ feitas de versões anteriores

**Decidido (autor, 2026-11-23): não traduzir de novo.** As traduções para C++
de Fortran, Go, Swift, Java, JavaScript e Julia (e talvez C#) partiram das
versões de julho das nossas implementações, que mudaram depois. Ficam como
foram entregues e são medidas contra o nosso C++ atual, que não mudou desde
julho. Consequências:

- **C1. Medição:** Java, JavaScript e Julia → C++ entram na frente B (B0 a
  B5). Fortran, Go, Swift e C# → C++ já fecham na análise nova.
- **C2. Texto:** a nota "tradução de uma versão anterior da nossa
  implementação (julho de 2026)" fica nas Tabelas 8 a 13 para os sete
  pares, sem o "to be translated again". O texto da 7.7 explica que o
  fechamento compara com o nosso C++, e que as diferenças das fontes antigas
  (por exemplo, os deslocamentos de T, C e W lidos de campos no nosso
  JavaScript e Julia de julho) entram como fatores herdados da fonte.
- **C3. Figura 7 e durações:** o lado to-C++ continua válido (mesma tradução,
  mesmo C++ de referência), salvo o F3 (método de tempo). A
  `tab:translation-durations` não muda.
- **C4. Pendente:** conferir se o C# de julho já tinha o `GCHandle` (o
  escopo inicial do `PLANO.md` diz que não; o comentário do artigo sobre a
  Tabela 7 diz que o C# não foi afetado), para saber se o C# → C++ leva a
  nota.

## Frente D. Figura 7, lado from-C++

O tempo relativo dos from-C++ Fortran, Go, Swift, Java e JavaScript usa as
nossas versões antigas. Conferir se os tempos das traduções nas 66
instâncias são do mesmo servidor e do mesmo método das nossas atuais
(`results/main/cleaned_data_2026-11-23.csv`); se forem, só recalcular a
razão; se não, medir as 5 traduções nas 66 instâncias (cerca de 4 h).

## Ordem proposta

1. A1 (dados da `p26`, já rodando) e C4.
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

(vazio)
