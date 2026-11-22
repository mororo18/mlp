# Plano: nova análise de fatores das traduções por LLM

Registrado em 2026-11-20. Este arquivo é o plano e também o diário do
experimento: todo problema que aparecer entra na seção "Problemas e soluções",
com a solução adotada.

## Objetivo

Refazer a análise de fatores das traduções por LLM (Seção 7 do artigo) contra
as implementações atuais do `mlp-main`, nas traduções dos dois sentidos, e
reconstruir com ela as tabelas de gaps (hoje Tabelas 8 a 12) e o texto das
Seções 7.5 e 7.6. A análise antiga foi feita sobre versões anteriores das
nossas implementações, em outra máquina, sem dados brutos nem scripts
(`doc/artigo-linguagens/claude/PROBLEMAS_METODOLOGICOS_TRADUCOES.md`, itens 4,
11, 12, 13, 24 e 26).

## Escopo atual (autor, 2026-11-20, substitui o escopo inicial abaixo)

A análise fica restrita aos casos em que a análise antiga não serve. Um par
entra se (1) a tradução ganha da nossa implementação na Figura 7, (2) a nossa
implementação mudou depois da análise antiga, ou (3) a análise antiga contradiz
a Figura 7 (sinal oposto, ou gap grande que ela não explica).

- **Variantes completas (Fases 2 a 4):** C++ → Java, Julia, Swift, Go e
  Fortran (critério 1; tempos relativos de 0,93 a 0,99), C++ → JavaScript
  (critério 2; a análise antiga foi feita sobre o nosso JavaScript anterior e
  não entra nas tabelas) e Rust → C++ (critério 3; 1,104 na Figura 7 contra
  0,93 na análise antiga).
- **Só gap de base (Fase 1):** C# → C++ (1,132 contra 1,00 a 1,02), Swift →
  C++ (1,147 contra 1,02 a 1,03), Fortran → C++ (1,30 contra 1,13 a 1,17) e
  Go → C++ (1,23 contra 1,11 a 1,15). Se o gap medido no mesmo experimento
  bater com a explicação antiga, a diferença era da máquina e da medição e
  eles não precisam de variantes. Se não bater, entram nas Fases 2 a 4.
- **Fora:** os outros pares, cuja análise antiga continua coerente com a
  Figura 7.
- **Sentido das variantes:** para os pares em que a tradução ganha, o sentido
  principal é o direto (a escolha da tradução no nosso código), e o inverso
  confirma. Nos outros, como antes.
- **Instâncias:** as mesmas quatro (berlin52, kroA100, gr202, lin318), 10
  execuções; não há linguagem interpretada no escopo atual.
- **ETA:** cerca de 12 a 16 h de trabalho e 4,5 a 5 h de servidor, em 2 a 3
  sessões.

## Escopo inicial (autor, 2026-11-20, substituído acima)

As decisões abaixo continuam valendo no que não contradiz o escopo atual
(correções só de medição, patches pelo assistente, local, push antes da Fase 1).

- **Traduções:** todas, nos dois sentidos (12 a partir de C++, mais LuaJIT
  sobre a tradução para Lua, e 11 para C++).
- **Gap de base:** a Figura 7 do artigo **não** é medida de novo. O gap de cada
  tradução é medido dentro deste experimento, nas mesmas instâncias das
  variantes.
- **Instâncias das variantes:** berlin52, kroA100, gr202 e lin318 (uma por
  faixa de tamanho), 10 execuções. Lua e Python interpretados: berlin52 e
  st70.
- **Correções nas traduções:** só de medição (relógio de parede em volta da
  busca, impressão fora do trecho medido), cada uma como patch e declarada no
  texto. O código da busca não muda.
- **Sentido das variantes:** a mudança inversa (a nossa escolha na tradução)
  para todos os fatores. A direta (a escolha da tradução no nosso código) só
  onde ela confirma um fator ou mede um defeito nosso (por exemplo, o
  `unsigned` do nosso C++).
- **Traduções para C++ feitas de fontes antigas** (Fortran, Go, Java,
  JavaScript, Swift, Julia, e C# antes do `GCHandle`): medidas como estão,
  contra o nosso C++ atual. O texto diz que partiram de versões anteriores.
  Retraduzir com a LLM fica fora deste experimento.
- **Patches:** escritos pelo assistente. O autor revisa a lista de fatores da
  Fase 2 antes das rodadas.
- **Local:** esta pasta, no `mlp-main`, que é o repositório público citado pelo
  artigo. O repositório do artigo recebe só as tabelas geradas
  (`results/factor_analysis/`), citando o commit deste repositório.
- **Push do `mlp-main`:** antes da Fase 1, para a rodada citar um commit
  público (hoje `852fe0c`, local).

## Estrutura desta pasta

```
llm-translations/
  PLANO.md               este arquivo (plano e diário)
  README.md              protocolo, versões e como reproduzir (na Fase 5)
  translations/<par>/    fontes como a LLM entregou (sem binários), RESULTS.md, CLAUDE.md
  prompt.txt             prompt inicial das traduções
  harness/               bench.py e core.py (cópia, com a origem registrada) + catalog.py
  variants/
    manifest.csv         uma linha por variante: id, base, patches, nota
    <par>/measure.patch  só medição
    <par>/<fator>.patch  mudança inversa (a nossa escolha na tradução)
    ours/<ling>/<fator>.patch   mudança direta (a escolha da tradução no nosso código)
  make_variants.py       monta cada variante (base + patches) e o catálogo do bench.py
  run.sh                 comando exato no servidor + log de carga + versões
  instances.txt
  runs/<data>/           CSVs brutos, log de carga, versões dos compiladores
  analyze.py             tabelas de gaps (CSV e LaTeX) a partir dos CSVs brutos
```

## Decisões técnicas fixadas

1. **Base das nossas implementações:** as pastas de linguagem deste
   repositório, no commit da rodada, sem cópia. O `run.sh` grava o hash e se
   recusa a rodar com mudanças não commitadas nessas pastas.
2. **Base das traduções:** copiadas uma vez para `translations/`, como
   entregues (de `claude_new/claude/`). A tradução de base do experimento é a
   entregue mais o `measure.patch`.
3. **Variante = base + lista ordenada de patches**, nunca cópia editada à mão.
   Cada patch faz uma única mudança e é gerado com `diff -u` a partir de uma
   cópia editada. Variantes combinadas são várias linhas de patch no
   manifesto, sem arquivo de patch próprio.
4. **Build e execução:** cada variante herda a configuração da sua base, com as
   flags do nosso `build.sh`, reaproveitando o catálogo do harness antigo
   (`benchmark_claude_3instancias.py`: `REFERENCIAS`, `PROJETOS`,
   `TOOLCHAINS`). Nenhuma variante tem script de build próprio.
5. **Medição:** todo programa imprime `wall clock (s):` em volta da busca, e o
   `bench.py` grava esse valor (ele usa o `wall clock` quando existe).
6. **Correção:** toda execução tem de reproduzir o custo esperado da instância;
   o `analyze.py` falha se alguma não reproduzir.
7. **Rodadas:** um `bench.py -t` por lote, com todas as variantes do lote
   intercaladas por instância (o `bench.py` põe as instâncias no laço de fora),
   `-n 10 -j 4`, `taskset -c 8,9,10,11` (duas execuções por socket, núcleos
   físicos distintos), Julia 1.12.6 e as demais toolchains do servidor. Log de
   carga a cada minuto, como nas rodadas do Julia.
   - Lote 1: compiladas e JIT, nas quatro instâncias.
   - Lote 2: Lua e Python interpretados, em berlin52 e st70.
8. **Análise:** os CSVs brutos nunca são editados. Gap = 100 x (média, nas
   instâncias, da razão variante / nossa implementação - 1), como nas tabelas
   atuais; o script também dá o intervalo por instância.
9. **Nomes:** `<par>.base`, `<par>.<fator>` (inversa) e `ours_<ling>.base`,
   `ours_<ling>.<fator>` (direta). `<par>` é o nome da pasta da tradução
   (`cpp_rust`, `julia_cpp`, ...).

## Fases

### Fase 0. Protocolo
- [x] Push do `mlp-main` e registro do commit da rodada (`c205463`).
- [x] Copiar as traduções (só fontes) para `translations/`, com `RESULTS.md`,
      `CLAUDE.md` e `prompt.txt`.
- [x] Copiar `bench.py` e `core.py` do servidor para `harness/`, registrando a
      origem (`script.zip` do coautor, sem alterações).
- [x] Conferir em cada uma das traduções do escopo (11): relógio usado, o que fica dentro
      do trecho medido, impressão de `wall clock (s):`. Escrever os
      `measure.patch` necessários (o do C++ → Julia já existe:
      `doc/artigo-linguagens/results/julia/cpp_julia-2026-11-20/main.jl.patch`).
- [x] Conferir que as nossas implementações imprimem `wall clock (s):` só em
      volta da busca.
- [x] `catalog.py` (como `harness/catalog_base.py` + catálogo gerado) e
      `make_variants.py`; teste de montagem e build de todas as bases (local e
      no servidor, `runs/smoke*`).

### Fase 1. Gap de base
- [x] Rodar `<par>.base` e as nossas implementações nas instâncias do lote,
      no servidor (`runs/base-20261120/`, commit `9e5c200`).
- [x] Gap de cada tradução; lista dos gaps que ainda existem (abaixo).

#### Fase 1: resultado (rodada de 2026-11-20, 22:46 a 00:05 UTC)

720 execuções (18 programas x 4 instâncias x 10), nenhum timeout, todos os
custos iguais aos da nossa implementação. Carga: os 4 processos da rodada e
os 7 do outro usuário em quase todas as leituras (67 de 80 com 11 processos
acima de 50%, 6 com 12). cv máximo por instância: 2,8% (Java), demais
abaixo de 1,1%. Gap = 100 x (média das razões nas 4 instâncias - 1):

| Par | Gap | Razões (berlin52, kroA100, gr202, lin318) | Figura 7 |
|---|---|---|---|
| C++ → Fortran | −0,1% | 1,005 0,999 0,994 1,000 | 0,986 |
| C++ → Go | −2,8% | 0,975 0,969 0,969 0,974 | 0,955 |
| C++ → Java | −11,1% | 0,872 0,805 0,951 0,930 | 0,928 |
| C++ → JavaScript | +12,7% | 1,140 1,195 1,184 0,990 | 1,147 |
| C++ → Julia | −7,7% (−3,4% sem berlin52) | 0,796 0,966 0,961 0,970 | 0,942 |
| C++ → Swift | −1,4% | 0,991 0,991 0,974 0,989 | 0,950 |
| Rust → C++ | +8,4% | 1,071 1,089 1,093 1,082 | 1,104 |
| C# → C++ | +12,3% | 1,140 1,134 1,132 1,087 | 1,132 |
| Swift → C++ | +12,2% | 1,129 1,133 1,137 1,089 | 1,147 |
| Fortran → C++ | +33,6% | 1,363 1,346 1,341 1,295 | 1,303 |
| Go → C++ | +28,8% | 1,257 1,279 1,318 1,299 | 1,232 |

Leitura:
- **Fortran:** sem gap; sai da análise (F1 e F2 descartados).
- **Go e Swift (from-C++):** a tradução ainda ganha, por pouco e de forma
  consistente (cv abaixo de 1%): 3% no Go, 1 a 2,6% no Swift.
- **Java e Julia:** a tradução ganha. No Julia, a diferença maior está no
  berlin52 (0,796), o que reforça o L1 (compilação dentro do tempo).
- **JavaScript:** a tradução perde de 14 a 20% até o gr202 e empata no lin318.
  No lin318 os dois lados levam cerca de 333 s, 6,5 vezes o gr202 (no C++ a
  razão é 3,1): um efeito do V8 nessa instância, a investigar.
- **Rust → C++:** a tradução perde 7 a 9%, como na Figura 7 (1,104), e não
  ganha como na análise antiga (0,93).
- **As quatro traduções to-C++ só de base não batem com a explicação
  antiga:** C# 12% (antes ~0), Swift 12% (antes 2 a 3%, passagem por
  valor), Fortran 34% e Go 29% (antes 13 a 17% e 11 a 15%, o ε). Pelo
  critério do escopo, passam a precisar de variantes (decisão do autor).

### Fase 2. Mapa das diferenças
- [x] Para cada par, listar as diferenças de projeto entre a tradução e a
      nossa implementação atual (atualiza a Tabela 7 e o mapa do apêndice). As
      seis linguagens que mudaram precisam ser feitas do zero.
- [x] Lista de fatores a testar (seção "Fase 2: lista de fatores proposta"): Tabela 7, candidatos dos relatórios antigos
      (`doc/tcc/relatorios-ia/`) e o que aparecer de novo. **Revisão do
      autor antes da Fase 3.**

### Fase 3. Variantes
- [x] Um patch por fator; manifesto; montagem e build de todas as variantes;
      conferência do custo em uma instância pequena antes de subir para o
      servidor.

### Fase 4. Rodadas
- [x] Lotes 1 e 2 no servidor. Estimativa: umas 60 variantes compiladas ou
      JIT, de 6 a 8 horas; Lua e Python interpretados, umas 2 horas.

### Fase 5. Análise e texto
- [ ] `analyze.py`, tabelas em `results/factor_analysis/<data>/` do artigo,
      texto da 7.5 e da 7.6, `README.md` desta pasta.

## Fase 5: plano das tabelas e da análise (2026-11-22, para revisão do autor)

### Situação de hoje no artigo (Seção 7, `main.tex` em `186064e`)

| Tabela | Rótulo | Origem dos números | Situação |
|---|---|---|---|
| Desvios do catálogo (7.4) | `tab:translation-deviations` | leitura do código | conferir célula a célula com o mapa da Fase 2 |
| Gaps das armadilhas, from/to (7.5) | `tab:translation-gaps-pitfalls-from/-to` | relatórios antigos (`gaps_pitfalls.csv`): outra máquina, berlin52/bier127/gr202, 1 a 10 execuções, versões antigas das nossas | refazer |
| Gaps dos outros fatores, from/to (7.6) | `tab:translation-gaps-other-from/-to` | relatórios antigos (`gaps_other.csv`) | refazer |
| Traduções mais rápidas (7.6) | `tab:translation-faster` | relatórios antigos + coluna de situação | substituir (ver decisão 3) |

Dependências fora desta fase: a Figura 7 (`fig:translation-runtime`) ainda é
da campanha antiga (pendência "unificar o método de tempo da Figura 7"), e o
JavaScript dela tem o efeito do Node 24 (P16 a P18).

### O que a análise nova cobre

Mesmo protocolo para tudo: servidor, berlin52, kroA100, gr202 e lin318, 10
execuções, comparação sempre dentro da mesma rodada; gap = 100 x (média das
razões nas 4 instâncias − 1), o método da Figura 7, mais o gap sem o berlin52
nas linguagens com JIT.

- **From-C++ (a tradução ganha; variante direta no nosso, inversa confirma):**
  Java, Julia, Swift, Go; JavaScript (a tradução perde; critério 2). Fortran
  sai (−0,1%, sem gap).
- **To-C++:** Rust, C#, Swift, Go, Fortran.
- **Fora do escopo (só análise antiga):** from-C++ Rust, C#, C, Lua, LuaJIT,
  Python, PyPy; to-C++ C, Python, Java, JavaScript, Julia.

### Esboço das tabelas novas (números provisórios, do diário acima)

Armadilhas do catálogo (7.5), to-C++ (gap sobre o nosso C++):

| Par | Gap | Armadilha | Gap com a nossa escolha |
|---|---|---|---|
| Rust → C++ | +8,4% | matriz de distâncias | −1,1% |
| C# → C++ | +12,5% | matriz de distâncias | −0,2% |
| Swift → C++ | +12,4% | matriz de distâncias | +1,5% |
| | | passagem por valor | +11,4% |
| Go → C++ | +28,7% | matriz de distâncias | +25,7% |
| Fortran → C++ | +34,1% | matriz de distâncias | +23,6% |

Armadilhas do catálogo (7.5), from-C++ (gap da tradução sobre o nosso; a
mudança é a escolha da tradução no nosso código quando ela ganha):

| Par | Gap | Armadilha | Gap depois da mudança |
|---|---|---|---|
| C++ → JavaScript | +12,8% | representação das subsequências (`Float64Array`) | +0,9% (I) |
| C++ → Swift | −1,5% | cópia da solução na perturbação | ≈ 0% (D) |
| C++ → Java | −10,6% | subsequências em três vetores | entra na combinação abaixo |
| C++ → Julia | −7,1% | compilação JIT no tempo medido (berlin52) | −3,3% (D: os dois com aquecimento) |

Outros fatores (7.6):

| Par | Gap | Fator | Gap depois da mudança |
|---|---|---|---|
| Go → C++ | +28,7% | ε na atualização do melhor movimento | +2,2% |
| | | ε + matriz | +1,1% |
| Fortran → C++ | +34,1% | ε | +8,9% |
| | | ε + matriz | −2,6% |
| Rust, C#, Swift, Go, Fortran → C++ | | índice `unsigned` de 32 bits (o nosso) | +44,0, +27,8, +13,5, +33,6, +47,3% (piora) |
| C++ → Java | −10,6% | cópias locais da rota (J2), acessores (J1) | a calcular no `analyze.py` (o nosso com J1 + J2 fica 6,6% mais rápido) |
| | | J1 + J2 + três vetores + aquecimento | +0,2% sem o berlin52 |
| C++ → Julia | −3,8% sem berlin52 | índice de cidade a partir de 1 | −0,4% (com aquecimento) |
| C++ → Go | −2,8% | acessor `at()` (G2) | −1,4% (resto: geração de código) |

E o controle do índice no nosso C++ (`size_t_counter`, +32,3%) e o P14 (só
subsequências +16,1%, com índices não negativos +1,1%) para o mecanismo.

### Decisões do autor

1. **Pares fora do escopo (análise antiga).** (a) Linhas antigas nas mesmas
   tabelas, marcadas com nota ("análise anterior: outra máquina, berlin52,
   bier127, gr202"), e o parágrafo de método descreve os dois experimentos;
   (b) só a análise nova nas tabelas e os achados antigos no texto, com a
   ressalva (perde a linha mais forte da seção, verificação de limites no
   C++ → Rust, 81% → 3%); (c) medir de novo esses pares (escrever as
   variantes: cerca de um dia de trabalho e 3 a 4 h de servidor).
   Recomendação: (a) agora, (c) se der tempo antes da submissão.
2. **Divisão 7.5 × 7.6.** Manter "armadilhas do catálogo" × "outros
   fatores". Classificação proposta: catálogo = matriz de distâncias,
   armazenamento das subsequências (inclui `Float64Array` e os três vetores
   do Java), passagem por valor/cópia, verificação de limites, tipos
   estáveis, compilação JIT no tempo medido (Seção 5); outros = ε, índice
   `unsigned`, acessores, cópias locais, índice a partir de 1, geração de
   código.
3. **Tabela das traduções mais rápidas.** A coluna de situação existia
   porque a análise antiga era sobre versões velhas das nossas; a nova é
   sobre as atuais. Proposta: tirar a tabela e pôr esses pares nas tabelas
   de gaps (linhas from-C++ com a mudança direta), com um parágrafo curto na
   7.6 sobre as traduções que ganham.
4. **Colunas.** Um só par de colunas, "Gap" e "Gap depois da mudança", com
   (D) ou (I) dizendo em qual programa entrou a mudança. Alternativa: duas
   tabelas separadas para os pares em que a tradução ganha.
5. **JavaScript no lin318.** Os dois lados do par estão no patamar lento do
   Node 24 (empate em cerca de 333 s). Proposta: medir o par e o
   `cpp_javascript.array` com `--no-concurrent-recompilation` nas 4
   instâncias (cerca de 2 h de servidor, depois das rodadas do P19) e usar
   esses números, com nota.

### Decisões do autor (2026-11-22)

1. Pares fora do escopo: **(a)**, linhas antigas nas tabelas com nota; (c)
   fica como pendência futura (F1).
2. Divisão 7.5 × 7.6 mantida, com a classificação proposta. O aquecimento
   do JIT (Java, Julia) não é escolha de nenhum dos dois programas: vai em
   linha própria, com o valor medido com os dois aquecidos e nota.
3. Tabela das traduções mais rápidas sai; os pares vão para as tabelas de
   gaps, com um parágrafo curto na 7.6.
   **Revisto pelo autor (2026-11-22):** a tabela fica, no formato das de
   gaps (duas colunas de mudança, aquecimento em linha própria, sem o
   berlin52 entre parênteses para Java e Julia). Feito no artigo (`8aec4d8`),
   com os números do `tables.py`; o texto da 7.6 sobre ela ainda é o antigo.
4. **Colunas: duas, sem marcador (D)/(I)**, que confunde a leitura. Formato:
   Par | Gap | Fator | Gap com a nossa escolha na tradução | Gap com a
   escolha da tradução no nosso código. As duas colunas dão o gap da
   tradução sobre o nosso depois da mudança, cada uma com a mudança num dos
   programas; "—" quando o sentido não foi medido. Nos pares com os dois
   sentidos (Go G2, Swift S3, JavaScript `Float64Array`, Java três vetores)
   a tabela mostra a confirmação nos dois.
5. JavaScript no lin318: medir o par com `--no-concurrent-recompilation` nas
   4 instâncias (rodada `runs/jsncr-20261122/`, encadeada depois das
   instâncias ruidosas) e usar esses números, com nota.

### Passos

1. **`tables.py`** (nesta pasta; o `analyze.py` continua sendo o resumo por
   rodada) — **feito em 2026-11-22** (`comparisons.csv`, 46 comparações,
   custos conferidos; o aquecimento do JIT ganhou a coluna `both`, mudança
   nos dois programas). As linhas do JavaScript trocam pelos números da
   rodada `jsncr` quando ela terminar. Plano original: uma lista de comparações (numerador,
   denominador, rodada, sentido D/I, linha da tabela); para cada uma, tempo
   médio, cv e razão por instância, gap e gap sem o berlin52; confere os
   custos; grava `comparisons.csv` e as linhas LaTeX das tabelas.
2. **Feito em 2026-11-22** (artigo `6f60713`): pasta com `tables.py`,
   `comparisons.csv`, README e `runs/` das 7 rodadas usadas; os antigos em
   `results/factor_analysis/earlier/`; mapa do `results/README.md`
   atualizado. Plano original: **Pasta do artigo** `results/factor_analysis/2026-11-22/`: os CSVs das
   rodadas usadas (só as colunas de tempo e custo), o `comparisons.csv`, o
   script e um README com rodada, commit e carga de cada comparação; os
   `gaps_*.csv` antigos vão para `results/factor_analysis/earlier/` (se a
   decisão 1 for (a)); atualizar o mapa do `results/README.md`.
3. **Feito em 2026-11-22** (`6f60713`): as quatro tabelas de gaps no
   formato de duas colunas, números novos sem marcador, linhas antigas com
   asterisco (análise anterior) e † (tradução de versão anterior); a última
   coluna da linha do índice usa o nosso C++ com índice de 64 bits
   (`ours_cplusplus.index64`, 5 comparações novas no `tables.py`). Falta: a
   conferência da tabela de desvios (7.4) e trocar as linhas do JavaScript
   pela rodada com a flag. Plano original: **Tabelas no `main.tex`** (com `git pull` antes): reconstruir as quatro
   de gaps a partir do `analyze.py`; tirar ou trocar a das mais rápidas
   (decisão 3); conferir a de desvios (7.4): passagem por valor do Swift →
   C++ continua, mas vale 1 ponto e não o gap; nota a sobre PyPy continua
   antiga.
4. **Parágrafo de método** da 7.5: protocolo novo (máquina, instâncias,
   execuções, mesma rodada, métrica, D/I) no lugar do atual (outra máquina,
   bier127).
5. **Texto da 7.5 e da 7.6, estrutura primeiro:** um esqueleto por
   parágrafo, cada afirmação com a célula da tabela que a sustenta, para
   revisão antes da prosa. Afirmações do texto atual que caem: "Swift → C++:
   a passagem por valor explica o gap", "C → C++: o índice com sinal cancela
   as matrizes aninhadas" (análise antiga), "o `unsigned` custa 4 a 6% ao
   nosso C++" (o nosso não paga nada, P4), "o resto do Julia não tem
   explicação" (P10/P12), o "Go: 3 a 4% de geração de código" (agora 1,4%).
6. **Parágrafo de fechamento** (o `\todo` da 7.6): lições práticas, com o
   índice `unsigned` (o compilador precisa provar que não é negativo), ε,
   matriz e as armadilhas de JIT (aquecimento, `Float64Array`, recompilação
   concorrente do V8).
7. **`README.md`** desta pasta: como montar as variantes, rodar no servidor e
   gerar as tabelas.
8. **TCC:** espelhar as tabelas na monografia depois que o artigo fechar.

## Fase 2: lista de fatores proposta (para revisão do autor)

Diferenças de projeto entre cada tradução e a nossa implementação atual
(`c205463`), com os candidatos dos relatórios antigos
(`doc/tcc/relatorios-ia/`). Sentido: **D** = direta (a escolha da tradução no
nosso código), **I** = inversa (a nossa escolha na tradução). Nos pares em
que a tradução ganha, a D é a principal. O gap de base de cada par está na
seção "Fase 1: resultado".

### C++ → Java (`cpp_java`)
- **J1. Acessores no laço quente.** O nosso Java lê as subsequências por
  `getSeq(i, j, k)` (com `Z` estático não final), a rota por
  `getSolut().get(i)` (`tArray`) e a matriz por `getCost(a, b)`; a tradução
  indexa direto `int[]` e `double[]`. D: acesso direto no nosso. I: os
  nossos acessores na tradução.
- **J2. Valores da rota copiados para variáveis locais** (`s_i`, `s_i_next`,
  ...) no nosso; a tradução não copia. Relatório antigo: tirar a cópia
  deixou o nosso 7,5% mais rápido; pô-la na tradução, 13 a 16% mais lento.
  D e I.
- **J3. Subsequências AoS (nosso) × três vetores SoA (tradução).** Relatório
  antigo: o AoS é melhor. I (conferência: a tradução ganharia com o AoS).
- **J4. `tArray` com cópias e deslocamento manual (nosso) × `int[]` com
  `System.arraycopy` (tradução).** Antigo: 1 a 5%. D.

### C++ → Julia (`cpp_julia`)
- **L1. Compilação JIT dentro do trecho medido.** O nosso Julia gasta cerca
  de 0,03 s a mais que a tradução entre 42 e 76 clientes (rodada de 20/11).
  Teste de medição, sem mudar código: imprimir o tempo de compilação
  acumulado (`Base.cumulative_compile_time_ns`) antes e depois da busca, nos
  dois programas.
- **L2. Campos lidos no laço quente.** O nosso lê `data.dimen` e
  `data.EPSILON` (campos de `mutable struct`) dentro dos laços; a tradução
  usa `dimen` local e `DBL_EPSILON` constante. D.
- **L3. Subsequências em `Array{Float64,3}` indexado `[T, j, i]` (nosso) ×
  vetor de `struct` imutável com `to_1D` (tradução).** Antigo: sem efeito
  mensurável. D e I.

### C++ → Swift (`cpp_swift`)
- **S1. Acessores** `s.getSeq(i, j)` e `data.getDistance(a, b)` (nosso) ×
  acesso direto `solut.seq[to1D(...)]` e `data.cost[a][b]` (tradução). D.
- **S2. Valores da rota copiados para locais** (`sIPrev`, `sI`, ...) no
  nosso. D (em Rust e Java isso piorou o nosso).
- **S3. Cópia da solução na perturbação** (`var s = s`, que copia `seq`
  inteira ao escrever). Antigo: 0,5 a 2%. D.
- **S4. `ContiguousArray` (nosso) × `Array` (tradução); laços `for` (nosso) ×
  `while` (tradução).** Antigo: `for` deixou o nosso mais rápido. Baixa
  prioridade.
- Diagnóstico sem variante: contar `cmov` nas três buscas dos dois binários
  (efeito de geração de código do LLVM, relatório antigo).

### C++ → Go (`cpp_go`)
- **G1. ε na atualização do melhor movimento:** o nosso Go não tem (é o único
  dos nossos sem ε). Antigo: pôr ε no nosso deu −3 a −4%. D.
- **G2. Acessor `solut.at(i, j)` (método com ponteiro, lê `s.dim`) × índice
  direto `solut.seq[to_1D(...)]`.** Antigo: −4 a −5% no nosso. D.
- O relatório antigo atribuiu o gap a ruído de geração de código do
  compilador Go; se G1 e G2 não confirmarem nos dois sentidos, fica assim.

### C++ → Fortran (`cpp_fortran`)
- Sem gap na base (−0,1%): fora.

### C++ → JavaScript (`cpp_javascript`)
- **JS1. Subsequências:** o nosso usa um `Array` comum criado com
  `Array(n)` (elementos com buracos), intercalado, com os offsets `T`, `C`,
  `W` lidos de `info` em tempo de execução e a variável `seq` global
  implícita; a tradução usa três `Float64Array` e constantes. D e I.
- **JS2. `to_1D` como `var` de função (nosso) × função declarada
  (tradução).** D.
- Os dois lados usam linhas `Float64Array` na matriz e ε.

### Rust → C++ (`rust_cpp`)
- **R1. Tipo do índice:** `size_t` de 64 bits (tradução) × `unsigned` de 32
  bits (nosso). Antigo: D −4,5 a −6,4%; I +48 a +54%. D e I.
- **R2. Matriz de distâncias:** vetor plano (tradução) × `double**`
  (nosso). Antigo: plano custa 4 a 8% no nosso. I.
- **R3. Rota `vector<size_t>` × `vector<int>`.** Antigo: 1 a 3%. I.

### Traduções to-C++ que eram só de base (se o autor aprovar as variantes)
Os gaps medidos (12 a 34%) não batem com a explicação antiga. Candidatos,
todos do relatório `REVERSE_CPP_PERFORMANCE_ANALYSIS.md` e do mapa do artigo;
sentido I (a nossa escolha na tradução), com D onde o fator é medido no
nosso C++:
- **Comuns às quatro:** índice com sinal ou de 64 bits (I, com o nosso
  `unsigned`; antigo: favorece a tradução em 4 a 5%) e o layout da matriz de
  distâncias (I, com o nosso `double**`; antigo: plano custa 4 a 8% e
  aninhado 5 a 7%).
- **C# → C++:** subsequências em `vector<double>` intercalado com offsets
  constantes; acessores `GetSeq`/`GetCost` (métodos de classe) no laço
  quente; `std::find` fora do laço (antigo: sem efeito).
- **Swift → C++:** passagem da solução por valor na perturbação e
  `localBest = current` (antigo: 2 a 4%); impressão do progresso em toda
  iteração (já fora do tempo, capturada em memória).
- **Fortran → C++:** melhor movimento sem ε (antigo: 13 a 17%); índices a
  partir de 1 com conversões para `size_t`; subsequências em colunas
  (column-major).
- **Go → C++:** melhor movimento sem ε (antigo: 9 a 13%); matriz em
  `vector<vector<double>>`; ponteiros para a solução e `tData` por valor.

## Decisões do autor sobre a Fase 2 (2026-11-21)

- Lista de fatores aprovada, sem o S4 (`ContiguousArray` e laços `for` no
  Swift).
- C#, Swift, Fortran e Go → C++ passam a ter variantes, só com os fatores já
  conhecidos (tipo do índice, layout da matriz e os específicos de cada uma).
- JavaScript no lin318 (cerca de 333 s nos dois lados): registrado, não
  investigado agora.
- O trabalho segue na sessão principal, sem subagente (custo de tokens).
- Antes da Fase 3: repetir o nosso Java nas quatro instâncias, com o mesmo
  protocolo (`run.sh` sem variantes, `3bf8628`), para conferir o gap de −11%
  do C++ → Java (`runs/java-20261121/`).

### Repetição do nosso Java (2026-11-21, 14:28 a 14:32 do servidor)

40 execuções, custos iguais aos da base, carga com os 4 processos da rodada e
os 7 do outro usuário. O nosso Java repete a rodada de base e a das figuras:

| Instância | Repetição (cv) | Base 20/11 | Figuras | Repetição ÷ base | Tradução ÷ repetição |
|---|---|---|---|---|---|
| berlin52 | 0,775 s (8,5%) | 0,785 s | 0,845 s | 0,987 | 0,883 |
| kroA100 | 3,079 s (1,6%) | 3,077 s | 3,071 s | 1,001 | 0,804 |
| gr202 | 22,93 s (0,4%) | 22,90 s | 22,86 s | 1,001 | 0,949 |
| lin318 | 65,09 s (0,9%) | 65,07 s | 65,87 s | 1,000 | 0,929 |

Gap com a repetição: −10,9% (base: −11,1%). O gap é da tradução, não da
medição do nosso Java; segue para a Fase 3 com J1 a J4.

## Fase 4: resultados

### Lote Java (2026-11-21, `runs/java-factors-20261121/`)

280 execuções (7 programas x 4 instâncias x 10), todos os custos certos;
carga com os 4 processos da rodada e os 7 do outro usuário (22 de 28
leituras com 11 acima de 50%). Razão sobre o nosso Java na mesma rodada
(berlin52, kroA100, gr202, lin318) e gap médio; o berlin52 tem cv de 4 a
12%, então vale também o gap sem ele:

| Variante | Razões | Gap | Gap sem berlin52 |
|---|---|---|---|
| `cpp_java.base` (tradução) | 0,907 0,820 0,951 0,926 | −9,9% | −10,1% |
| `ours_java.no_local_copies` (J2, D) | 0,934 0,970 0,953 0,944 | −5,0% | −4,4% |
| `ours_java.direct_access` (J1, D) | 1,056 0,958 0,981 0,978 | −0,7% | −2,7% |
| `ours_java.z_final` (J1a, D) | 0,996 1,030 1,030 1,026 | +2,0% | +2,9% |
| `ours_java.arraycopy_shift` (J4, D) | 1,055 1,000 1,002 1,009 | +1,6% | +0,4% |
| `cpp_java.aos` (J3, I) | 0,863 0,862 0,939 0,956 | −9,5% | −8,1% |

Leitura:
- **J2 explica cerca de metade do gap:** sem as cópias locais da rota, o
  nosso Java fica 4 a 5% mais rápido, em todas as instâncias (o relatório
  antigo dava 7,5%).
- **J1 (acessores) ajuda pouco:** 2 a 4% a partir do kroA100.
- **H-J1a refutada:** o `Z` como `static final` deixou o nosso Java 3% mais
  lento, e não mais rápido. Fica registrado como resultado contrário à
  hipótese; o mecanismo não foi investigado.
- **J4 sem efeito.**
- **J3:** o AoS na tradução muda pouco no total (−9,5% contra −9,9%), mas
  piora a tradução no kroA100 (0,862 contra 0,820) e no lin318 (0,956
  contra 0,926) e melhora no gr202.
- **Falta explicar:** J1 e J2 somados dariam uns 7%, contra 10% do gap; o
  kroA100 (0,82) continua o caso mais distante. Próximo teste: a variante
  combinada `ours_java.direct_no_locals` (J1 + J2) e o diagnóstico de
  aquecimento do JIT (H-J5).

### Lote principal (2026-11-21, 15:10 a 18:19 do servidor, `runs/factors-20261121/`)

1.560 execuções (39 programas x 4 instâncias x 10), commit `d8b35cf`, todos
os custos certos, nenhum timeout. Carga: 167 de 190 leituras com 11
processos acima de 50% (os 4 da rodada e os 7 do outro usuário), 12 com 12.
cv até 1,6%, exceto duas execuções isoladas: uma de 122,7 s entre nove de
cerca de 336 s em `ours_javascript.no_globals` no lin318, e o nosso Julia no
berlin52 com dois patamares (4 execuções em 0,134 s e 6 em 0,169 s). Os
gaps de base repetem a Fase 1 (diferença de até 0,7 ponto).

Tempos do nosso lado (média, s): C++ 0,093, 1,201, 11,59, 36,01; Go 0,226,
2,965, 28,49, 80,84; Swift 0,123, 1,665, 16,65, 49,59; Julia 0,155, 1,712,
16,27, 47,57; JavaScript 0,421, 4,337, 43,08, 336,6.

Razões (berlin52, kroA100, gr202, lin318) sobre a nossa implementação da
linguagem de destino:

**To-C++ (sentido I, a nossa escolha na tradução):**

| Variante | Razões | Gap | Sem berlin52 |
|---|---|---|---|
| `rust_cpp.base` | 1,072 1,095 1,090 1,077 | +8,4% | +8,7% |
| `rust_cpp.row_matrix` | 0,977 0,993 0,989 0,997 | −1,1% | −0,7% |
| `rust_cpp.index32` | 1,395 1,468 1,475 1,422 | +44,0% | +45,5% |
| `csharp_cpp.base` | 1,145 1,138 1,127 1,091 | +12,5% | +11,9% |
| `csharp_cpp.row_matrix` | 0,995 0,998 1,005 0,994 | −0,2% | −0,1% |
| `csharp_cpp.index32` | 1,256 1,287 1,314 1,254 | +27,8% | +28,5% |
| `swift_cpp.base` | 1,129 1,137 1,135 1,094 | +12,4% | +12,2% |
| `swift_cpp.row_matrix` | 1,011 1,014 1,025 1,009 | +1,5% | +1,6% |
| `swift_cpp.route_copy` | 1,117 1,142 1,112 1,086 | +11,4% | +11,3% |
| `swift_cpp.index32` | 1,146 1,143 1,143 1,107 | +13,5% | +13,1% |
| `go_cpp.base` | 1,255 1,285 1,310 1,301 | +28,7% | +29,8% |
| `go_cpp.epsilon` | 1,023 1,022 1,026 1,018 | +2,2% | +2,2% |
| `go_cpp.row_matrix` | 1,220 1,260 1,275 1,272 | +25,7% | +26,9% |
| `go_cpp.index32` | 1,323 1,341 1,343 1,336 | +33,6% | +34,0% |
| `fortran_cpp.base` | 1,365 1,353 1,335 1,310 | +34,1% | +33,3% |
| `fortran_cpp.epsilon` | 1,138 1,107 1,074 1,038 | +8,9% | +7,3% |
| `fortran_cpp.row_matrix` | 1,211 1,241 1,241 1,251 | +23,6% | +24,4% |
| `fortran_cpp.index32` | 1,501 1,503 1,471 1,416 | +47,3% | +46,3% |
| `ours_cplusplus.index64` (D) | 0,995 1,002 0,990 1,000 | −0,3% | −0,3% |

**From-C++ (sentido D, a escolha da tradução no nosso):**

| Variante | Razões | Gap | Sem berlin52 |
|---|---|---|---|
| `cpp_go.base` | 0,973 0,968 0,972 0,975 | −2,8% | −2,8% |
| `ours_go.direct_index` (G2) | 0,990 0,985 0,985 0,986 | −1,4% | −1,5% |
| `ours_go.epsilon` (G1) | 1,009 1,014 1,013 1,001 | +0,9% | +0,9% |
| `cpp_swift.base` | 0,999 0,992 0,963 0,985 | −1,5% | −2,0% |
| `ours_swift.route_copy` (S3) | 0,989 0,994 0,976 0,984 | −1,4% | −1,5% |
| `ours_swift.no_local_copies` (S2) | 1,001 0,997 0,983 0,988 | −0,8% | −1,1% |
| `ours_swift.direct_access` (S1) | 0,997 0,999 0,992 0,998 | −0,4% | −0,4% |
| `cpp_julia.base` | 0,830 0,966 0,959 0,960 | −7,1% | −3,8% |
| `ours_julia.struct_vector` (L3) | 0,870 1,012 0,999 0,997 | −3,1% | +0,3% |
| `ours_julia.hoist_fields` (L2) | 1,077 1,002 1,001 1,001 | +2,0% | +0,2% |
| `cpp_javascript.base` | 1,137 1,196 1,192 0,987 | +12,8% | +12,5% |
| `ours_javascript.float64array` (JS1a) | 1,163 1,194 1,164 1,044 | +14,1% | +13,4% |
| `ours_javascript.const_offsets` (JS1b) | 0,907 0,905 0,901 0,970 | −7,9% | −7,5% |
| `ours_javascript.no_globals` (JS1c) | 0,998 1,004 1,007 0,936 | −1,4% | −1,8% |
| `ours_javascript.function_decl` (JS2) | 0,996 1,003 1,002 0,999 | 0,0% | +0,1% |

Leitura:
- **Matriz de distâncias plana (to-C++):** com o nosso `double**`, o gap
  some no Rust (+8,4% → −1,1%), no C# (+12,5% → −0,2%) e no Swift
  (+12,4% → +1,5%). É o fator principal desses três pares. O relatório
  antigo dizia o contrário para o Rust (o plano custaria 4 a 8% no nosso).
- **ε (Go e Fortran → C++):** explica quase todo o gap do Go (+28,7% →
  +2,2%) e três quartos do Fortran (+34,1% → +8,9%). A matriz em linhas
  ajuda os dois (3 e 10 pontos). Próximo teste: os dois juntos
  (`go_cpp.epsilon_row_matrix`, `fortran_cpp.epsilon_row_matrix`, conferidos
  localmente).
- **Índice de 32 bits sem sinal:** truncar o índice para `unsigned` piora
  todas as traduções (Rust +36 pontos, Fortran +13, C# +15, Go +5, Swift
  +1), mas usar 64 bits no nosso C++ não muda nada (−0,3%). O R1 antigo
  repete no sentido I (+48 a +54% antes, +44% agora) e não repete no D
  (antes −4,5 a −6,4%). Ver H-R1b.
- **Swift → C++:** a cópia da solução na perturbação (`route_copy`) vale só
  1 ponto; o gap é da matriz.
- **Go (from-C++):** G2 explica metade do gap (−1,4% de −2,8%). G1
  refutado: o ε deixa o nosso Go 1% mais lento (o relatório antigo dava −3
  a −4%).
- **Swift (from-C++):** S3 explica o gap sozinho (−1,4% de −1,5%); S1 e S2
  dão menos de 1%. É o H-S3b (copy-on-write das subsequências), com efeito
  pequeno.
- **Julia:** nenhuma variante explica o gap sem o berlin52 (−3,8%). L2 sem
  efeito; L3 só muda o berlin52, onde o nosso Julia tem dois patamares (ver
  acima), então não conta. Fica o L1 (compilação JIT), a medir.
- **JavaScript:** o `Float64Array` no nosso reproduz o gap (+14,1% contra
  +12,8%, com o mesmo perfil por instância, menor no lin318). As constantes
  literais nos offsets vão no sentido oposto (−7,9%): a tradução ganha com
  elas e perde mais com o `Float64Array`. H-JS3 refutada (ver abaixo).
  JS1c e JS2 sem efeito (a execução isolada de 122,7 s puxa o lin318 do
  JS1c para 0,936). O lin318 do nosso JavaScript continua lento: 336 s,
  7,8 vezes o gr202 (C++: 3,1).

### Lote das combinações (2026-11-21, 18:32 a 18:48 do servidor, `runs/combo-20261121/`)

200 execuções (5 programas x 4 instâncias x 10), commit `ba85d9f`, todos os
custos certos. Carga: 11 de 15 leituras com 11 processos acima de 50%. O
nosso C++ e o nosso Java repetem os lotes anteriores (diferença abaixo de
1%, exceto o berlin52 do Java, 0,803 s contra 0,775 s).

| Variante | Razões | Gap | Sem berlin52 |
|---|---|---|---|
| `go_cpp.epsilon_row_matrix` (I) | 1,017 1,011 1,008 1,008 | +1,1% | +0,9% |
| `fortran_cpp.epsilon_row_matrix` (I) | 0,981 0,979 0,967 0,969 | −2,6% | −2,8% |
| `ours_java.direct_no_locals` (J1 + J2, D) | 0,920 0,919 0,950 0,945 | −6,6% | −6,2% |

Leitura:
- **Go → C++:** ε e matriz em linhas juntos levam o gap de +28,7% a +1,1%.
- **Fortran → C++:** os dois juntos levam o gap de +34,1% a −2,6%: a
  tradução passa a ganhar do nosso por pouco.
- **Java:** J1 + J2 juntos dão −6,2% sem o berlin52, perto da soma dos
  dois separados (−7,1%). No gr202 e no lin318, a variante chega à razão
  da tradução (0,950 e 0,945 contra 0,951 e 0,926 no lote Java): o gap
  dessas duas instâncias fica explicado. No kroA100 não fica (0,919 contra
  0,820). O que falta explicar está no kroA100 e no berlin52, as duas
  instâncias curtas, o que reforça a H-J5 (aquecimento do JIT).

## Pendências de investigação (lista de trabalho, 2026-11-21)

Pedido do autor: tratar os itens em ordem; problema ou hipótese nova que
possa esperar vai para o fim da lista.

- [x] **P1. Java nas instâncias curtas (H-J5).** J1 + J2 explicam o gr202 e
  o lin318, mas não o kroA100 (0,919 contra 0,820) nem o berlin52.
- [x] **P2. Julia sem o berlin52 (L1).** Nenhuma variante explica os −3,8%;
  medir a compilação dentro do trecho medido. Também: os dois patamares do
  nosso Julia no berlin52 (0,134 s e 0,169 s).
- [x] **P3. Metade do gap do C++ → Go.** G2 explica −1,4% de −2,8%; G1
  refutado.
- [x] **P4. Índice `unsigned` (H-R1b).** Custa de 1% a 36% nas traduções
  para C++ e nada no nosso C++.
- [x] **P5. `Float64Array` (H-JS4).** Confirmar na inversa: as três
  subsequências da tradução como `Array`.
- [x] **P6. Inversas dos pares em que a tradução ganha.** O escopo pede a
  inversa para confirmar a direta; só o Java (J3) tem. Go (G2: o acessor
  `at()` na tradução) e Swift (S3: a cópia da solução na tradução); o
  JavaScript é o P5.
- [x] **P7. lin318 do JavaScript (H-JS5).** 7,8 vezes o gr202 nos dois lados
  (C++: 3,1). Adiado pelo autor em 2026-11-21; entra por último.
- [x] **P8. Mecanismo do `index32` em C# → C++ e Fortran → C++.** Nesses
  dois, o patch trunca uma conta feita em `int` (a mesma forma do nosso
  `to_1D`) e ainda custa 13 a 15 pontos no servidor; localmente, o C# → C++
  só perde 3% em ciclos (instruções +11%). Ver a explicação do P4.
- [x] **P9. Mesmo caminho de busca nos outros pares.** O vetor `rnd` tem
  exatamente o número de sorteios do caminho de referência; no Go, o nosso
  e a tradução consomem todos (P3). Conferir os demais da mesma forma (em
  JavaScript, ler além do fim não dá erro, ver P7).
- [x] **P10. Julia: gap de −3,3% em regime.** Depois do aquecimento, a
  tradução continua 3,1 a 3,6% mais rápida do kroA100 ao lin318 (P2), e
  nenhuma variante (L2, L3) explica. Reler os laços quentes dos dois.
- [x] **P11. Java no kroA100 depois do aquecimento.** Com aquecimento e
  J1 + J2, sobram 7% no kroA100 (0,927); no gr202 e no lin318, nada
  (1,012 e 0,991). No berlin52, uma busca de aquecimento não basta (o nosso
  continua em 0,31 s, a tradução cai para 0,19 s).
- [x] **P12. Índice de cidade a partir de 1 no Julia (H-L4).** Ver o P10.
- [x] **P13. Node 24 nas rodadas das figuras.** O lin318 lento é do Node 24
  (P7). Conferir com que Node foi medido o JavaScript das figuras e se a
  frase do artigo "JavaScript is the only language whose gap grows with
  instance size" é efeito da versão (instâncias acima de cerca de 300
  clientes).
- [x] **P14. `index32` no C# → C++: subsequências ou matriz?** O controle
  com o deslocamento do campo fora do truncamento não mudou nada no
  servidor (P8). Os dois patches truncam também o índice da matriz de
  distâncias; separar os dois.
- [x] **P15. JavaScript das figuras não se repete.** No P13, a pr299 leva
  125 s (Node 24) e 148 s (Node 22) no servidor, contra 348 s em todas as
  10 execuções das figuras; a pr439 das figuras tem dois patamares (4
  execuções em cerca de 370 s, 5 em cerca de 1.100 s). Decidir com o autor
  a versão do Node da campanha nova e medir o JavaScript de novo em todas as
  instâncias antes de manter a frase sobre o gap crescer com o tamanho.
  **Decisão do autor (2026-11-21): Node 24** (o 24.18 do `flake.nix`).
  Rodada `runs/js66-20261121/`: o nosso JavaScript nas 66 instâncias
  (`instances66.txt`, a lista `inst-list-66` do servidor), 10 execuções.
  Resultado: se repete (ver a rodada do P15).
- [x] **P16. Dois patamares do Node 24.** Na rodada do P15, a pr299 tem uma
  execução em 125 s e nove em cerca de 348 s, e a pr439 uma em 350 s e nove
  em cerca de 970 s, com o mesmo caminho de busca (o custo é o mesmo em
  todas). O patamar lento sai por sorteio de execução, não da instância.
  Hipótese H-JS6: a compilação em segundo plano do V8 (Maglev/TurboFan
  concorrentes) termina em momentos diferentes e às vezes deixa a busca no
  código pior; com 4 processos em 4 núcleos, as threads do compilador
  disputam CPU (no diagnóstico do P13, 3 processos em 4 núcleos, as três
  execuções caíram no patamar rápido). Teste: a pr299 com 10 execuções,
  `-j 4`, com `--no-concurrent-recompilation` (ou `--single-threaded`) e sem.
  Para o texto: decidir com o autor como reportar o JavaScript (média com os
  dois patamares, ou mediana) e se o lin318 e o linhp318 também têm um
  patamar rápido que não apareceu em 10 execuções.
  Resultado: H-JS6 confirmada (ver a rodada do P16).
- [x] **P17. C da d493 nas figuras.** Nos dados das figuras, o C leva 403 s
  na d493 (dispersão máx/mín de 1,31), contra 176 s na pcb442; todas as
  linguagens ficam com a razão sobre o C menor nela (Rust 1,02 contra 1,35
  nas outras quatro instâncias acima de 400, Fortran 0,80, Java 1,25,
  JavaScript 1,66). O C da d493 (e talvez o da pcb442, dispersão 1,16) está
  inflado e puxa para baixo as razões do grupo acima de 100 clientes. Medir
  de novo o C (e o C++) nas cinco instâncias acima de 400 clientes, junto
  com a campanha completa nova.
  Resultado: confirmado, e o efeito vai além do C (ver a rodada do P17 e o
  P19).
- [x] **P18. Compilação concorrente nas outras instâncias.** O P16 mostrou
  o patamar lento na pr299 vindo da recompilação concorrente do V8. Conferir
  se o lin318, o linhp318 e a pr439 também caem no patamar rápido com
  `--no-concurrent-recompilation`, e se a flag muda o tempo numa instância
  sem pico (gr202). Rodada `runs/p18-20261122/`, 5 execuções por conjunto.
  Para o texto, decidir com o autor: manter o Node 24 padrão (e explicar os
  picos) ou medir o JavaScript com a flag.
  Resultado: os quatro picos são o mesmo efeito (ver a rodada do P18).
- [ ] **P19. Instâncias ruidosas nas figuras pesam no gap acima de 100.**
  A métrica da figura (razão das médias do grupo) é dominada pelas
  instâncias maiores, e a d493 sozinha tem cerca de 20% do tempo do C acima
  de 100 clientes. Nos dados das figuras, a d493 e a pcb442 têm dispersão
  alta em várias linguagens (d493: C 1,31, C++ 1,28, Fortran 1,55, Java
  1,20; pcb442: Fortran 1,31, Java 1,25, C++ 1,19, C 1,16); Julia, Lua,
  Python, Rust e Swift ficam limpas (até 1,04 na d493). Tirar só a d493 muda
  o gap acima de 100 de quase todas: Rust 27,3 → 36,8%, Fortran 2,9 → 11,6%,
  Swift 44,9 → 59,5%, Java 74,8 → 93,6%, Julia 11,2 → 28,9%, JavaScript
  303,8 → 365,7%. Precisa medir de novo todas as linguagens nas instâncias
  ruidosas (pelo menos pcb442 e d493; com dispersão acima de 1,10 em cinco
  ou mais linguagens: pcb442, pr439, d493, gr120) ou fazer a campanha
  completa nova (pendência antiga). Estimativa só das AOT e JIT nessas
  quatro: cerca de 7 h de servidor. Decidir com o autor, junto com a flag do
  JavaScript (P18).
  **Decisão do autor (2026-11-22): só as instâncias ruidosas, e o
  JavaScript nas duas configurações.** Rodadas encadeadas: AOT e JIT (menos
  JavaScript) na d493, gr120, pcb442 e pr439 (`instances_noisy4.txt`,
  `runs/noisy4-20261122/`), depois a variante
  `ours_javascript.no_concurrent_recompilation` nas 66 instâncias
  (`runs/js66ncr-20261122/`); o JavaScript padrão é o do P15.

### Andamento das pendências

**P1 e P2 (2026-11-21).** Variantes de diagnóstico `warm`: uma busca inteira
fora do tempo antes da busca medida, no mesmo processo (o nosso Java recria
o `GILS_RVND`; nos outros, `rnd_index` volta ao início). Se o gap das
instâncias curtas sumir com o aquecimento, a causa é o JIT. Rodada
`runs/warm-20261121/`, com `cpp_java.base` e `cpp_julia.base` de novo para
comparar na mesma rodada.

Resultado (`runs/warm-20261121/`, 18:59 a 19:46 do servidor, 360
execuções, commit `3e7b545`, todos os custos certos; 38 de 48 leituras com
11 processos acima de 50%). Razões por instância (berlin52, kroA100, gr202,
lin318):

| Comparação | Razões | Gap sem berlin52 |
|---|---|---|
| Java: tradução / nosso, sem aquecimento | 0,891 0,807 0,951 0,918 | −10,8% |
| Java: tradução / nosso, com aquecimento | 0,609 0,872 0,971 0,930 | −7,6% |
| Java: nosso com aquecimento / sem | 0,407 0,761 0,967 0,989 | −9,4% |
| Java: tradução com aquecimento / sem | 0,278 0,821 0,987 1,001 | −6,3% |
| Java: tradução / nosso J1 + J2, ambos com aquecimento | 0,632 0,927 1,012 0,991 | −2,3% |
| Julia: tradução / nosso, sem aquecimento | 0,797 0,965 0,950 0,966 | −4,0% |
| Julia: tradução / nosso, com aquecimento | 0,968 0,964 0,966 0,969 | −3,3% |
| Julia: nosso com aquecimento / sem | 0,822 1,003 0,982 1,002 | −0,4% |
| Julia: tradução com aquecimento / sem | 0,997 1,001 0,999 1,005 | +0,2% |

Leitura:
- **Java (H-J5 confirmada em parte):** o aquecimento do JIT pesa nos dois
  programas até o kroA100 (o nosso perde 24% no kroA100 por ele, a
  tradução 18%) e quase nada do gr202 em diante. Com os dois aquecidos, o
  gap do kroA100 cai de 19% para 13%; com aquecimento e J1 + J2, para 7%.
  No gr202 e no lin318, J1 + J2 explicam tudo. O resto do kroA100 e o
  berlin52 vão para o P11.
- **Julia (L1 confirmada só no berlin52):** a compilação dentro do trecho
  medido custa ao nosso Julia 0,028 s no berlin52 (0,161 → 0,133 s) e nada
  na tradução, que já compila antes; isso explica toda a diferença do
  berlin52 (0,797 → 0,968) e os dois patamares (com aquecimento, cv 0,2%).
  Nas outras instâncias, o aquecimento não muda nada: os −3,3% são de
  regime e vão para o P10.

**P10 (2026-11-21, local, `perf` com jitdump do Julia, kroA100).** Por
função, a tradução gasta menos ciclos em `search_reinsertion` (−14%),
`search_swap` (−10%) e `search_two_opt` (−6%) e mais em
`update_subseq_info_matrix` (+7%); no total da busca, −7,7% (no servidor,
−3,3%). No `code_native` do `search_reinsertion`, o laço interno tem 29 a
40 instruções no nosso e 26 a 33 na tradução. As instruções a mais do
nosso: mais ponteiros de indução (o índice 3D `[campo, j, i]`) e o ajuste
do índice de cidade, que começa em 1 (`cost[s[k], s[i]]`, refeito a cada
iteração com `movabs` + `add`, por falta de registrador); a tradução guarda
as cidades a partir de 0 e soma 1 dentro de `dist`, o que vira deslocamento
fixo. Como o L3 (subsequências como vetor de `struct`) não mudou nada no
servidor, o índice 3D sozinho não explica. Hipótese H-L4 (o índice de
cidade a partir de 1 na matriz de distâncias) no fim da lista (P12).

**P11 (2026-11-21).** A inversa do J3 (AoS na tradução) piorou a tradução
no kroA100 (0,862 contra 0,820) e no lin318; a direta nunca foi testada.
Variantes `ours_java.soa` (J1 + J2 com as subsequências em três vetores,
como na tradução) e `ours_java.soa_warm` (com aquecimento), conferidas
localmente.

**P12 (2026-11-21).** Variante `ours_julia.zero_based`: a construção
devolve a rota com as cidades a partir de 0 e todo `data.cost[s[a], s[b]]`
da busca vira `data.cost[s[a]+1, s[b]+1]`, como na tradução (17 acessos);
`ours_julia.zero_based_warm` com aquecimento. Custos conferidos localmente.

**P9 (2026-11-21, local, kroA100 salvo onde indicado).** Mesmo caminho de
busca (consumo igual ao tamanho do vetor `rnd`) conferido no Go (berlin52,
kroA100, gr202; P3), no JavaScript (gr202 e lin318; P7), no Java (kroA100:
50.313 nos dois) e no Julia (kroA100: 50.313 nos dois). O resto do Java no
kroA100 (P11) e o gap do Julia em regime (P10) não vêm de caminho
diferente. Traduções para C++ (kroA100): Rust, C#, Fortran, Go e Swift →
C++ consomem os 50.313 sorteios. Só o Swift nativo (o nosso e o C++ →
Swift) ficou sem conferir, por falta do compilador na máquina local; o gap
dele já está explicado pelo S3.

### Rodada das pendências P3 a P8 (`runs/p3p8-20261121/`, 19:50 a 21:21 do servidor)

640 execuções (16 programas x 4 instâncias x 10), commit `0f24a3e`, todos os
custos certos; 81 de 92 leituras com 11 processos acima de 50%. Razões
(berlin52, kroA100, gr202, lin318):

| Comparação | Razões | Gap | Sem berlin52 |
|---|---|---|---|
| `cpp_go.base` / nosso Go | 0,971 0,967 0,968 0,979 | −2,9% | −2,9% |
| `ours_go.style` (G2 + campos + RVND) / nosso Go | 0,995 0,985 0,975 0,987 | −1,5% | −1,8% |
| `cpp_go.at_accessor` (G2 inversa) / `cpp_go.base` | 1,032 1,019 1,013 1,013 | +1,9% | +1,5% |
| `cpp_swift.full_copy` (S3 inversa) / `cpp_swift.base` | 1,012 1,006 1,018 1,012 | +1,2% | +1,2% |
| `cpp_javascript.base` / nosso JavaScript | 1,138 1,194 1,196 0,954 | +12,1% | +11,5% |
| `cpp_javascript.array` / nosso JavaScript | 1,012 1,027 1,063 0,933 | +0,9% | +0,8% |
| `cpp_javascript.array` / `cpp_javascript.base` | 0,890 0,860 0,888 0,979 | −9,6% | −9,1% |
| `ours_cplusplus.size_t_counter` / nosso C++ | 1,263 1,349 1,362 1,318 | +32,3% | +34,3% |
| `rust_cpp.index32` / nosso C++ | 1,369 1,465 1,476 1,425 | +43,4% | +45,6% |
| `rust_cpp.index32_int` / nosso C++ | 1,142 1,164 1,148 1,134 | +14,7% | +14,8% |
| `csharp_cpp.index32` / nosso C++ | 1,231 1,284 1,294 1,256 | +26,6% | +27,8% |
| `csharp_cpp.index32_field` / `csharp_cpp.index32` | 0,993 1,005 1,006 0,994 | −0,1% | +0,2% |

Leitura:
- **P3 (Go):** a ordem dos campos e a lista do RVND não somam nada ao G2
  (−1,5% contra −1,4% do G2 sozinho). A outra metade do gap (cerca de
  1,4%) fica como efeito de geração de código do compilador Go, como
  concluía o relatório antigo; o caminho de busca é o mesmo (P9).
- **P6 (inversas):** o acessor `at()` na tradução a deixa 1,5 a 1,9% mais
  lenta (G2 confirmado nos dois sentidos); a cópia da solução inteira no
  C++ → Swift, 1,2% (S3 confirmado nos dois sentidos).
- **P5 (JavaScript, H-JS4 confirmada):** com as três subsequências em
  `Array`, a tradução fica 9,6% mais rápida e empata com o nosso (+0,9%).
  O `Float64Array` explica o gap do C++ → JavaScript nos dois sentidos.
- **P4 (índice):** o contador `size_t` truncado pelo `to_1D(int, ...)`
  deixa o nosso C++ 32% mais lento no servidor (24% de ciclos na máquina
  local): o mecanismo está confirmado. No Rust → C++, fazer a conta em
  `int` tira dois terços do custo no servidor (+43% → +15%), embora quase
  não mude o número de instruções na máquina local; o resto vem dos
  contadores `size_t` truncados (`(int)col`). Conclusão para o texto: o
  custo é de truncar para 32 bits um índice calculado a partir de
  contadores de 64 bits; com contadores `int`, como no nosso C++, o
  `unsigned` não custa nada.
- **P8 (C#):** tirar o deslocamento do campo do truncamento não muda nada
  no servidor (−0,1%), apesar de tirar um terço das instruções a mais na
  máquina local. Hipótese refutada; a causa no C# → C++ fica aberta (P14).

**P7 (2026-11-21, servidor, `runs/node-20261121/`).** O nosso JavaScript,
uma execução por versão, gr202 e lin318 em paralelo:

| Node | gr202 | lin318 | lin318 / gr202 |
|---|---|---|---|
| 20.19.6 (nixos-25.05) | 59,7 s | 168,9 s | 2,8 |
| 22.20.0 (nixos-25.05) | 49,7 s | 143,6 s | 2,9 |
| 24.18.0 (flake das rodadas) | 42,5 s | 333,1 s | 7,8 |

O lin318 lento é uma regressão do Node 24: ele é o mais rápido no gr202 e
2,3 vezes mais lento que o Node 22 no lin318. Não é da CPU nem do código
(H-JS5 resolvida). Consequência para o artigo: P13.

**P13 (2026-11-21).** Nos dados das figuras (`results/main/cleaned_data.csv`
do artigo, rodadas de 16 a 18/11), o JavaScript fica entre 3,3 e 4,0 vezes o
C em todas as instâncias acima de 150 clientes, menos quatro: pr299 (9,48),
linhp318 (9,58), lin318 (9,50) e pr439 (6,49); a rd400, maior, fica em 3,52.
É o mesmo padrão do Node 24 no lin318 (P7). Sem essas quatro, o JavaScript
fica 228,8% mais lento que o C acima de 100 clientes (com elas na mediana
das demais, 3,67 vezes o C: 234,8%), contra 303,8% na figura e 282,8% até
100 clientes. Ou seja, a frase do artigo de que o gap do JavaScript cresce
com o tamanho da instância vem dessas quatro instâncias. O servidor não tem
Node fora do container, e a versão usada nas rodadas de 16 a 18/11 não está
registrada (pendência antiga do TODO do artigo). Próximo passo:
`diag/node_versions.sh` com Node 22 e 24 na pr299, na pr439 e na rd400
(controle), depois da rodada do P11 e do P12.

**P14 (2026-11-21, local, `perf`, kroA100, média de 3).** Separando os dois
truncamentos do `index32` do C# → C++: só nas subsequências, 14,95 bilhões
de instruções e 5,52 bilhões de ciclos; só na matriz, 12,42 e 4,87; base,
12,94 e 5,05. O custo é todo das subsequências. No laço de
`update_subseq_info_matrix`, a variante relê da memória o T e o C do
elemento anterior e calcula um endereço para cada campo (`lea` + 3, + 4,
+ 5 em 32 bits), o mesmo sintoma do Rust → C++ (P4). Mecanismo: o
`unsigned` só sai de graça quando o compilador prova que o índice não é
negativo. No nosso C++, o passo do índice (`dimen`) é o mesmo valor do
limite do laço, e daí sai a prova; no C# → C++, o passo é o `d` de outro
objeto (`tSolution`), que pode ser negativo para o compilador. Teste:
declarar `i`, `j`, `k` e `d` não negativos (`__builtin_unreachable`) volta a
12,97 bilhões de instruções e 5,29 bilhões de ciclos. Variantes
`csharp_cpp.index32_seq`, `csharp_cpp.index32_matrix` e
`csharp_cpp.index32_seq_assume` para a próxima rodada. Conclusão para o
texto (junto com o P4): truncar o índice para 32 bits só não custa quando o
compilador prova que ele é não negativo e calculado de contadores de 32
bits; o nosso C++ satisfaz as duas condições por acaso, pelo jeito como os
laços estão escritos.

### Rodada do P11 e do P12 (`runs/p11p12-20261121/`, 21:36 a 22:46 do servidor)

560 execuções (14 programas x 4 instâncias x 10), commit `ef8bbd8`, todos os
custos certos; 63 de 71 leituras com 11 processos acima de 50%. Razões
(berlin52, kroA100, gr202, lin318):

| Comparação | Razões | Gap | Sem berlin52 |
|---|---|---|---|
| `cpp_java.base` / nosso Java | 0,901 0,801 0,949 0,925 | −10,6% | −10,8% |
| `ours_java.direct_no_locals` (J1 + J2) / nosso Java | 0,984 0,917 0,946 0,939 | −5,3% | −6,6% |
| `ours_java.soa` (J1 + J2 + SoA) / nosso Java | 0,947 0,865 0,951 0,925 | −7,8% | −8,6% |
| tradução / J1 + J2, ambos com aquecimento | 0,615 0,919 1,011 0,988 | −11,7% | −2,7% |
| tradução / J1 + J2 + SoA, ambos com aquecimento | 0,622 0,992 1,013 1,001 | −9,3% | +0,2% |
| `cpp_julia.base` / nosso Julia | 0,810 0,966 0,963 0,965 | −7,4% | −3,5% |
| `ours_julia.zero_based` / nosso Julia | 0,817 0,963 0,966 0,973 | −7,0% | −3,3% |
| tradução / nosso, ambos com aquecimento | 0,965 0,964 0,961 0,963 | −3,7% | −3,8% |
| tradução / `zero_based`, ambos com aquecimento | 0,982 1,005 0,995 0,989 | −0,7% | −0,4% |

Leitura:
- **P11 (Java, explicado):** as subsequências em três vetores (a direta do
  J3) deixam o nosso Java mais rápido no kroA100 (0,865 contra 0,917 só com
  J1 + J2) e não mudam o gr202 e o lin318. Com aquecimento, J1 + J2 + SoA
  chegam à tradução do kroA100 ao lin318 (0,992, 1,013, 1,001). O gap do
  C++ → Java fica assim: J2 (cópias locais da rota), J1 (acessores), SoA no
  kroA100 e aquecimento do JIT nas instâncias curtas. O berlin52 continua
  0,62 com uma busca de aquecimento só: o JIT ainda não terminou de
  compilar o nosso programa nessa instância.
- **P12 (Julia, H-L4 confirmada):** guardar as cidades a partir de 0 e somar
  1 no acesso à matriz deixa o nosso Julia 2,7 a 3,7% mais rápido do
  kroA100 ao lin318, o mesmo tamanho do gap; com aquecimento, a tradução e
  a variante empatam (−0,4%). O gap do C++ → Julia fica explicado: a
  compilação JIT dentro do tempo medido no berlin52 (P2) e o índice de
  cidade a partir de 1 nas outras instâncias (o P10 tinha achado o
  mecanismo no assembly).

**P13 (2026-11-21, servidor, `runs/node2-20261121/`, uma execução por
versão, três instâncias em paralelo).**

| Instância | Node 22.20 | Node 24.18 | Figuras (média de 10) |
|---|---|---|---|
| pr299 | 147,7 s | 125,1 s | 348,2 s |
| pr439 | 410,3 s | 348,9 s | 742,2 s |
| rd400 | 562,9 s | 491,2 s | 535,4 s |

O Node 24 não fica lento na pr299 nem na pr439, ao contrário do lin318
(P7). Nas figuras, a pr299 leva 348 s nas 10 execuções, 2,8 vezes o tempo
de agora, e a pr439 tem dois patamares (4 execuções em cerca de 370 s, que
é o tempo de agora, e 5 em cerca de 1.100 s). Então só o lin318 (e
provavelmente o linhp318) se explica pelo Node 24; a pr299 e a pr439 das
figuras não se repetem com nenhuma das duas versões, e a versão usada
naquelas rodadas não está registrada. A conclusão do P13 continua: a frase
de que o gap do JavaScript cresce com o tamanho vem de instâncias com
tempos que não se repetem. Fica o P15.

### Rodada do P14 (`runs/p14-20261121/`, 23:21 a 23:35 do servidor)

240 execuções (6 programas x 4 instâncias x 10), commit `fe245f4`, todos os
custos certos; 14 de 15 leituras com 11 processos acima de 50%. Razões sobre
`csharp_cpp.base` (berlin52, kroA100, gr202, lin318):

| Variante | Razões | Gap |
|---|---|---|
| `csharp_cpp.index32` (subsequências e matriz) | 1,104 1,129 1,148 1,150 | +13,3% |
| `csharp_cpp.index32_seq` (só subsequências) | 1,143 1,163 1,174 1,163 | +16,1% |
| `csharp_cpp.index32_matrix` (só matriz) | 0,976 0,980 0,984 0,983 | −1,9% |
| `csharp_cpp.index32_seq_assume` (só subsequências, índices não negativos) | 1,007 1,009 1,012 1,014 | +1,1% |

Leitura: o servidor confirma o que a máquina local mostrou. O custo do
`index32` no C# → C++ é todo do índice das subsequências; truncar o índice
da matriz deixa a tradução 2% mais rápida. Declarar os índices e o `d` não
negativos tira quase todo o custo (+16,1% → +1,1%). Mecanismo confirmado:
o `unsigned` custa quando o compilador não prova que o índice é não
negativo (C# → C++) ou quando ele vem de contadores de 64 bits truncados
(Rust → C++, P4).

### Rodada do P15 (`runs/js66-20261121/`, 23:36 de 21/11 a 03:27 de 22/11 do servidor)

O nosso JavaScript nas 66 instâncias, Node 24.18 (`flake.nix`), 10
execuções, `-j 4`; 660 execuções, todos os custos iguais aos do C; carga
de cerca de 11 durante toda a rodada (como na do P14). Comparação com o
JavaScript das figuras (`results/main/cleaned_data.csv`), com o C das
figuras como base nas duas colunas:

| Instância | C (figuras) | JavaScript, figuras | JavaScript, Node 24 | Razão, figuras | Razão, Node 24 |
|---|---|---|---|---|---|
| a280 | 44,6 s | 146,8 s | 145,9 s | 3,29 | 3,27 |
| pr299 | 36,7 s | 348,2 s | 326,1 s | 9,48 | 8,88 |
| lin318 | 35,4 s | 336,1 s | 337,9 s | 9,50 | 9,55 |
| linhp318 | 46,6 s | 446,7 s | 452,2 s | 9,58 | 9,70 |
| rd400 | 152,2 s | 535,4 s | 494,3 s | 3,52 | 3,25 |
| fl417 | 135,6 s | 525,8 s | 444,5 s | 3,88 | 3,28 |
| gr431 | 120,7 s | 483,9 s | 403,7 s | 4,01 | 3,34 |
| pr439 | 114,3 s | 742,2 s | 909,6 s | 6,49 | 7,96 |
| pcb442 | 176,1 s | 665,3 s | 550,3 s | 3,78 | 3,12 |
| d493 | 403,5 s | 964,3 s | 671,0 s | 2,39 | 1,66 |

Até a a280, as duas medições do JavaScript batem (diferença abaixo de 2%
em quase todas as instâncias). Os quatro picos das figuras se repetem com
o Node 24: lin318 e linhp318 nas 10 execuções, pr299 em nove (uma em
125 s) e pr439 em nove (uma em 350 s, e cerca de 970 s nas outras; nas
figuras foram 4 em cerca de 370 s e 5 em cerca de 1.100 s). As execuções
isoladas do P13 (125 s na pr299, 349 s na pr439) caíram no patamar
rápido; daí a impressão de que as figuras não se repetiam. Conclusão: o
JavaScript das figuras foi, muito provavelmente, medido com o Node 24, e os
picos são do Node 24, não de ruído daquelas rodadas. Novo: os dois
patamares (P16). Acima de 400 clientes, o JavaScript de agora fica 8 a 30%
mais rápido que o das figuras, que tinham dispersão de 1,07 a 1,13 nessas
instâncias (agora 1,01 a 1,04).

Gap sobre o C (média dos tempos do grupo, a métrica da figura):

| Grupo | Figuras | Node 24 |
|---|---|---|
| até 100 clientes (28) | 282,8% | 283,0% |
| acima de 100 (38) | 303,8% | 273,0% |

Por faixa, a mediana da razão sobre o C com o Node 24 é 4,53 até 100
clientes, 3,70 de 101 a 300, 9,55 de 301 a 400 (três instâncias, os três
picos) e 3,28 de 401 a 500. Consequência para o texto: com a medição nova,
o gap do JavaScript não cresce com o tamanho da instância (283% contra
273%); a frase do artigo vem de dois efeitos que não são do tamanho, os
picos do Node 24 em quatro instâncias e o JavaScript mais lento das
figuras acima de 400 clientes. O C de base ainda é o das figuras, com a
d493 inflada (P17).

### Rodada do P16 (`runs/p16-20261122/`, 03:35 a 04:25 do servidor)

`diag/node_flags.sh`: o nosso JavaScript na pr299, Node 24.18, 10
execuções de cada conjunto de flags, intercaladas, 4 em paralelo; custo
certo nas 30.

| Flags | Execuções em cerca de 125 s | Execuções em cerca de 348 s | Média |
|---|---|---|---|
| padrão | 2 | 8 | 302,1 s |
| `--single-threaded` | 10 | 0 | 125,7 s |
| `--no-concurrent-recompilation` | 10 | 0 | 125,9 s |

H-JS6 confirmada: o patamar lento vem da recompilação concorrente do V8
(otimização do TurboFan em segundo plano). Desligando só ela, as 10
execuções caem no patamar rápido; `--single-threaded`, que desliga também
GC e marcação concorrentes, não muda mais nada. Com a otimização em
segundo plano, o código otimizado entra num momento que varia de execução
para execução, e em 8 de 10 a busca fica num código pior até o fim.
Próximo: P18 (as outras instâncias com pico).

### Rodada do P17 (`runs/p17-20261122/`, 04:25 a 05:20 do servidor)

C e C++ nas cinco instâncias acima de 400 clientes, 10 execuções, mesmo
`run.sh` das outras rodadas; custos certos nas 100; carga mediana de 11.
Média (dispersão máx/mín):

| Instância | C, figuras | C, agora | C++, figuras | C++, agora | JavaScript (Node 24) / C de agora |
|---|---|---|---|---|---|
| rd400 | 152,2 s (1,07) | 148,7 s (1,05) | 152,3 s (1,07) | 150,5 s (1,03) | 3,32 |
| fl417 | 135,6 s (1,05) | 133,7 s (1,03) | 137,6 s (1,05) | 135,0 s (1,04) | 3,33 |
| gr431 | 120,7 s (1,02) | 121,3 s (1,02) | 121,2 s (1,02) | 122,7 s (1,03) | 3,33 |
| pcb442 | 176,1 s (1,16) | 167,8 s (1,03) | 186,4 s (1,19) | 168,8 s (1,03) | 3,28 |
| d493 | 403,5 s (1,31) | 204,2 s (1,04) | 423,5 s (1,28) | 205,6 s (1,03) | 3,29 |

O C e o C++ da d493 nas figuras levam o dobro do tempo de agora; a pcb442
fica 5 a 10% acima. As outras três batem com as figuras (até 2%). Com o C
de agora, o JavaScript fica em 3,3 vezes o C nas cinco, o mesmo nível de
101 a 300 clientes. Com o JavaScript do Node 24 e o C de agora nessas
cinco, o gap do JavaScript acima de 100 clientes vai a 335,6% (283,0% até
100); sem os quatro picos, 237,7%. O C da d493 inflado puxava para baixo o
gap acima de 100 de todas as linguagens medidas limpas nela (P19).

### Rodada do P18 (`runs/p18-20261122/`, 05:20 a 06:25 do servidor)

`diag/node_flags.sh`, Node 24.18, 5 execuções por conjunto, intercaladas, 4
em paralelo; custo certo em todas. Tempos (s):

| Instância | Padrão | `--no-concurrent-recompilation` | Com a flag / C |
|---|---|---|---|
| gr202 | 43 43 43 43 42 | 42 42 43 42 42 | 3,73 |
| lin318 | 336 334 334 333 333 | 122 122 122 123 124 | 3,46 |
| linhp318 | 440 445 442 443 443 | 162 162 162 162 162 | 3,48 |
| pr439 | 960 974 966 970 973 | 355 352 356 354 353 | 3,10 |

(pr299 do P16: 125,9 s com a flag, 3,43 vezes o C.) Os quatro picos são o
mesmo efeito: com a recompilação concorrente desligada, o JavaScript fica
entre 3,1 e 3,5 vezes o C nas quatro, como nas instâncias vizinhas, e a
flag não muda a instância sem pico (gr202, −1,5%). No lin318 e no linhp318
o patamar lento saiu nas 15 execuções com a configuração padrão (10 do P15
e 5 daqui). Com os tempos da flag nos quatro picos e o C de agora acima de
400 clientes, o gap do JavaScript acima de 100 clientes fica em 235,9%,
contra 283,0% até 100: o gap cai com o tamanho, como nas outras linguagens
com JIT. Decisão pendente com o autor (P18 e P19).

**P3 (2026-11-21).** O resto do laço quente difere em detalhes: ordem dos
campos da subsequência (T, W, C no nosso; T, C, W na tradução), lista do
RVND (o nosso copia de uma lista-base com `make` + `copy`), `swap` e
`reverse` como funções sobre o vetor e o índice `uint` do `rnd`. Os dois
fazem o mesmo caminho de busca: instrumentados localmente, consomem
exatamente os 17.313, 50.313 e 108.601 sorteios do berlin52, do kroA100 e
do gr202 (o tamanho do vetor `rnd`). Variante `ours_go.style` (G2 + ordem
dos campos + lista do RVND) para a próxima rodada; se não chegar à razão da
tradução, o resto fica como efeito de geração de código.

**P4 (2026-11-21, local, `perf`).** No Rust → C++ (kroA100), o `index32`
sobe as instruções de 12,1 para 15,1 bilhões (+25%) e os ciclos em 33%; o
tempo de `update_subseq_info_matrix` passa de 0,81 para 1,86 bilhão de
ciclos. No assembly desse laço, a versão de base mantém o T e o C do
elemento anterior em registrador e avança o ponteiro 24 bytes por iteração;
com o `index32`, ela relê os dois da memória a cada iteração (a escrita de
uma iteração vira leitura da seguinte) e recalcula o endereço com três
`lea`. A causa não é a conta em 64 bits: a variante de controle
`rust_cpp.index32_int`, que faz a conta em `int` como o nosso `to_1D`,
fica igual (14,8 bilhões de instruções). A causa é o contador do laço:
truncar um contador de 64 bits (`size_t j`) para 32 bits tira do índice a
forma linear no contador, e o compilador deixa de provar que `(i, j-1)` e
`(i, j)` são vizinhos. Confirmação no nosso C++: trocar só o `j` de
`update_subseq_info_matrix` para `size_t` (que o `to_1D(int, ...)` trunca)
sobe as instruções em 14% e os ciclos em 24%. Com contador `int` (o nosso),
o `unsigned` não custa nada, o que explica o `index64` sem efeito. Variantes
`ours_cplusplus.size_t_counter` e `rust_cpp.index32_int` para a próxima
rodada. Consequência para o texto: a recomendação não é "índice `unsigned`",
e sim não misturar contador de 64 bits com índice de 32.

**P8 (2026-11-21, local, `perf`).** No C# → C++ (kroA100), o `index32` sobe
as instruções em cerca de 10% em todas as buscas, não num laço só
(`search_reinsertion` 7,3 → 8,0 bilhões, `update_subseq_info_matrix` 2,5 →
3,0). O índice ali é `(unsigned)(i * (d + 1) * 3 + j * 3 + k)`: o
deslocamento do campo (`k`, de T, C e W) fica dentro do truncamento, e cada
campo precisa do próprio cálculo de endereço. Variante de controle
`csharp_cpp.index32_field` (`3 * (size_t)(unsigned)(i * (d + 1) + j) + k`,
o `k` fora, como os campos da `struct` do nosso C++): 13,89 bilhões de
instruções, contra 12,94 na base e 14,41 no `index32`, então explica cerca
de um terço do acréscimo. Os ciclos locais não serviram (máquina ocupada);
o tempo sai da próxima rodada.

**P7 (2026-11-21, local).** O lin318 lento é do servidor: na máquina local
(Node 20.16, Intel de 4 núcleos), o nosso JavaScript leva 198,7 s no gr202 e
448,2 s no lin318 (razão 2,3; a tradução, 212,3 e 474,7 s), contra 7,8 no
servidor (Node 24.18, Xeon Silver 4316). O caminho de busca é o de
referência nos dois programas e nas duas instâncias (108.601 e 110.584
sorteios, o tamanho do vetor `rnd`), não há laço de desotimização
(`--trace-deopt`: 9 a 14 desotimizações, todas de aquecimento, nas duas
instâncias), e o perfil por função é o mesmo (`search_reinsertion` 62%,
`update_subseq_info_matrix` 14%, `search_two_opt` 12%, `search_swap` 10%).
IPC e faltas de cache e de TLB por instrução também são parecidos nas duas
instâncias. Falta separar versão do Node e CPU: `diag/node_versions.sh` roda
o nosso JavaScript nas duas instâncias com Node 20, 22 e 24 no servidor.

## Hipóteses novas da Fase 3 (achados durante a escrita das variantes)

Pedido do autor (2026-11-21): registrar aqui, como hipóteses, as diferenças
novas que aparecerem ao ler o código, mesmo fora da lista da Fase 2.

### Java
- **H-J1a. `Z` não constante no nosso `getSeq`** (refutada, ver Fase 4). Em `tSolution.java`,
  `private static int Z = 3` não é `final`; o C2 só dobra constantes de
  campos `static final`, então todo acesso às subsequências lê `Z` da
  memória e multiplica. Variante `ours_java.z_final`.
- **H-J5. Aquecimento do JIT.** O gap se concentra nas instâncias pequenas
  (berlin52 0,88, kroA100 0,80) e cai nas grandes (gr202 0,95, lin318
  0,93). O nosso Java tem mais classes e chamadas no caminho quente
  (`tArray`, `tSolution`, `tInfo`, acessores), o que pode atrasar a
  compilação e a inlining do C2; nas execuções curtas, esse atraso pesa
  mais. Teste possível sem mudar código: contar e cronometrar as
  compilações (`-Xlog:jit+compilation`) nos dois programas, ou medir o
  tempo de cada busca local. Fica para depois das variantes J1 a J4.
- **H-J6. JIT e acessores:** ver H-J5; as variantes J1 e J1a dizem se o
  custo é dos acessores em si ou do `Z` lido da memória.
- **Sem efeito esperado, só registro:** a tradução aloca um vetor novo em
  cada `reinsert` (o nosso usa um buffer fixo), o que joga contra ela; e
  aceita o melhor movimento com `cost_best < cost - DBL_EPSILON` no swap e
  na reinserção, enquanto o nosso compara sem ε (custos idênticos nas
  quatro instâncias).
- O `index` de `update_subseq_info_matrix` na tradução (herdado do nosso
  C++) não muda o laço: a recomputação é completa nos dois lados.

### Swift
- **H-S3b. Copy-on-write do vetor de subsequências.** No nosso Swift,
  `localBest = current` e `var s = s` no `perturb` fazem as duas soluções
  dividirem o mesmo `seq`; a escrita seguinte (`rvnd` ou
  `updateSubseqMatrix`) copia o vetor inteiro (n² estruturas de 24 bytes)
  a cada perturbação e a cada melhora. A tradução copia só a rota e o
  custo (`solutionCpy`). É o S3 da lista, com o mecanismo identificado.
  Variante `ours_swift.route_copy`. A tradução Swift → C++ faz o mesmo em
  C++ com cópia profunda (`perturb(Solution s)` por valor e
  `localBest = current`): variante `swift_cpp.route_copy`.
- O nosso Swift guarda a matriz de distâncias em `[[Double]]` (linhas
  separadas); a Swift → C++ guarda em vetor plano (herdado de uma versão
  anterior do nosso Swift).

### JavaScript
- **H-JS3. `Array(n)` grande em modo dicionário.** O V8 cria `Array(n)` com
  elementos esparsos (modo dicionário) quando n passa de um limite
  (cerca de 100 mil). O `seq` do nosso JavaScript tem 3(n+1)² posições:
  31 mil no kroA100, 124 mil no gr202 e 305 mil no lin318. Se o array
  não voltar ao modo rápido ao ser preenchido, os acessos ficam muito
  mais lentos nas instâncias grandes. A variante `float64array` testa isso
  junto com o tipo dos elementos. Não explica o lin318 lento nos dois
  lados (a tradução usa `Float64Array`). **Refutada** (lote principal): o
  `Float64Array` deixa o nosso 4 a 19% mais lento, e não mais rápido.
- **H-JS4. `Float64Array` mais lento que `Array` de `double`.** Um `Array`
  preenchido só com números fica em modo de elementos `double` sem caixa
  no V8; o acesso a `Float64Array` passa pelas verificações de array
  tipado. Explicaria o gap do C++ → JavaScript. Teste possível: a inversa,
  as três subsequências da tradução como `Array`.
- **H-JS5. lin318 lento nos dois lados.** O lin318 leva 7,8 vezes o gr202
  no nosso JavaScript e na tradução (C++: 3,1). Como os dois usam layouts
  diferentes (`Array` intercalado e três `Float64Array`), o efeito não vem
  do vetor de subsequências. Não investigado (decisão do autor).

### C++ (índice)
- **H-R1b. O custo do índice de 32 bits depende do contexto.** O nosso C++
  calcula `to_1D` em `int` e devolve `unsigned`, o mesmo truncamento das
  variantes `index32`, e trocar por `size_t` não muda nada; nas traduções,
  o truncamento custa de 1% (Swift → C++) a 36% (Rust → C++). Hipótese: o
  `unsigned` só atrapalha quando o compilador precisa refazer a extensão de
  zero do índice dentro de um laço em que, com 64 bits, ele transformaria a
  multiplicação em incremento de ponteiro (indexação plana `row * n + col`
  com `n` lido de membro); no nosso, o índice é recalculado em cada acesso
  de todo jeito. Teste possível: ler o assembly dos laços do `rust_cpp`
  com e sem o `index32`, e uma variante `index32` sobre o `row_matrix`.

### To-C++
- **Fortran → C++:** a matriz de distâncias fica em colunas (column-major)
  e os índices partem de 1, com conversões para `size_t` em cada acesso;
  além do ε (variantes `epsilon`, `index32`, `row_matrix`).
- **Go → C++:** a matriz fica em `vector<vector<double>>` e a solução é
  acessada por ponteiro (`solut->at(i, j)`), sem ε.
- **C# → C++:** subsequências e matriz em vetores planos com índice `int`
  calculado em cada acesso (`GetSeq`, `GetCost`).

## Restrições

- Servidor q-dual (`ssh q-dual`, via `jump`, com fail2ban): logins únicos e
  cuidadosos. CPUs 0 a 6 e 40 a 46 ocupados por outro usuário.
- Não alterar o código das nossas implementações neste experimento (só esta
  pasta). Commits com assunto `<escopo>: <resumo>`, sem corpo, sem push sem
  pedido do autor.
- Os testes rodam no servidor, não na máquina local (que só monta e confere).

## Problemas e soluções

Registro de cada problema encontrado no caminho e da solução adotada.

1. **`bench.py` fora do repositório (2026-11-20).** O `bench.py` e o `core.py`
   usados nas rodadas do Julia estão só no servidor (`~/mlp-test/mlp/benchmarks/`,
   sem histórico no git); vieram do `script.zip` de um coautor. **Solução:**
   copiados sem alteração para `harness/` (sha256 `5a6447eb…e468` do
   `bench.py`, `a5ae0445…ac15` do `core.py`). O `core.py` deduz a raiz do
   repositório pela posição do arquivo, então o `run.sh` copia os dois para
   `<repo>/benchmarks/` (pasta não versionada) antes de rodar. O catálogo do
   harness antigo foi copiado sem alteração para `harness/catalog_base.py`
   (sha256 `2ce5e558…7001`, de
   `doc/artigo-linguagens/results/translations/claude/benchmark_claude_3instancias.py`).
2. **Servidor sem as toolchains (2026-11-20).** A conta no q-dual só tem
   `g++` 11.4, Java 11 sem `javac` e o Julia de `~/opt`; faltam `gfortran`,
   `go`, `javac`, `node`, `swift`, `cargo`, e não há Nix. O podman (3.4.4,
   sem root) funciona. **Solução:** rodar tudo num container
   `docker.io/nixos/nix` (digest `sha256:7462805…d5470`) com
   `nix develop` do `flake.nix`/`flake.lock` deste repositório, que declara
   todas as toolchains (as mesmas que o repositório usa); o `/nix` fica num
   volume (`fa-nix`) para não baixar tudo de novo. O `taskset -c 8,9,10,11` é
   aplicado ao `podman run`, e o `bench.py` herda essa afinidade (ele fixa
   cada execução num dos núcleos permitidos). O log de carga roda fora do
   container. Versões registradas em `runs/<data>/versions.txt`.
3. **Impressões e relógios das traduções (2026-11-20).** Nove das onze
   traduções imprimem o custo (e, em algumas, a solução) dentro da função
   cronometrada; o `swift_cpp` imprime também o progresso em toda iteração de
   `Imax`; o `rust_cpp` só mede tempo de CPU; o `swift_cpp` só imprime `TIME:`;
   as traduções para Java, Swift, Go, Fortran e JavaScript medem o tempo de
   relógio em milissegundos inteiros. **Solução:** um `measure.patch` por
   tradução (dez ao todo; o `fortran_cpp` não precisa): a função da busca
   passa a devolver a melhor solução (ou o custo), e o `main` imprime depois
   de parar o relógio; o tempo de relógio passa a ser calculado em
   nanossegundos; o `rust_cpp` ganha um `steady_clock` e a linha
   `wall clock (s):`. No `swift_cpp`, a saída da busca é capturada em
   memória durante a busca e impressa depois (a formatação dessas linhas, de
   microssegundos, fica dentro do tempo). Nada no código da busca muda. As
   onze bases compilaram e reproduziram o custo do berlin52 (143721) na
   máquina local.
4. **Go sem `go.mod` (2026-11-20).** A tradução para Go (`cpp_go/go_src`) não
   é um módulo, e o `go build` atual recusa o pacote. **Solução:** o catálogo
   gerado pelo `make_variants.py` define `GO111MODULE=off` no ambiente de
   build.
5. **Resolução de milissegundos em três implementações nossas (2026-11-20).**
   O nosso C++ (`duration_cast<milliseconds>`), o nosso JavaScript (`Date`) e
   o nosso Fortran (`system_clock` com a taxa padrão) medem o tempo de relógio
   em milissegundos. Este experimento não altera o código das nossas
   implementações. **Solução:** aceito e registrado; o efeito é de no máximo
   0,5 ms por execução, perto de 1% só no berlin52 do C++ (cerca de 65 ms);
   a análise dá também o gap sem o berlin52 para conferir.
6. **Teste de build parou sem mensagem (2026-11-20).** No primeiro teste de
   todas as bases, o `run.sh` parou dentro do container ao gravar as
   versões: o `swift --version` do flake imprime antes um aviso
   (`glibc not found for 'x86_64-pc-linux-gnu'`), o `| head -1` fechou o
   pipe, o `swift` recebeu SIGPIPE e o `set -o pipefail` encerrou o script.
   O Swift em si funciona (versão 5.10.1, a do `nixpkgs-stable` do flake; o
   `swift build` do nosso Swift compila). **Solução:** o bloco de versões
   roda com `pipefail` desligado, e o aviso é filtrado (`2db8ec3`).
7. **`rust_cpp` sem binário no caminho de execução (2026-11-20).** O catálogo
   antigo compila a tradução Rust → C++ em `./solve`, mas a executa como
   `./cpp/solve`; no harness antigo isso só funcionava porque a pasta da
   tradução guardava um `cpp/solve` compilado antes, que não copiamos por
   ser binário. **Solução:** o catálogo gerado pelo `make_variants.py` faz o
   build das traduções em C++ gravar o binário no caminho de execução.
8. **Interrupção pelo limite de uso da API (2026-11-20 a 21).** O subagente
   parou no meio da Fase 1, enquanto a rodada de base rodava no servidor.
   **Solução:** a rodada seguiu sozinha (lançada com `nohup`) e terminou
   inteira (720 execuções, `FIM`); ao retomar, conferi o servidor e o
   repositório, trouxe os resultados para `runs/base-20261120/` e continuei
   da Fase 1. Nada foi perdido nem repetido.
9. **Verificação local do L1 (JIT do Julia) inconclusiva (2026-11-20).** Medir
   o tempo de compilação por fora (`Base.cumulative_compile_time_ns` em volta
   de `include`) não isola o trecho cronometrado: inclui o carregamento, e a
   tradução só executa o `main` quando roda como programa principal.
   **Solução:** fica pendente, como sonda dentro do trecho cronometrado
   (variante de medição na Fase 3).
10. **Variantes sobre as nossas implementações (2026-11-21).** O
    `make_variants.py` só copiava traduções. **Solução:** base
    `ours/<pasta>` no manifesto copia a pasta da nossa implementação (sem
    binários) e usa o build e a execução de `REFERENCIAS` do catálogo; o
    `analyze.py` compara `ours_<pasta>.<fator>` com `<pasta>`.
11. **Patches gerados por script (2026-11-21).** As trocas repetitivas
    (acessores, cópias locais, índices) foram feitas por expressões
    regulares sobre uma cópia, e o patch sai do `diff -ruN`. Os scripts
    ficam em `variants/scripts/` (`WORK=<pasta> python gen_<ling>.py`
    regenera os patches); cada script confere quantas trocas fez e falha se
    sobrar um acesso antigo.
12. **Teste local sem Swift nem Nix (2026-11-21).** Esta máquina não tem
    Swift nem Nix, e o Python do sistema não tem `pandas`/`psutil` (exigidos
    pelo `core.py`). **Solução:** venv só com `pandas` e `psutil` no
    scratchpad para o `bench.py` local (`-n 1`, berlin52); as três
    variantes do nosso Swift são conferidas no servidor (o build falha antes
    de qualquer execução se houver erro, e o custo é conferido pelo
    `analyze.py`).
13. **R3 descartado (2026-11-21).** Trocar a rota do Rust → C++ para
    `vector<int>` não compila com uma troca só: o `size_t` aparece nas
    assinaturas de todas as funções. Efeito antigo de 1 a 3%; fora.
14. **Patch errado pego no teste (2026-11-21).** No `swift_cpp.route_copy`,
    a troca que tirava o `return s` do `perturb` casou também com o fim de
    `construction`, que termina igual; o programa abortou ("stack
    smashing"). **Solução:** a troca passou a valer só a partir de
    `void perturb(`; o teste local reproduz o custo.
15. **`c/solve` modificado barrou a rodada (2026-11-22).** O binário do C é
    versionado; a compilação da rodada do P17 o reescreveu no checkout do
    servidor, e o `run.sh` seguinte recusou rodar ("uncommitted changes in
    the implementations"). **Solução:** `git checkout -- c/solve` antes de
    cada `run.sh` no script encadeado (`~/mlp-test/fa/p19.sh`).

16. **`c/solve` entrou num commit (2026-11-22).** Um `git commit -a` levou o
    binário do C recompilado localmente (`3d3d4db`+1). **Solução:** commit
    que restaura o binário anterior; commits daqui em diante só com
    `git add` dos arquivos.
17. **Variante com flag mediria o Node padrão (2026-11-22).** O comando de
    execução das variantes vem do catálogo (`catalog_base.py`: `node
    main.js`), não do `run_node.sh`; o patch do `run_node.sh` da variante
    `ours_javascript.no_concurrent_recompilation` não teria efeito, e a
    rodada das 66 instâncias com a flag (encadeada no `p19.sh`) mediria a
    configuração padrão de novo. Pego antes de a rodada começar.
    **Solução:** o `make_variants.py` passa a usar o comando do `run_*.sh`
    quando um patch da variante altera ou cria esse script (nenhuma variante
    antiga mexe nele, conferido); a cadeia `p19.sh` foi parada depois da
    rodada das instâncias ruidosas (que continuou) e trocada pela `p20.sh`
    (par do JavaScript com a flag nas 4 instâncias, depois o nosso nas 66).
    Teste local: os três programas com a flag dão o custo certo no berlin52
    e o catálogo gerado roda `node --no-concurrent-recompilation main.js`.

## Flag do JavaScript no repositório (autor, 2026-11-22)

O autor pediu a flag no script que roda o JavaScript e a informação no
artigo. Feito:
- `javascript/run_node.sh`: `node --no-concurrent-recompilation main.js`
  (commit `10e381f`); `harness/catalog_base.py`: `REF_NODE` com a flag, então
  o nosso JavaScript e todas as traduções para JavaScript rodam com ela.
- As variantes `*.no_concurrent_recompilation` saíram (a flag virou o
  padrão); o `make_variants.py` continua usando o `run_*.sh` alterado por
  patch, se aparecer outro caso.
- Servidor: a cadeia `p20.sh` foi parada antes de começar e trocada pela
  `p21.sh`: depois das instâncias ruidosas, o par (`javascript`,
  `cpp_javascript.base`, `cpp_javascript.array`) nas 4 instâncias
  (`runs/jsncr-20261122/`) e o nosso JavaScript nas 66
  (`runs/js66ncr-20261122/`), ambos já com a flag do repositório.
- Artigo (`ffee501`): parágrafo na Seção 6, depois do protocolo, com o
  mecanismo e os números do P16 e do P18, e um `\todo` dizendo que as
  figuras da seção ainda são da medição sem a flag; no protocolo da 7.1, as
  traduções rodam "with the same build and runtime settings as our own
  (Table compilation-modes and, for JavaScript, the flag described in
  Section 6)".
- Consequência para a decisão "os dois" (P18): a configuração principal das
  figuras passa a ser com a flag (`js66ncr`); a medição sem a flag (P15) fica
  como análise à parte, citada no parágrafo.

## Pendências futuras (registradas em 2026-11-22)

- **F1.** Medir de novo os pares fora do escopo (from-C++ Rust, C#, C, Lua,
  LuaJIT, Python, PyPy; to-C++ C, Python, Java, JavaScript, Julia) com o
  protocolo novo, para tirar as linhas antigas das tabelas (decisão 1c):
  escrever as variantes, cerca de um dia, e 3 a 4 h de servidor.
- **F2.** Seção 6 com os dados novos: trocar as instâncias ruidosas
  (`runs/noisy4-20261122/`) e o JavaScript (Node 24 padrão, P15; com a flag,
  `runs/js66ncr-20261122/`) nos dados das figuras, refazer as Figuras 4 a 6 e
  os mapas de calor e rever a frase "JavaScript is the only language whose
  gap grows with instance size" (P13, P15, P18, P19). Avisar o autor antes de
  mexer no artigo.
- **F3.** Figura 7: unificar o método de tempo (pendência antiga) e decidir
  se ela é medida de novo com o protocolo da análise nova.
- **F4.** Hipótese H-JS7: no teste local, os três programas JavaScript com a
  flag gastaram mais tempo de CPU no berlin52 (1,35 a 1,75 s) que o nosso
  sem a flag (0,61 s, outra hora). Conferir nas instâncias curtas da rodada
  `jsncr` contra o P15; se a flag pesa nelas, dizer no texto.
- **F7.** `mlp-standard` e `mlp-naive` ainda rodam `node main.js` sem a
  flag. Perguntar ao autor se a flag vale para as três (é configuração do
  runtime, não otimização do código).
- **F5.** Espelhar na monografia (`doc/tcc`) as tabelas e o texto da Seção 7
  depois que o artigo fechar.
- **F6.** Pendências antigas que continuam: bump do submódulo no repositório
  do TCC; referência do Nix na bibliografia; campanha completa nova (menos
  Julia), se o F2 não bastar; métrica dos mapas de calor (`\todo` do
  apêndice).

