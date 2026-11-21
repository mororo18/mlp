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
- [ ] Um patch por fator; manifesto; montagem e build de todas as variantes;
      conferência do custo em uma instância pequena antes de subir para o
      servidor.

### Fase 4. Rodadas
- [ ] Lotes 1 e 2 no servidor. Estimativa: umas 60 variantes compiladas ou
      JIT, de 6 a 8 horas; Lua e Python interpretados, umas 2 horas.

### Fase 5. Análise e texto
- [ ] `analyze.py`, tabelas em `results/factor_analysis/<data>/` do artigo,
      texto da 7.5 e da 7.6, `README.md` desta pasta.

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
