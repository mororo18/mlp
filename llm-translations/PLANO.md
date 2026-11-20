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
- [ ] `catalog.py` e `make_variants.py`; teste local de montagem e build de
      todas as bases.

### Fase 1. Gap de base
- [ ] Rodar `<par>.base` e as nossas implementações nas instâncias do lote
      (lotes 1 e 2), no servidor.
- [ ] Gap de cada tradução; lista dos gaps que ainda existem.

### Fase 2. Mapa das diferenças
- [ ] Para cada par, listar as diferenças de projeto entre a tradução e a
      nossa implementação atual (atualiza a Tabela 7 e o mapa do apêndice). As
      seis linguagens que mudaram precisam ser feitas do zero.
- [ ] Lista de fatores a testar: Tabela 7, candidatos dos relatórios antigos
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
