# Python — Implementation Notes

Decision log for `mlp-main/python`. Format per `code_optimization.md`
step 7. Only `main` is optimized (see root `CLAUDE.md`) — this file
documents investigation and decisions for this branch only.

## Contexto: dois interpretadores, dois arquivos

`main` já mantém dois arquivos de entrada desde antes desta investigação:
`main.py` (usado por `run_python3.sh`, interpretador CPython) e
`main_pypy.py` (usado por `run_pypy.sh`, PyPy). Isso é relevante porque a
pergunta "qual estrutura de dados `seq` é melhor?" não tem uma resposta
única — depende de qual interpretador está rodando.

## Achado inicial que motivou a investigação (2026-07-11)

Comparação `main.py` vs `mlp-standard/python/main.py`: as duas usam a
mesma estrutura aninhada `List[List[List[float]]]` — nenhuma diferença
estrutural entre as variantes pro caminho CPython. Levantado como
possível bug (seção 6b do `doc/tcc/EXPERIMENTO.md`): "Python em `main`
nunca teve, de fato, uma estrutura otimizada". Um rewrite inicial pra
achatar `seq` (lista 1D com indexação manual, replicando o padrão já
usado por JavaScript em `main`) foi escrito e validado por corretude
(COST batendo), mas **revertido antes de aplicar** — a pedido do autor,
pra seguir o mesmo processo empírico usado na otimização do Lua
(medir → perfilar → pesquisar → microbenchmark → só então aplicar) em
vez de assumir que "achatar" é sempre a escolha certa.

## Achado que já existia e não tinha sido visto: `main_pypy.py`

Ao investigar, descoberto que `main_pypy.py` **já usa uma estrutura
achatada** (`seq` como lista 1D, indexada via `to_1D(i, j, k, size)`) —
exatamente a técnica que o rewrite revertido acima também usava,
desenvolvida independentemente em algum momento anterior a esta sessão.
O arquivo tem bastante código morto comentado e um docstring alertando
sobre "movimento avaliado de maneira incorreta" (aparentemente um bug já
resolvido, não reproduzido nos testes desta investigação). Validado
como correto: `COST` bate com o valor esperado em `burma14` (20315) e
`att48` (209320), rodando via `pypy3 main_pypy.py`.

Ou seja: **antes desta investigação**, `main` já tinha, sem
documentação, uma estrutura otimizada só pro caminho PyPy — o caminho
CPython (`main.py`) nunca recebeu o mesmo tratamento, o que pareceu
(incorretamente, ver abaixo) uma lacuna a corrigir.

## Profiling (2026-07-12)

`python3 -m cProfile -s cumulative main.py`, instância `att48` (48
nós), `main.py` original (lista aninhada):

| função | tottime | cumtime | % do total |
|---|---|---|---|
| `search_reinsertion` | 11.89s | 14.94s | 69% |
| `update_subseq_info_matrix` | 5.85s | 5.85s | 27% |
| `search_swap` | 1.94s | 2.98s | 14% |
| `search_two_opt` | 1.73s | 2.90s | 13% |

Mesmo padrão do Lua: as funções que tocam `seq` pesadamente dominam o
tempo total (~96% combinado). Confirma que a estrutura de dados de
`seq` é o alvo certo de investigação — não é preciso procurar overhead
em outro lugar.

## Pesquisa externa (2026-07-12)

Buscado material sobre lista aninhada vs. lista achatada vs.
`array.array` em Python — a maioria dos resultados genéricos de busca
foca em desempenho de *achatamento* (a operação de transformar uma
estrutura, ex. `itertools.chain`) e em NumPy, não no padrão de acesso
repetido que o GILS-RVND faz. O achado mais relevante veio de testar
diretamente (ver microbenchmark abaixo) em vez de generalizar de
fontes sobre outros workloads: números de int/float em CPython são
objetos alocados no heap ("boxed"), e a maior parte dos inteiros usados
numa indexação `i*size+j` (tipicamente > 256) fica fora do cache de
inteiros pequenos do CPython — cada operação aritmética de índice cria
um novo objeto `PyLong`. Isso não é uma preocupação em PyPy, cujo JIT
consegue especializar essas operações pra inteiros de máquina sem
boxing, uma vez que o trace é compilado.

## Microbenchmark (2026-07-12): nested vs. flat (via função) vs. flat (inline) vs. array.array

Metodologia: acesso "windowed" (simula o padrão de `search_reinsertion`
— índice base aleatório + janela de até 30 posições, não acesso
totalmente aleatório), em escala **leve** (n=48, como `att48`,
~56 KB) e **intensiva** (n=299, como `pr299`, ~2.1 MB — maior que o L2
mas ainda dentro do L3 de 6 MB desta máquina). Scripts descartados após
uso (não commitados — resultado documentado aqui). Somatório de
checagem (`chk`) idêntico em todas as variantes em cada rodada,
confirmando que todas fazem o mesmo trabalho.

**CPython, n=48, 2.000.000 repetições:**

| estrutura | tempo | relativo |
|---|---|---|
| **nested (atual)** | **8.53s** | **1.00x** |
| flat-list (chamada `to_1D()`) | 36.36s | 4.26x mais lento |
| flat-array (chamada `to_1D()`) | 45.46s | 5.33x mais lento |
| flat-list (aritmética inline) | 30.30s | 3.55x mais lento |
| flat-array (aritmética inline) | 36.79s | 4.31x mais lento |

**CPython, n=299, 1.000.000 repetições:**

| estrutura | tempo | relativo |
|---|---|---|
| **nested (atual)** | **7.44s** | **1.00x** |
| flat-list (chamada) | 27.32s | 3.67x mais lento |
| flat-array (chamada) | 30.99s | 4.16x mais lento |
| flat-list (inline) | 21.62s | 2.90x mais lento |
| flat-array (inline) | 27.29s | 3.67x mais lento |

**PyPy, n=48, 2.000.000 repetições:**

| estrutura | tempo | relativo |
|---|---|---|
| nested | 0.972s | 2.68x mais lento |
| flat-list (chamada) | 0.490s | 1.35x mais lento |
| flat-array (chamada) | 1.628s | 4.48x mais lento |
| **flat-list (inline)** | **0.363s** | **1.00x** |
| flat-array (inline) | 1.010s | 2.78x mais lento |

**PyPy, n=299, 1.000.000 repetições:**

| estrutura | tempo | relativo |
|---|---|---|
| nested | 1.741s | 6.72x mais lento |
| **flat-list (chamada)** | **0.259s** | **1.00x** |
| flat-array (chamada) | 1.207s | 4.66x mais lento |
| flat-list (inline) | 0.425s | 1.64x mais lento |
| flat-array (inline) | 0.738s | 2.85x mais lento |

**Conclusões do microbenchmark:**
1. **Achatar a estrutura piora CPython em 3-5x**, nas duas escalas —
   consistente, não é ruído. A causa raiz (confirmada testando a versão
   inline, que remove o overhead de chamada de função e ainda assim
   fica 2.9-3.6x mais lenta que nested) é o custo de alocar novos
   `PyLong` a cada aritmética de índice, não o overhead de chamada de
   função em si (embora a chamada também contribua — call é sempre
   pior que inline dentro de cada família).
2. **Achatar melhora PyPy em 2.7-6.7x** — o JIT elimina o custo de
   boxing que domina em CPython, e ganha adicionalmente por evitar 3
   níveis de indireção de ponteiro (nested) em favor de 1 acesso
   direto.
3. **`array.array` nunca venceu**, em nenhuma combinação — o custo de
   empacotar/desempacotar `double` a cada acesso supera qualquer
   ganho de armazenamento compacto, tanto em CPython quanto em PyPy
   (PyPy já otimiza bem listas nativas de float via suas "list
   strategies" internas; não precisa de `array.array` pra evitar
   boxing).
4. A conclusão **não muda entre escala leve e intensiva** — não é uma
   questão de caber ou não em cache, é uma questão de custo por
   operação em cada interpretador.

## Validação em escala real (2026-07-12)

Confirmando que o microbenchmark isolado bate com o comportamento do
algoritmo completo, instância `att48`:

| arquivo | interpretador | TIME |
|---|---|---|
| `main.py` (nested) | python3 | 21.85s |
| `main_pypy.py` (flat) | python3 | **timeout em 2min** (>5.5x mais lento, confirma a regressão) |
| `main.py` (nested) | pypy3 | 2.68s |
| `main_pypy.py` (flat) | pypy3 | **2.20s** (~18% mais rápido — ganho real, menor que o microbenchmark isolado porque o algoritmo completo inclui trabalho que não toca `seq`, ex. `construction`/`sort`) |

## Conclusão: nenhuma mudança de código necessária

`main.py` (nested, CPython) e `main_pypy.py` (flat, PyPy) **já estão,
cada um, na estrutura correta pro seu interpretador** — confirmado
empiricamente em duas escalas e em escala real de algoritmo completo,
não só suposto. O achado original ("Python não tem estrutura otimizada")
estava certo sobre o *fato* (as duas variantes CPython são idênticas)
mas errado sobre a *interpretação* (não é uma lacuna — é porque a
estrutura ótima pra CPython é a mesma que a estrutura genérica de
`standard`, ao contrário de C/C++/Rust/etc., onde achatar sempre ajuda).

**Consequência pro texto do TCC**: Python é um caso genuinamente
diferente das outras 10 linguagens — não há ganho de reestruturação de
dados disponível pro CPython nesse workload, o que é em si um
resultado interessante pra RQ3 da ESTRUTURA.md ("as técnicas de
otimização transferem entre linguagens?" — resposta parcial: não,
nem sempre, e o motivo mecanístico é rastreável ao custo de boxing de
inteiros do CPython). Vale registrar isso explicitamente na tabela de
metodologia como "verificado, sem alteração" em vez de deixar como
item em aberto.

## Diferença de API do main_pypy.py: revisado e confirmado inócuo

Nota original desta seção dizia isso como algo não verificado ("cosmético,
não afeta corretude nem performance"). Revisão pedida — em vez de só
reformular o texto, refiz a checagem com o mesmo rigor do resto deste
documento (não aceitar nada como garantido sem validar).

**Diferença de API**: `main_pypy.py` não aceita a flag `verbose` que
`main.py` aceita (`grep` por `verbose|'-v'|sys.argv` não encontra nenhum
tratamento de argumento de linha de comando no arquivo) e sempre imprime
progresso de forma incondicional dentro de `GILS_RVND` — 4 prints por
iteração de `Imax` (10 iterações) mais 2 no final (`COST`, `SOLUTION`),
totalizando ~42 chamadas a `print()` por execução.

**Por que isso não afeta o harness**: `run_bm.py` nunca passa `-v`/
`--verbose` a nenhuma linguagem (`grep -n "verbose" run_bm.py` sem
match) e `run_pypy.sh` chama `pypy3 main_pypy.py` sem argumentos — logo
a ausência da flag é irrelevante na prática: o harness nunca exercita
modo verbose em nenhuma linguagem. O parsing de `get_COST`/`get_TIME`
usa `l.find("COST")`/`l.find("TIME")` (case-sensitive, primeira linha
que casa) — o único print de `main_pypy.py` fora do padrão oficial que
se aproxima é `"\tCurrent best solution cost: {}"`, com "cost" em
minúsculo, que não colide com a busca por `"COST"`. Não há colisão de
substring com nenhum dos ~42 prints ativos.

**Por que os prints incondicionais não são um viés de medição**: como
esses prints ficam dentro de `GILS_RVND`, que é a região cronometrada
(`start = time.time()` antes da chamada, `time.time() - start` depois),
uma hipótese razoável era que eles inflassem o `TIME:` reportado do
PyPy de forma sistemática (viés de metodologia, não só estética) — vale
medir em vez de assumir. Isolei o custo real dos ~42 `print()`/`format()`
idênticos aos do arquivo, com saída redirecionada a arquivo (mesmo
padrão de `subprocess.Popen(..., stdout=tempf)` que `run_bm.py` usa),
em duas escalas:

| Tamanho de `SOLUTION` | Overhead medido (5x pypy3) |
|---|---|
| 300 elementos (~pr299) | 0.0025 – 0.0029 s |
| 18512 elementos (maior instância do dataset, d18512) | 0.0031 – 0.0032 s |

Overhead de ~2.5-3.2 milissegundos por execução, independente da escala
da instância (dominado pelo custo fixo das ~42 chamadas, não pelo
tamanho da lista impressa). Comparado a runtimes reais medidos nesta
sessão (~1.2-1.6s em att48/48 nós, ~422s em pr299/299 nós), isso é
<0.25% mesmo no caso mais rápido — na prática desprezível, não um viés
que mude conclusões de desempenho. **Tentativa inicial de medir isso via
comparação end-to-end (N=8 em pr299, com/sem prints) foi abandonada por
ser desproporcional**: uma única execução em pr299 sob PyPy leva ~7 min,
tornando N=16 execuções (~112 min) caro demais pra uma pergunta que a
medição isolada acima já resolve com muito mais precisão e em segundos.

**Código morto**: os blocos de `print()` "ativos" à primeira vista nas
linhas ~354-356 (dentro de `search_reinsertion`) e ~539-545 (fim do
arquivo, prints de `t_swap`/`t_reinsertion`/`t_or_opt2`/`t_or_opt3`/
`t_two_opt`/`t_seq`) estão de fato dentro de blocos `"""..."""`
(string literal solta, não docstring formal) — ou seja, são
comprovadamente inertes (nunca executados), não apenas comentados
visualmente com `#`. Confirmado lendo o código-fonte diretamente.

**Conclusão**: a caracterização original ("cosmético, não afeta
corretude nem performance") estava correta, mas não tinha sido validada
— agora está, com números. Nenhuma mudança de código é necessária em
`main_pypy.py`.
