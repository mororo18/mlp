import Foundation

struct SeqInfo {
    var T: Double = 0.0
    var C: Double = 0.0
    var W: Double = 0.0
}

enum Neighborhood: Int {
    case reinsertion = 1
    case orOpt2
    case orOpt3
    case swapNb
    case twoOpt
}

struct TData {
    var cost: [[Double]] = []
    var dimen: Int = 0
    var rnd: [Int] = []
    var rndIndex: Int = 0
}

struct TSolution {
    var s: [Int] = []
    var seq: [SeqInfo] = []
    var cost: Double = Double.greatestFiniteMagnitude
}

@inline(__always)
func to1D(_ x: Int, _ y: Int, _ size: Int) -> Int {
    return x * (size + 1) + y
}

func solutionInit(_ data: TData) -> TSolution {
    var solut = TSolution()
    solut.s = Array(repeating: 0, count: data.dimen + 1)
    solut.seq = Array(repeating: SeqInfo(), count: (data.dimen + 1) * (data.dimen + 1))
    solut.cost = Double.greatestFiniteMagnitude
    return solut
}

func solutionCpy(_ src: TSolution, _ tgt: inout TSolution) {
    tgt.s = src.s
    tgt.cost = src.cost
}

func rTable(_ i: Int) -> Double {
    let table: [Double] = [0.00, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12,
                            0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25]
    return table[i]
}

func printS(_ s: [Int]) {
    var out = ""
    for v in s {
        out += "\(v + 1) "
    }
    print(out)
}

func partition(_ arr: inout [Int], _ left: Int, _ right: Int, _ data: TData, _ r: Int) -> Int {
    let pivot = arr[right]
    var i = left - 1
    for j in left..<right {
        if data.cost[r][arr[j]] < data.cost[r][pivot] {
            i += 1
            arr.swapAt(i, j)
        }
    }
    arr.swapAt(i + 1, right)
    return i + 1
}

func quicksort(_ arr: inout [Int], _ left: Int, _ right: Int, _ data: TData, _ r: Int) {
    if left < right {
        let pivot = partition(&arr, left, right, data, r)
        quicksort(&arr, left, pivot - 1, data, r)
        quicksort(&arr, pivot + 1, right, data, r)
    }
}

func sortByCost(_ arr: inout [Int], _ r: Int, _ data: TData) {
    quicksort(&arr, 0, arr.count - 1, data, r)
}

func construct(_ alpha: Double, _ data: inout TData) -> [Int] {
    var s = [0]
    s.reserveCapacity(data.dimen + 1)

    var cL = Array(repeating: 0, count: data.dimen - 1)
    for i in 1..<data.dimen {
        cL[i - 1] = i
    }

    var r = 0
    while !cL.isEmpty {
        sortByCost(&cL, r, data)

        let index = data.rnd[data.rndIndex]
        data.rndIndex += 1

        let c = cL[index]
        s.append(c)
        r = c
        cL.remove(at: index)
    }

    s.append(0)

    return s
}

func swapElems(_ vec: inout [Int], _ i: Int, _ j: Int) {
    vec.swapAt(i, j)
}

func reverseRange(_ vec: inout [Int], _ i: Int, _ j: Int) {
    vec[i...j].reverse()
}

func reinsert(_ vec: inout [Int], _ i: Int, _ j: Int, _ pos: Int) {
    let seq = Array(vec[i...j])
    if pos < i {
        vec.removeSubrange(i...j)
        vec.insert(contentsOf: seq, at: pos)
    } else {
        vec.insert(contentsOf: seq, at: pos)
        vec.removeSubrange(i...j)
    }
}

func updateSubseqInfoMatrix(_ solut: inout TSolution, _ data: TData, _ index: Int = 0) {
    var from = index

    for i in 0..<(data.dimen + 1) {
        let k = 1 - i - (i == 0 ? 1 : 0)
        let t = i == from

        solut.seq[to1D(i, i, data.dimen)].T = 0.0
        solut.seq[to1D(i, i, data.dimen)].C = 0.0
        solut.seq[to1D(i, i, data.dimen)].W = (i == 0) ? 0.0 : 1.0

        if i + 1 < data.dimen + 1 {
            for j in (i + 1)..<(data.dimen + 1) {
                let jPrev = j - 1

                solut.seq[to1D(i, j, data.dimen)].T = data.cost[solut.s[jPrev]][solut.s[j]] + solut.seq[to1D(i, jPrev, data.dimen)].T
                solut.seq[to1D(i, j, data.dimen)].C = solut.seq[to1D(i, j, data.dimen)].T + solut.seq[to1D(i, jPrev, data.dimen)].C
                solut.seq[to1D(i, j, data.dimen)].W = Double(j + k)
            }
        }
        from += t ? 1 : 0
    }

    solut.cost = solut.seq[to1D(0, data.dimen, data.dimen)].C
}

func searchSwap(_ solut: inout TSolution, _ data: TData) -> Bool {
    var costBest = Double.greatestFiniteMagnitude
    var I = 0
    var J = 0

    var i = 1
    while i < data.dimen - 1 {
        let iPrev = i - 1
        let iNext = i + 1

        let costConcat1 = solut.seq[to1D(0, iPrev, data.dimen)].T + data.cost[solut.s[iPrev]][solut.s[iNext]]
        let costConcat2 = costConcat1 + solut.seq[to1D(i, iNext, data.dimen)].T + data.cost[solut.s[i]][solut.s[iNext + 1]]

        var costNew = solut.seq[to1D(0, iPrev, data.dimen)].C +
            solut.seq[to1D(i, iNext, data.dimen)].W * costConcat1 + data.cost[solut.s[iNext]][solut.s[i]] +
            solut.seq[to1D(iNext + 1, data.dimen, data.dimen)].W * costConcat2 + solut.seq[to1D(iNext + 1, data.dimen, data.dimen)].C

        if costNew < costBest {
            costBest = costNew - Double.ulpOfOne
            I = i
            J = iNext
        }

        var j = iNext + 1
        while j < data.dimen {
            let jNext = j + 1
            let jPrev = j - 1

            let cc1 = solut.seq[to1D(0, iPrev, data.dimen)].T + data.cost[solut.s[iPrev]][solut.s[j]]
            let cc2 = cc1 + data.cost[solut.s[j]][solut.s[iNext]]
            let cc3 = cc2 + solut.seq[to1D(iNext, jPrev, data.dimen)].T + data.cost[solut.s[jPrev]][solut.s[i]]
            let cc4 = cc3 + data.cost[solut.s[i]][solut.s[jNext]]

            costNew = solut.seq[to1D(0, iPrev, data.dimen)].C +
                cc1 +
                solut.seq[to1D(iNext, jPrev, data.dimen)].W * cc2 + solut.seq[to1D(iNext, jPrev, data.dimen)].C +
                cc3 +
                solut.seq[to1D(jNext, data.dimen, data.dimen)].W * cc4 + solut.seq[to1D(jNext, data.dimen, data.dimen)].C

            if costNew < costBest {
                costBest = costNew - Double.ulpOfOne
                I = i
                J = j
            }

            j += 1
        }

        i += 1
    }

    if costBest < solut.cost - Double.ulpOfOne {
        swapElems(&solut.s, I, J)
        updateSubseqInfoMatrix(&solut, data, I)
        return true
    }

    return false
}

func searchTwoOpt(_ solut: inout TSolution, _ data: TData) -> Bool {
    var costBest = Double.greatestFiniteMagnitude
    var I = 0
    var J = 0

    var i = 1
    while i < data.dimen - 1 {
        let iPrev = i - 1

        var revSeqCost = solut.seq[to1D(i, i + 1, data.dimen)].T

        var j = i + 2
        while j < data.dimen {
            let jNext = j + 1

            revSeqCost += data.cost[solut.s[j - 1]][solut.s[j]] * (solut.seq[to1D(i, j, data.dimen)].W - 1.0)

            let cc1 = solut.seq[to1D(0, iPrev, data.dimen)].T + data.cost[solut.s[j]][solut.s[iPrev]]
            let cc2 = cc1 + solut.seq[to1D(i, j, data.dimen)].T + data.cost[solut.s[jNext]][solut.s[i]]

            let costNew = solut.seq[to1D(0, iPrev, data.dimen)].C +
                solut.seq[to1D(i, j, data.dimen)].W * cc1 + revSeqCost +
                solut.seq[to1D(jNext, data.dimen, data.dimen)].W * cc2 + solut.seq[to1D(jNext, data.dimen, data.dimen)].C

            if costNew < costBest {
                costBest = costNew - Double.ulpOfOne
                I = i
                J = j
            }

            j += 1
        }

        i += 1
    }

    if costBest < solut.cost - Double.ulpOfOne {
        reverseRange(&solut.s, I, J)
        updateSubseqInfoMatrix(&solut, data)
        return true
    }

    return false
}

func searchReinsertion(_ solut: inout TSolution, _ data: TData, _ opt: Int) -> Bool {
    var costBest = Double.greatestFiniteMagnitude
    var I = 0
    var J = 0
    var POS = 0

    var i = 1
    var j = opt + i - 1
    while i < data.dimen - opt + 1 {
        let jNext = j + 1
        let iPrev = i - 1

        var k = 0
        while k < iPrev {
            let kNext = k + 1

            let cc1 = solut.seq[to1D(0, k, data.dimen)].T + data.cost[solut.s[k]][solut.s[i]]
            let cc2 = cc1 + solut.seq[to1D(i, j, data.dimen)].T + data.cost[solut.s[j]][solut.s[kNext]]
            let cc3 = cc2 + solut.seq[to1D(kNext, iPrev, data.dimen)].T + data.cost[solut.s[iPrev]][solut.s[jNext]]

            let costNew = solut.seq[to1D(0, k, data.dimen)].C +
                solut.seq[to1D(i, j, data.dimen)].W * cc1 + solut.seq[to1D(i, j, data.dimen)].C +
                solut.seq[to1D(kNext, iPrev, data.dimen)].W * cc2 + solut.seq[to1D(kNext, iPrev, data.dimen)].C +
                solut.seq[to1D(jNext, data.dimen, data.dimen)].W * cc3 + solut.seq[to1D(jNext, data.dimen, data.dimen)].C

            if costNew < costBest {
                costBest = costNew - Double.ulpOfOne
                I = i
                J = j
                POS = k
            }

            k += 1
        }

        k = i + opt
        while k < data.dimen {
            let kNext = k + 1

            let cc1 = solut.seq[to1D(0, iPrev, data.dimen)].T + data.cost[solut.s[iPrev]][solut.s[jNext]]
            let cc2 = cc1 + solut.seq[to1D(jNext, k, data.dimen)].T + data.cost[solut.s[k]][solut.s[i]]
            let cc3 = cc2 + solut.seq[to1D(i, j, data.dimen)].T + data.cost[solut.s[j]][solut.s[kNext]]

            let costNew = solut.seq[to1D(0, iPrev, data.dimen)].C +
                solut.seq[to1D(jNext, k, data.dimen)].W * cc1 + solut.seq[to1D(jNext, k, data.dimen)].C +
                solut.seq[to1D(i, j, data.dimen)].W * cc2 + solut.seq[to1D(i, j, data.dimen)].C +
                solut.seq[to1D(kNext, data.dimen, data.dimen)].W * cc3 + solut.seq[to1D(kNext, data.dimen, data.dimen)].C

            if costNew < costBest {
                costBest = costNew - Double.ulpOfOne
                I = i
                J = j
                POS = k
            }

            k += 1
        }

        i += 1
        j += 1
    }

    if costBest < solut.cost - Double.ulpOfOne {
        reinsert(&solut.s, I, J, POS + 1)
        updateSubseqInfoMatrix(&solut, data, I < POS + 1 ? I : POS + 1)
        return true
    }

    return false
}

func rvnd(_ solut: inout TSolution, _ data: inout TData) {
    var neighbdList: [Neighborhood] = [.swapNb, .twoOpt, .reinsertion, .orOpt2, .orOpt3]

    while !neighbdList.isEmpty {
        let index = data.rnd[data.rndIndex]
        data.rndIndex += 1

        let neighbd = neighbdList[index]
        var improveFlag = false

        switch neighbd {
        case .reinsertion:
            improveFlag = searchReinsertion(&solut, data, Neighborhood.reinsertion.rawValue)
        case .orOpt2:
            improveFlag = searchReinsertion(&solut, data, Neighborhood.orOpt2.rawValue)
        case .orOpt3:
            improveFlag = searchReinsertion(&solut, data, Neighborhood.orOpt3.rawValue)
        case .swapNb:
            improveFlag = searchSwap(&solut, data)
        case .twoOpt:
            improveFlag = searchTwoOpt(&solut, data)
        }

        if improveFlag {
            neighbdList = [.swapNb, .twoOpt, .reinsertion, .orOpt2, .orOpt3]
        } else {
            neighbdList.remove(at: index)
        }
    }
}

func perturb(_ solut: TSolution, _ data: inout TData) -> [Int] {
    var s = solut.s
    var aStart = 1
    var aEnd = 1
    var bStart = 1
    var bEnd = 1

    while (aStart <= bStart && bStart <= aEnd) || (bStart <= aStart && aStart <= bEnd) {
        aStart = data.rnd[data.rndIndex]; data.rndIndex += 1
        aEnd = aStart + data.rnd[data.rndIndex]; data.rndIndex += 1

        bStart = data.rnd[data.rndIndex]; data.rndIndex += 1
        bEnd = bStart + data.rnd[data.rndIndex]; data.rndIndex += 1
    }

    if aStart < bStart {
        reinsert(&s, bStart, bEnd - 1, aEnd)
        reinsert(&s, aStart, aEnd - 1, bEnd)
    } else {
        reinsert(&s, aStart, aEnd - 1, bEnd)
        reinsert(&s, bStart, bEnd - 1, aEnd)
    }

    return s
}

func gilsRvnd(_ Imax: Int, _ Iils: Int, _ data: inout TData, _ verbose: Bool) {
    var solutPartial = solutionInit(data)
    var solutCrnt = solutionInit(data)
    var solutBest = solutionInit(data)

    for i in 0..<Imax {
        let aux = data.rnd[data.rndIndex]
        data.rndIndex += 1

        let alpha = rTable(aux)

        if verbose {
            print("[+] Search \(i + 1)")
            print("\t[+] Constructing..")
        }

        solutCrnt.s = construct(alpha, &data)
        updateSubseqInfoMatrix(&solutCrnt, data)

        solutionCpy(solutCrnt, &solutPartial)
        if verbose {
            print("\t[+] Looking for the best Neighbor..")
            print(String(format: "\t    Construction Cost: %.3lf", solutPartial.cost))
        }

        var iterILS = 0
        while iterILS < Iils {
            rvnd(&solutCrnt, &data)
            if solutCrnt.cost < solutPartial.cost - Double.ulpOfOne {
                solutionCpy(solutCrnt, &solutPartial)
                iterILS = 0
            }

            solutCrnt.s = perturb(solutPartial, &data)
            updateSubseqInfoMatrix(&solutCrnt, data)
            iterILS += 1
        }

        if solutPartial.cost < solutBest.cost - Double.ulpOfOne {
            solutionCpy(solutPartial, &solutBest)
        }

        if verbose {
            print("\tCurrent best cost: \(solutBest.cost)")

            var out = "SOLUCAO: "
            for v in solutBest.s {
                out += "\(v) "
            }
            print(out)
        }
    }
    print(String(format: "COST: %.2lf", solutBest.cost))
}

// MARK: - Entry point

let argv = CommandLine.arguments
let argc = argv.count

var Imax = 10
var Iils = 0
var verbose = false

for i in 1..<argc {
    if argv[i] == "-v" || argv[i] == "--verbose" {
        verbose = true
    }
}

var data = TData()

var rnd: [Int] = []
var cost: [[Double]] = []

data.dimen = loadData(matrix: &cost, rnd: &rnd)
data.cost = cost
data.rnd = rnd
data.rndIndex = 0

_ = solutionInit(data)

srand(UInt32(clock()))

Iils = data.dimen < 100 ? data.dimen : 100

let t1 = DispatchTime.now()
let start = clock()
gilsRvnd(Imax, Iils, &data, verbose)
let cpuTime = Double(clock() - start) / Double(CLOCKS_PER_SEC)
let t2 = DispatchTime.now()
let duration = (t2.uptimeNanoseconds - t1.uptimeNanoseconds) / 1_000_000

let res = Double(duration) / 10e2
print("TIME: \(cpuTime)")
print("wall clock (s): \(res)")

print("Tamanho RND \(rnd.count)")

exit(0)
