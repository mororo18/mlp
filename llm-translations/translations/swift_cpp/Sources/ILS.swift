import Foundation

private let R = [0.00,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09,0.10,0.11,0.12,0.13,0.14,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25]
private let moves = [Move.swap, Move.twoOpt, Move.orOpt(1), Move.orOpt(2), Move.orOpt(3)]
private let EPSILON = 1e-16

enum Move {
	case swap
	case twoOpt
	case orOpt(Int)
}

extension ContiguousArray {

	mutating func reinsert(blockStart: Int, blockSize: Int, at: Int) {
		if at > blockStart {
	        self.reverse(blockStart, at)
	        self.reverse(blockStart, at-blockSize)
	        self.reverse(at-blockSize+1, at)
		} else {
	        self.reverse(at+1, blockStart+blockSize)
	        self.reverse(at+1, at+1+blockSize)
	        self.reverse(at+1+blockSize, blockStart+blockSize)
		}
	}
    
    mutating func reverse(_ i: Int, _ j: Int) {
		var i = i
		var j = j
        while i < j {
            self.swapAt(i, j)
            i += 1
            j -= 1
        }
    }
}

func sort(arr: inout ContiguousArray<Int>, r: Int, data: Data){
	quickSort(&arr, 0, arr.count-1, data, r)
}

func quickSort(_ arr: inout ContiguousArray<Int>, _ left: Int, _ right: Int, _ data: Data, _ r: Int) {

	if left >= right {
		return
	}

	let pivot = partition(&arr, left, right, data, r)
	quickSort(&arr, left, pivot-1, data, r)
	quickSort(&arr, pivot+1, right, data, r)

}

func partition(_ arr: inout ContiguousArray<Int>, _ left: Int, _ right: Int, _ data: Data, _ r: Int) -> Int {
	let pivot = arr[right]
	var i = left-1

	for j in left..<right {
		if data.getDistance(r, arr[j]) < data.getDistance(r, pivot) {
			i += 1
			arr.swapAt(i, j)
		}
	}

	arr.swapAt(i+1, right)
	return i + 1
}

func construction(data: inout Data) -> Solution {

	var s = Solution()
	s.sequence = [0]

	let _ = data.rndCrnt() // Mock alpha rnd choice
	var cl = ContiguousArray(1..<data.dimension)
	var r = 0
	while !cl.isEmpty {
		sort(arr: &cl, r: r, data: data)
		let index = data.rndCrnt()

		let c = cl[index]
		r = c
		s.sequence.append(c)
		cl.remove(at: index)
	}
	s.sequence.append(0)
	s.updateSubseqMatrix(data: data)

	return s
}

func perturb(s: Solution, data: inout Data) -> Solution {
	var s = s
	var aStart = 1;
    var aEnd = 1;
    var bStart = 1;
    var bEnd = 1;

    while (aStart <= bStart &&  bStart <= aEnd) || (bStart <= aStart && aStart <= bEnd) {
        aStart = data.rndCrnt();
        aEnd = aStart + data.rndCrnt();

        bStart = data.rndCrnt();
        bEnd = bStart + data.rndCrnt();
    }

	if (bStart < aStart) {
		(aStart, bStart) = (bStart, aStart)
		(aEnd, bEnd) = (bEnd, aEnd)
	}

	let aBlockSize = aEnd - aStart
	let bBlockSize = bEnd - bStart


	s.sequence.reinsert(blockStart: bStart, blockSize: bBlockSize,at: aEnd-1)
	s.sequence.reinsert(blockStart: aStart, blockSize: aBlockSize, at: bEnd-1)

	s.updateSubseqMatrix(data: data)

	return s
}

func searchSwap(s: inout Solution, data: Data) -> Bool {

	var costBest = Double.infinity
	var bestI = -1
	var bestJ = -1

	for i in 1..<data.dimension-1 {
		let iPrev = i-1
		let iNext = i+1
		let sIPrev = s.sequence[iPrev]
		let sI = s.sequence[i]
		let sINext = s.sequence[iNext]

		let costConcat1 =               s.getSeq(0, iPrev).t + data.getDistance(sIPrev, sINext)
		let costConcat2 = costConcat1 + s.getSeq(i, iNext).t + data.getDistance(sI, s.sequence[i+2])

		let costNew = s.getSeq(0, iPrev).c
					+ s.getSeq(i, iNext).w                * costConcat1 + data.getDistance(sINext, sI)
					+ s.getSeq(iNext+1, data.dimension).w * costConcat2 + s.getSeq(iNext+1, data.dimension).c

		if costNew < costBest {
			costBest = costNew - EPSILON
			bestI = i
			bestJ = iNext
		}
		
		for j in iNext+1..<data.dimension {
			let jNext = j+1
			let jPrev = j-1
			let sJPrev = s.sequence[jPrev]
			let sJ = s.sequence[j]
			let sJNext = s.sequence[jNext]

			let costConcat1 = s.getSeq(0, iPrev).t + data.getDistance(sIPrev, sJ)
			let costConcat2 = costConcat1 + data.getDistance(sJ, sINext)
			let costConcat3 = costConcat2 + s.getSeq(iNext, jPrev).t + data.getDistance(sJPrev, sI)
			let costConcat4 = costConcat3 + data.getDistance(sI, sJNext)

			let costNew = s.getSeq(0, iPrev).c
						+ costConcat1
						+ s.getSeq(iNext, jPrev).w * costConcat2 + s.getSeq(iNext, jPrev).c
						+ costConcat3
						+ s.getSeq(jNext, data.dimension).w * costConcat4 + s.getSeq(jNext, data.dimension).c

			if costNew < costBest {
				costBest = costNew - EPSILON
				bestI = i
				bestJ = j
			}

		}
	}

	if costBest < s.cost {
		s.sequence.swapAt(bestI, bestJ)
		s.updateSubseqMatrix(data: data)
		return true
	}
	
	return false
}

func searchTwoOpt(s: inout Solution, data: Data) -> Bool {

	var costBest = Double.infinity
	var bestI = -1
	var bestJ = -1

	for i in 1..<data.dimension-1 {
		let iPrev = i-1
		let sIPrev = s.sequence[iPrev]
		let sI = s.sequence[i]
		var revSeqCost = s.getSeq(i, i+1).t

		for j in i+2..<data.dimension {
			let jNext = j+1
			let sJPrev = s.sequence[j-1]
			let sJ = s.sequence[j]
			let sJNext = s.sequence[jNext]

			revSeqCost += data.getDistance(sJPrev, sJ) * (s.getSeq(i, j).w - 1.0)

			let costConcat1 = s.getSeq(0, iPrev).t + data.getDistance(sJ, sIPrev)
			let costConcat2 = costConcat1 + s.getSeq(i, j).t + data.getDistance(sJNext, sI)

			let costNew = s.getSeq(0, iPrev).c
						+ s.getSeq(i, j).w * costConcat1 + revSeqCost
						+ s.getSeq(jNext, data.dimension).w * costConcat2 + s.getSeq(jNext, data.dimension).c

			if costNew < costBest {
				costBest = costNew - EPSILON
				bestI = i
				bestJ = j
			}

		}
	}

	if costBest < s.cost {
		s.sequence.reverse(bestI, bestJ)
		s.updateSubseqMatrix(data: data)
		return true
	}
	
	return false
}

func searchOrOpt(s: inout Solution, blockSize: Int, data: Data) -> Bool {

    var costBest = Double.infinity;
	var bestI = -1
	var bestPos = -1

	for i in 1...data.dimension-blockSize {
        let j = blockSize + i - 1
        let iPrev = i-1
        let jNext = j+1
        let sI = s.sequence[i]
        let sJ = s.sequence[j]
        let sIPrev = s.sequence[iPrev]
        let sJNext = s.sequence[jNext]

		for k in 0..<iPrev {
            let kNext = k+1
            let sK = s.sequence[k]
            let sKNext = s.sequence[kNext]
            
            let costConcat1 = s.getSeq(0, k).t + data.getDistance(sK, sI)
            let costConcat2 = costConcat1 + s.getSeq(i, j).t + data.getDistance(sJ, sKNext)
            let costConcat3 = costConcat2 + s.getSeq(kNext, iPrev).t + data.getDistance(sIPrev, sJNext)

            let costNew = s.getSeq(0, k).c
                         + s.getSeq(i, j).w                  * costConcat1 + s.getSeq(i, j).c
                         + s.getSeq(kNext, iPrev).w          * costConcat2 + s.getSeq(kNext, iPrev).c
                         + s.getSeq(jNext, data.dimension).w * costConcat3 + s.getSeq(jNext, data.dimension).c

            if(costNew < costBest){
                costBest = costNew - EPSILON
                bestI = i
                bestPos = k
            }
        }

		for k in i+blockSize..<data.dimension {
            let kNext = k+1
            let sK = s.sequence[k]
            let sKNext = s.sequence[kNext]

            let costConcat1 = s.getSeq(0, iPrev).t + data.getDistance(sIPrev, sJNext)
            let costConcat2 = costConcat1 + s.getSeq(jNext, k).t + data.getDistance(sK, sI)
            let costConcat3 = costConcat2 + s.getSeq(i, j).t + data.getDistance(sJ, sKNext)

            let costNew = s.getSeq(0, iPrev).c
                    + s.getSeq(jNext, k).w              * costConcat1 + s.getSeq(jNext, k).c
                    + s.getSeq(i, j).w                  * costConcat2 + s.getSeq(i, j).c
                    + s.getSeq(kNext, data.dimension).w * costConcat3 + s.getSeq(kNext, data.dimension).c

            if costNew < costBest {
                costBest = costNew - EPSILON
                bestI = i
                bestPos = k
            }
        }
    }

    if costBest < s.cost {
        s.sequence.reinsert(blockStart: bestI, blockSize: blockSize, at: bestPos)
        s.updateSubseqMatrix(data: data)
        return true
    }

	return false
}

func rvnd(s: inout Solution, data: inout Data) {
	
	var nl = moves

	while !nl.isEmpty {
		let index = data.rndCrnt()
		let move = nl[index]
		var improved = false
		switch move {
			case .swap: improved = searchSwap(s: &s, data: data)
			case .twoOpt: improved = searchTwoOpt(s: &s, data: data)
			case .orOpt(let blockSize): improved = searchOrOpt(s: &s, blockSize: blockSize, data: data)
		}

		if improved {
			nl = moves
			continue
		}

		nl.remove(at: index)
	}
}

func gils_rvnd(iMax: Int, iIls: Int, data: inout Data) {
	var bestOfAll = Solution()
	bestOfAll.cost = Double.infinity

	for i in 0..<iMax {
		print("[+] Local Search", i);

		var localBest = construction(data: &data)
		var current = localBest

		print("\t[+] Constructing Initial Solution.. ");
		print("\t\tSequence:", current.sequence)
		print("\t\tCost:", current.cost)

		print("\t[+] Looking for the best Neighbor..")
		var iterIls = 0
		while iterIls < iIls {

			rvnd(s: &current, data: &data)
		
			if current.cost < localBest.cost {
				localBest = current
				iterIls = 0
			}

			current = perturb(s: localBest, data: &data)

			iterIls += 1
		}

		if localBest.cost < bestOfAll.cost {
			bestOfAll = localBest
		}
		print("\tCurrent Best Cost", bestOfAll.cost)
	}

	print("COST:", bestOfAll.cost)
}
