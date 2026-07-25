enum infoIndex: Int {
    case T = 0, C, W
}

struct SubseqInfo {
	var t = 0.0
	var c = 0.0
	var w = 0.0
}

struct Solution {
	var sequence = [Int]()
	var cost = 0.0

	// Subseq matrix
	var seq = [[[Double]]]()
	var dimension = 0

	// This will copy the subseq info
	// possible variations would be implementing CoW for SubseqInfo
	// or making it a class, benchmark is needed
	func getSeq(_ i: Int, _ j: Int, _ k: Int) -> Double {
		return seq[i][j][k]
	}

	mutating func updateSubseqMatrix(data: Data) {
		dimension = sequence.count

		if seq.count == 0 {
			seq = Array(repeating: Array(repeating: Array(repeating:0, count: 3), count: dimension), count: dimension)
		}

		for i in 0..<dimension {
			let k = 1 - i - (i != 0 ? 0 : 1)

			seq[i][i][infoIndex.T.rawValue] = 0.0
			seq[i][i][infoIndex.C.rawValue] = 0.0
			seq[i][i][infoIndex.W.rawValue] = (i != 0 ? 1.0 : 0.0)

			for j in i+1..<dimension {
				let jPrev = j-1
				seq[i][j][infoIndex.T.rawValue] = data.getDistance(sequence[jPrev], sequence[j]) + self.getSeq(i , jPrev, infoIndex.T.rawValue)
				seq[i][j][infoIndex.C.rawValue] = self.getSeq(i,j,infoIndex.T.rawValue) + self.getSeq(i, jPrev, infoIndex.C.rawValue)
				seq[i][j][infoIndex.W.rawValue] = Double(j+k)
			}
		}

		cost = self.getSeq(0, dimension-1, infoIndex.C.rawValue)
	}

}
