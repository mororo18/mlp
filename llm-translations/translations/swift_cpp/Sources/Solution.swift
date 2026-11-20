struct SubseqInfo {
	var t = 0.0
	var c = 0.0
	var w = 0.0
}

struct Solution {
	var sequence = ContiguousArray<Int>()
	var cost = 0.0

	// Subseq matrix
	var seq = ContiguousArray<SubseqInfo>()
	var dimension = 0

	// This will copy the subseq info
	// possible variations would be implementing CoW for SubseqInfo
	// or making it a class, benchmark is needed
	func getSeq(_ i: Int, _ j: Int) -> SubseqInfo {
		return seq[i * dimension + j]
	}

	mutating func updateSubseqMatrix(data: Data) {
		dimension = sequence.count

		if seq.count == 0 {
			seq = ContiguousArray(repeating: SubseqInfo(), count: dimension*dimension)
		}

		for i in 0..<dimension {
			let k = 1 - i - (i != 0 ? 0 : 1)

			seq[i * dimension + i].t = 0.0
			seq[i * dimension + i].c = 0.0
			seq[i * dimension + i].w = (i != 0 ? 1.0 : 0.0)

			for j in i+1..<dimension {
				let jPrev = j-1
				seq[i * dimension + j].t = data.getDistance(sequence[jPrev], sequence[j]) + self.getSeq(i , jPrev).t
				seq[i * dimension + j].c = self.getSeq(i,j).t + self.getSeq(i, jPrev).c
				seq[i * dimension + j].w = Double(j+k)
			}
		}

		cost = self.getSeq(0, dimension-1).c
	}

}
