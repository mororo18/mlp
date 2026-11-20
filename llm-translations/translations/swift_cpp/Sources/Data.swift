import Foundation

// Check if class has better performance over struct
struct Data {
	var dimension = 0
	var matrix = ContiguousArray<Double>()
	var rnd = ContiguousArray<Int>()
	var rndIdx = 0

	func getDistance(_ i: Int, _ j: Int) -> Double {
		return matrix[i*dimension + j]
	}

	mutating func rndCrnt() -> Int {
		let value = rnd[rndIdx]
		rndIdx += 1
		return value
	}

	init() {
		let url = URL(fileURLWithPath: "../distance_matrix")
		let contents = (try? String(contentsOf: url))!
	
		let lines = contents.components(separatedBy: "\n")


		let dimension = Int(lines[0])!
		self.dimension = dimension
		self.matrix = ContiguousArray(repeating: 0, count: (dimension*dimension))
		for i in 0...dimension-2 {


			self.matrix[i*dimension + i] = 0
			var j = i+1
			var values = lines[i+1].components(separatedBy: " ")
			values.removeLast()
			for value in values {

				let converted = Double(value)!
				self.matrix[i*dimension + j] = converted
				self.matrix[j*dimension + i] = converted
				j += 1
			}
		}

		let rndCount = Int(lines[dimension+3])!
		for i in 0..<rndCount {
			let linesIdx = dimension + 4 + i
			self.rnd.append(Int(lines[linesIdx])!)
		}
	}
}

