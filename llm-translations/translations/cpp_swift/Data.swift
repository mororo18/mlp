import Foundation

func matrixAllocate(_ sz: Int) -> [[Double]] {
    return Array(repeating: Array(repeating: 0.0, count: sz), count: sz)
}

func loadData(matrix: inout [[Double]], rnd: inout [Int]) -> Int {
    var dimension = 0

    guard let content = try? String(contentsOfFile: "../distance_matrix", encoding: .utf8) else {
        return dimension
    }

    var lines = content.split(separator: "\n", omittingEmptySubsequences: false).map(String.init)
    if let last = lines.last, last.isEmpty {
        lines.removeLast()
    }

    guard !lines.isEmpty else { return dimension }

    dimension = Int(lines[0]) ?? 0
    matrix = matrixAllocate(dimension)

    var flag = false
    var rndSize = -1
    var i = 0

    for lineIndex in 1..<lines.count {
        var line = lines[lineIndex]
        let j0 = i + 1
        var j = j0

        while let spaceIndex = line.firstIndex(of: " ") {
            let value = Double(line[line.startIndex..<spaceIndex]) ?? 0.0
            line = String(line[line.index(after: spaceIndex)...])
            matrix[i][j] = value
            matrix[j][i] = value
            j += 1
        }

        i += 1

        if flag {
            if rndSize < 0 {
                rndSize = Int(line) ?? 0
                rnd.reserveCapacity(rndSize)
            } else {
                rnd.append(Int(line) ?? 0)
            }
        }

        if line == "RND" { flag = true }
    }

    return dimension
}

class Data {
    private var dimen: Int = 0
    private var cost: [[Double]] = []
    private var rnd: [Int] = []

    func load() {
        dimen = loadData(matrix: &cost, rnd: &rnd)
    }

    func getRndVec() -> [Int] {
        return rnd
    }

    func getDimen() -> Int {
        return dimen
    }

    func getCostPtr() -> [[Double]] {
        return cost
    }
}
