let verbose = CommandLine.arguments.contains("-v") || CommandLine.arguments.contains("--verbose")

var data = Data()

let clock = ContinuousClock()

let iIls = data.dimension < 100 ? data.dimension : 100

var cost = 0.0
let time = clock.measure {
	cost = gils_rvnd(iMax: 10, iIls: iIls, data: &data, verbose: verbose)
}
let timeInSeconds = Double(time.components.seconds) + Double(time.components.attoseconds) / 1e18

print("COST:", cost)
print("TIME:", String(format: "%.6f", timeInSeconds))
print("wall clock (s):", String(format: "%.6f", timeInSeconds))
