var data = Data()

let clock = ContinuousClock()

let iIls = data.dimension < 100 ? data.dimension : 100

let time = clock.measure {
	gils_rvnd(iMax: 10, iIls: iIls, data: &data)
}
let timeInSeconds = Double(time.components.seconds) + Double(time.components.attoseconds) / 1e18

print("TIME:", String(format: "%.6f", timeInSeconds))
print("wall clock (s):", String(format: "%.6f", timeInSeconds))
