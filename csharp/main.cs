using System;
using System.Diagnostics;

namespace MLP {
    class main {
        static void Main(string[] args){
            bool verbose = Array.Exists(args, a => a == "-v" || a == "--verbose");

            if (verbose) {
                Console.WriteLine("Hello World!");
            }
            GILS_RVND tsp = new GILS_RVND();

            Process proc = Process.GetCurrentProcess();
            TimeSpan cpuStart = proc.TotalProcessorTime;

            DateTime start = DateTime.Now;
            tsp.solve(verbose);
            DateTime end = DateTime.Now;

            proc.Refresh();
            TimeSpan cpuTime = proc.TotalProcessorTime - cpuStart;

            TimeSpan ts = (end - start);
            Console.WriteLine("TIME: "+ cpuTime.TotalSeconds);

            Console.WriteLine("wall clock (s): "+ ts.TotalSeconds);
        }
    }
}
