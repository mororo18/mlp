using System;
using System.Diagnostics;

namespace MLP {
    class main {
        static void Main(){
            Console.WriteLine("Hello World!");
            GILS_RVND tsp = new GILS_RVND();

            long s = Stopwatch.GetTimestamp();
            
            DateTime start = DateTime.Now;
            tsp.solve();
            DateTime end = DateTime.Now;

            long e = Stopwatch.GetTimestamp();

            TimeSpan ts = (end - start);
            Console.WriteLine("TIME: "+ ts.TotalMilliseconds/10e2);

            Console.WriteLine("Elapsed Time is {0} ticks", (e - s)/10e6);
        }
    }
}
