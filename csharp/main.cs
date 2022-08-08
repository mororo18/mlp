using System;
using System.Diagnostics;
using System.Threading;

namespace MLP {
    class main {
        static void Main(){
            //Thread.Sleep(1000 * 10);
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
