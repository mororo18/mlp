using System;
using System.Diagnostics;

namespace MLP {
    class main {
        static void Main(){
            Console.WriteLine("Hello World!");
            GILS_RVND tsp = new GILS_RVND();
            
            DateTime start = DateTime.Now;
            tsp.solve();
            DateTime end = DateTime.Now;

            TimeSpan ts = (end - start);
            Console.WriteLine("TIME: "+ ts.TotalMilliseconds/10e2);
        }
    }
}
