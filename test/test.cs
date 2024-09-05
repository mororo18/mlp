
using System.Diagnostics;
using System;

namespace Test {
    class test {
        static void Main(){
            int dimension = 100;
            double total = 0;

            double [][] jagged = new double [dimension][];
            for(int i = 0; i < dimension; i++){
                jagged[i] = new double [dimension];
            }

            double [,] multi = new double [dimension, dimension];

            long m = -Stopwatch.GetTimestamp();
            for(int k = 0; k < dimension*dimension*dimension; k++){
                for(int i = 0 ; i < dimension; i++){
                    for(int j = 0; j < dimension; j++){
                        multi[i, j] = i*j;
                    }
                }
            }
            m += Stopwatch.GetTimestamp();

            long ja = -Stopwatch.GetTimestamp();
            for(int k = 0; k < dimension*dimension*dimension; k++){
                for(int i = 0 ; i < dimension; i++){
                    for(int j = 0; j < dimension; j++){
                        jagged[i][j] = i*j;
                    }
                }
            }

            ja += Stopwatch.GetTimestamp();
            Console.WriteLine("WRITE");
            Console.WriteLine("Jagged " + ja/10e6);
            Console.WriteLine("Multi " + m/10e6);

            m = 0;
            m = -Stopwatch.GetTimestamp();
            for(int k = 0; k < dimension*dimension*dimension; k++){
                for(int i = 0 ; i < dimension; i++){
                    for(int j = 0; j < dimension; j++){
                        total += multi[i, j];// = i*j;
                    }
                }
            }
            m += Stopwatch.GetTimestamp();


            total = 0;
            ja = 0;
            ja = -Stopwatch.GetTimestamp();
            for(int k = 0; k < dimension*dimension*dimension; k++){
                for(int i = 0 ; i < dimension; i++){
                    for(int j = 0; j < dimension; j++){
                        total += jagged[i][j];// = i*j;
                    }
                }
            }

            ja += Stopwatch.GetTimestamp();

            Console.WriteLine("READ");
            Console.WriteLine("Jagged " + ja/10e6);
            Console.WriteLine("Multi " + m/10e6);
        }
    }
}
