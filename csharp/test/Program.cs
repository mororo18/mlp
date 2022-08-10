
using System.Diagnostics;
using System;

namespace Test {
    class test {
        const int Z = 3;
        const int dimension = 2000;

        static int to_1D(int i, int j, int k) {
            return Z* (i*dimension +j) + k;
        }

        static void Main(){
            double total = 0;
            Random rand = new Random();
            int [] index = new int[dimension];
            long time;

            const float D = 1e9f;

            // random index-es array
            for (int i = 0; i < dimension; i++) {
                index[i] = rand.Next(dimension);
            }

            {
            double [,] matrix = new double [dimension, dimension];
            double [][] jagg2d = new double [dimension][];
            double [] flat2d = new double [dimension * dimension];


            for(int i = 0; i < dimension; i++){
                jagg2d[i] = new double [dimension];
            }

            // RAND ==> Teste com acesso aleatorio de matriz buscando o "pior caso";
            //          -> Notar idexacao invertida;
            //
            //
            Console.Write("Teste com acesso ");
#if RAND
            Console.Write("ALEATORIO");
#else
            Console.Write("SEQUENCIAL");
#endif
            Console.Write(" Array2D\n");

            /*===============================================================*/
            time = -Stopwatch.GetTimestamp();
            for(int i = 0 ; i < dimension; i++){
                for(int j = 0; j < dimension; j++){
#if RAND
                    total += matrix[index[j], index[i]];
#else
                    total += matrix[i, j];
#endif
                }
            }
            matrix = null;
            time += Stopwatch.GetTimestamp();
            Console.WriteLine("MatrixStd\t" + time/D);
            /*===============================================================*/



            /*===============================================================*/
            time = -Stopwatch.GetTimestamp();
            for(int i = 0 ; i < dimension; i++){
                for(int j = 0; j < dimension; j++){
#if RAND
                    total += jagg2d[index[j]][index[i]];
#else
                    total += jagg2d[i][j];
#endif
                }
            }
            jagg2d = null;

            time += Stopwatch.GetTimestamp();
            Console.WriteLine("MatrixJagg\t" + time/D);
            /*===============================================================*/



            /*===============================================================*/
            time = -Stopwatch.GetTimestamp();
            for(int i = 0 ; i < dimension; i++){
                for(int j = 0; j < dimension; j++){
#if RAND
                    total += flat2d[index[j] *  dimension + index[i]];
#else
                    total += flat2d[i *  dimension + j];
#endif
                }
            }
            flat2d = null;
            time += Stopwatch.GetTimestamp();
            Console.WriteLine("MatrixFlat\t" + time/D);
            /*===============================================================*/
            }


            Console.Write("Teste com acesso ");
#if RAND
            Console.Write("ALEATORIO");
#else
            Console.Write("SEQUENCIAL");
#endif
            Console.Write(" Array3D\n");


            
            double [,,] arr3d = new double [dimension, dimension, Z];
            double [][][] jagg3d = new double [dimension][][];
            double [] flat3d = new double [dimension*dimension*Z];
            double [][,] test = new double [dimension][,];

            for(int i = 0; i < dimension; i++){
                jagg3d[i] = new double [dimension][];
                test[i] = new double[dimension, Z];
                for(int j = 0; j < dimension; j++){
                    jagg3d[i][j] = new double [Z];
                }
            }

            time = -Stopwatch.GetTimestamp();
            for(int i = 0 ; i < dimension; i++){
                for(int j = 0; j < dimension; j++){
                    for(int k = 0; k < Z; k++){
#if RAND
                        total += arr3d[index[j], index[i], k];
#else
                        total += arr3d[i,j,k];
#endif
                    }
                }
            }
            time += Stopwatch.GetTimestamp();
            Console.WriteLine("Array3D\t\t" + time/D);

            time = -Stopwatch.GetTimestamp();
            for(int i = 0 ; i < dimension; i++){
                for(int j = 0; j < dimension; j++){
                    for(int k = 0; k < Z; k++){
#if RAND
                        total += jagg3d[index[j]][index[i]][k];
#else
                        total += jagg3d[i][j][k];
#endif
                    }
                }
            }
            time += Stopwatch.GetTimestamp();
            Console.WriteLine("Jagged3D\t" + time/D);

            time = -Stopwatch.GetTimestamp();
            for(int i = 0 ; i < dimension; i++){
                for(int j = 0; j < dimension; j++){
                    for(int k = 0; k < Z; k++){
#if RAND
                        total += flat3d[to_1D(index[j], index[i], k)];
#else
                        total += flat3d[to_1D(i, j, k)];
#endif
                    }
                }
            }
            time += Stopwatch.GetTimestamp();
            Console.WriteLine("Flat3D\t\t" + time/D);

            time = -Stopwatch.GetTimestamp();
            for(int i = 0 ; i < dimension; i++){
                for(int j = 0; j < dimension; j++){
                    for(int k = 0; k < Z; k++){
#if RAND
                        total += test[index[j]][index[i], k];
#else
                        total += test[i][j, k];
#endif
                    }
                }
            }
            time += Stopwatch.GetTimestamp();
            Console.WriteLine("MatrixFlat\t" + time/D);


        }
    }
}
