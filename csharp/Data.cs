using System;
using static System.IO.File;

namespace MLP {
    class Data {
        private int dimension;
        private double [,] matrix;
        private int [] rnd;
        private static double [] c;
        private int rnd_index;

        public Data(){
            dimension = 0;
            rnd_index = 0;
        }

        public void loadData(){
            string [] file = ReadAllLines("../distance_matrix");
            int file_line = 0;
            dimension = Int32.Parse(file[file_line++]);
            matrix = new double [dimension, dimension];

            for(int i = 1; i < dimension; i++){
                //Console.WriteLine(file[i]);
                file_line++;
                int j = i;
                while(file[i].IndexOf(" ") != -1){
                    int index = file[i].IndexOf(" ");
                    matrix[i-1, j] = Convert.ToDouble(file[i].Substring(0, index));
                    matrix[j, i-1] = matrix[i-1, j];
                    matrix[i-1, i-1] = 0.0;
                    file[i] = file[i].Substring(index+1);
                    //Console.WriteLine(matrix[i-1, j]);
                    j++;
                }
            }
            matrix[dimension-1, dimension-1] = 0.0;
            
            file_line++;
            file_line++;
            file_line++;
            int rnd_size = Int32.Parse(file[file_line]);
            file_line++;
            rnd = new int [rnd_size];
            for (int i = 0; i < rnd_size; i++) {
                rnd[i] = Int32.Parse(file[file_line++]);
            }

            c = new double [dimension*dimension];
            for (int i = 0; i < dimension; i++) {
                for (int j = i; j < dimension; j++) {
                    c[to_1D(i,j)]  = getDistance(i,j);
                    c[to_1D(j,i)]  = getDistance(j,i);
                }
            }
        
        }

        private int to_1D(int i, int j) {
            return dimension * i + j;
        }

        public double getDistance(int i, int j) {
            return matrix[i, j];
        }

        public int [] GetRnd() {return rnd;}

        public int GetDimen() {return dimension;}
        public double GetCost(int i, int j) {
            //var f = ref c[0];
            //double ret = Unsafe.Add(ref c[0], to_1D(i, j));
            return c[to_1D(i,j)];
        }
        
        public int GetRndCrnt() {return rnd[rnd_index++];}
    }
}