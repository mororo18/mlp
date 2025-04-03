using System;
using static System.IO.File;

namespace MLP {
    class Data {
        private int dimension;
        private double [,] matrix;
        private int [] rnd;

        public Data(){
            dimension = 0;
        }

        public void loadData(){
            string [] file = ReadAllLines("../distance_matrix");
            int file_line = 0;
            dimension = Int32.Parse(file[file_line++]);
            matrix = new double [dimension, dimension];

            for(int i = 1; i < dimension; i++){
                file_line++;
                int j = i;
                while(file[i].IndexOf(" ") != -1){
                    int index = file[i].IndexOf(" ");
                    matrix[i-1, j] = Convert.ToDouble(file[i].Substring(0, index));
                    matrix[j, i-1] = matrix[i-1, j];
                    matrix[i-1, i-1] = 0.0;
                    file[i] = file[i].Substring(index+1);
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
        }

        public double getDistance(int i, int j) {
            return matrix[i, j];
        }

        public int getDimension() {
            return dimension;
        }

        public int [] GetRnd() {return rnd;}
    }
}
