using System;
using static System.IO.File;

namespace MLP {
    public class Data {
        private int dimension;
        private double [,] matrix;
        public Data(){
            Console.WriteLine("Opa meu abencoado!");
            dimension = 0;
        }

        public void loadData(){
            string [] file = ReadAllLines("../distance_matrix");
            dimension = Int32.Parse(file[0]);
            matrix = new double [dimension, dimension];

            for(int i = 1; i < dimension; i++){
                Console.WriteLine(file[i]);
                int j = i;
                while(file[i].IndexOf(" ") != -1){
                    int index = file[i].IndexOf(" ");
                    matrix[i-1, j] = Convert.ToDouble(file[i].Substring(0, index));
                    matrix[j, i-1] = matrix[i-1, j];
                    file[i] = file[i].Substring(index+1);
                    //Console.WriteLine(matrix[i-1, j]);
                    j++;
                }
            }

            
            for(int i = 0; i < dimension; i++){
                for(int j = i+1; j < dimension; j++){
                    Console.Write(matrix[i, j] + " ");
                }
                Console.WriteLine();
            }
        }

        public double getDistance(int i, int j){
            return matrix[i, j];
        }

        public int getDimension(){
            return dimension;
        }
    }
}
