using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;

namespace MlpSolver
{
    public static class Data
    {
        private static double[][] MatrixAllocate(int sz)
        {
            double[][] ptr = new double[sz][];

            for (int i = 0; i < sz; i++)
            {
                ptr[i] = new double[sz];
            }
            return ptr;
        }

        // Equivalent of std::stoi: parses the leading integer of the string.
        private static int Stoi(string s)
        {
            int i = 0;
            while (i < s.Length && char.IsWhiteSpace(s[i])) i++;

            int start = i;
            if (i < s.Length && (s[i] == '+' || s[i] == '-')) i++;
            while (i < s.Length && char.IsDigit(s[i])) i++;

            return int.Parse(s.Substring(start, i - start), CultureInfo.InvariantCulture);
        }

        // Equivalent of std::stod: parses the leading floating point value of the string.
        private static double Stod(string s)
        {
            int i = 0;
            while (i < s.Length && char.IsWhiteSpace(s[i])) i++;

            int start = i;
            if (i < s.Length && (s[i] == '+' || s[i] == '-')) i++;
            while (i < s.Length && char.IsDigit(s[i])) i++;
            if (i < s.Length && s[i] == '.')
            {
                i++;
                while (i < s.Length && char.IsDigit(s[i])) i++;
            }
            if (i < s.Length && (s[i] == 'e' || s[i] == 'E'))
            {
                int save = i;
                i++;
                if (i < s.Length && (s[i] == '+' || s[i] == '-')) i++;
                if (i < s.Length && char.IsDigit(s[i]))
                {
                    while (i < s.Length && char.IsDigit(s[i])) i++;
                }
                else
                {
                    i = save;
                }
            }

            return double.Parse(s.Substring(start, i - start), CultureInfo.InvariantCulture);
        }

        public static int LoadData(out double[][] matrix, List<int> rnd)
        {
            int dimension = 0;
            matrix = null;

            if (File.Exists("../distance_matrix"))
            {
                using (StreamReader file = new StreamReader("../distance_matrix"))
                {
                    string line;

                    line = file.ReadLine();
                    dimension = Stoi(line);
                    matrix = MatrixAllocate(dimension);

                    bool flag = false;
                    int rnd_size = -1;
                    int i = 0;
                    while ((line = file.ReadLine()) != null)
                    {
                        int j = i + 1;
                        while (line.IndexOf(" ", StringComparison.Ordinal) != -1)
                        {
                            int pos = line.IndexOf(" ", StringComparison.Ordinal);
                            double value = Stod(line.Substring(0, pos));
                            line = line.Substring(pos + 1);
                            matrix[i][j] = value;
                            matrix[j][i] = value;
                            j++;
                        }

                        i++;

                        if (flag)
                        {
                            if (rnd_size < 0)
                            {
                                rnd_size = Stoi(line);
                                rnd.Capacity = rnd_size;
                            }
                            else
                            {
                                rnd.Add(Stoi(line));
                            }
                        }

                        if (line == "RND") flag = true;
                    }
                }
            }

            return dimension;
        }
    }

    public class DataInstance
    {
        private int dimen;
        private double[][] cost;
        private List<int> rnd = new List<int>();

        public void Load()
        {
            dimen = Data.LoadData(out cost, rnd);
        }

        public List<int> GetRndVec()
        {
            return rnd;
        }

        public int GetDimen()
        {
            return dimen;
        }

        public double[][] GetCostPtr()
        {
            return cost;
        }
    }
}
