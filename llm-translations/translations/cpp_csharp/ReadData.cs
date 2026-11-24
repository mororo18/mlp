using System;
using System.Globalization;
using System.IO;

namespace MlpSolver
{
    // Port of readData.cpp (TSPLIB instance reader).
    public static class ReadDataFile
    {
        // Token reader, equivalent to the whitespace-delimited "in >> x" of std::ifstream.
        private class TokenReader
        {
            private readonly StreamReader reader;

            public TokenReader(StreamReader reader)
            {
                this.reader = reader;
            }

            public string ReadToken()
            {
                int c;
                do
                {
                    c = reader.Read();
                    if (c < 0) throw new EndOfStreamException("arquivo terminou inesperadamente");
                } while (char.IsWhiteSpace((char)c));

                string token = "";
                while (c >= 0 && !char.IsWhiteSpace((char)c))
                {
                    token += (char)c;
                    c = reader.Read();
                }

                return token;
            }

            public int ReadInt()
            {
                return int.Parse(ReadToken(), CultureInfo.InvariantCulture);
            }

            public double ReadDouble()
            {
                return double.Parse(ReadToken(), NumberStyles.Float, CultureInfo.InvariantCulture);
            }
        }

        public static void ReadDataInstance(string[] argv, out int Dimension, out double[][] Mdist)
        {
            if (argv.Length < 1)
            {
                Console.WriteLine("\nFaltando parametros\n");
                Console.WriteLine(" ./exec [Instancia] ");
                Environment.Exit(1);
            }

            if (argv.Length > 1)
            {
                Console.WriteLine("\nMuitos parametros\n");
                Console.WriteLine(" ./exec [Instancia] ");
                Environment.Exit(1);
            }

            int N;
            string arquivo = "", ewt;

            string instancia = argv[0];

            if (!File.Exists(instancia))
            {
                Console.WriteLine("arquivo nao pode ser aberto\n");
                Environment.Exit(1);
            }

            StreamReader file = new StreamReader(instancia);
            TokenReader inp = new TokenReader(file);

            while (arquivo != "DIMENSION:" && arquivo != "DIMENSION")
            {
                arquivo = inp.ReadToken();
            }

            if (arquivo == "DIMENSION") arquivo = inp.ReadToken();

            N = inp.ReadInt();

            while (arquivo != "EDGE_WEIGHT_TYPE:" && arquivo != "EDGE_WEIGHT_TYPE")
            {
                arquivo = inp.ReadToken();
            }
            if (arquivo == "EDGE_WEIGHT_TYPE") arquivo = inp.ReadToken();

            ewt = inp.ReadToken();

            double[] x = new double[N + 1];
            double[] y = new double[N + 1];

            // Alocar matriz 2D
            double[][] dist = new double[N + 1][];

            for (int i = 0; i < N + 1; i++)
            {
                dist[i] = new double[N + 1];
            }

            if (ewt == "EXPLICIT")
            {
                while (arquivo != "EDGE_WEIGHT_FORMAT:" && arquivo != "EDGE_WEIGHT_FORMAT")
                {
                    arquivo = inp.ReadToken();
                }

                string ewf;
                if (arquivo == "EDGE_WEIGHT_FORMAT") arquivo = inp.ReadToken();
                ewf = inp.ReadToken();

                if (ewf == "FUNCTION")
                {
                    Console.WriteLine("FUNCTION - Nao suportado!");
                }

                else if (ewf == "FULL_MATRIX")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int i = 1; i < N + 1; i++)
                    {
                        for (int j = 1; j < N + 1; j++)
                        {
                            dist[i][j] = inp.ReadDouble();
                        }
                    }
                }

                else if (ewf == "UPPER_ROW")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int i = 1; i < N; i++)
                    {
                        for (int j = i + 1; j < N + 1; j++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }

                    for (int i = 1; i < N + 1; i++)
                    {
                        dist[i][i] = 0;
                    }
                }

                else if (ewf == "LOWER_ROW")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int i = 2; i < N + 1; i++)
                    {
                        for (int j = 1; j < i; j++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }

                    for (int i = 1; i < N + 1; i++)
                    {
                        dist[i][i] = 0;
                    }
                }

                else if (ewf == "UPPER_DIAG_ROW")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int i = 1; i < N + 1; i++)
                    {
                        for (int j = i; j < N + 1; j++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }
                }

                else if (ewf == "LOWER_DIAG_ROW")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int i = 1; i < N + 1; i++)
                    {
                        for (int j = 1; j <= i; j++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }
                }

                else if (ewf == "UPPER_COL")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int j = 2; j < N + 1; j++)
                    {
                        for (int i = 1; i < j; i++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }

                    for (int i = 1; i < N + 1; i++)
                    {
                        dist[i][i] = 0;
                    }
                }

                else if (ewf == "LOWER_COL")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int j = 1; j < N; j++)
                    {
                        for (int i = j + 1; i < N + 1; j++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }

                    for (int i = 1; i < N + 1; i++)
                    {
                        dist[i][i] = 0;
                    }
                }

                else if (ewf == "UPPER_DIAG_COL")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int j = 1; j < N + 1; j++)
                    {
                        for (int i = 1; i <= j; i++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }
                }

                else if (ewf == "LOWER_DIAG_COL")
                {
                    while (arquivo != "EDGE_WEIGHT_SECTION")
                    {
                        arquivo = inp.ReadToken();
                    }

                    // Preencher Matriz Distancia
                    for (int j = 1; j < N + 1; j++)
                    {
                        for (int i = j; i < N + 1; j++)
                        {
                            dist[i][j] = inp.ReadDouble();
                            dist[j][i] = dist[i][j];
                        }
                    }
                }
            }

            else if (ewt == "EUC_2D")
            {
                while (arquivo != "NODE_COORD_SECTION")
                {
                    arquivo = inp.ReadToken();
                }
                // ler coordenadas
                int tempCity;
                for (int i = 1; i < N + 1; i++)
                {
                    tempCity = inp.ReadInt();
                    x[i] = inp.ReadDouble();
                    y[i] = inp.ReadDouble();
                }

                // Calcular Matriz Distancia (Euclidiana)
                for (int i = 1; i < N + 1; i++)
                {
                    for (int j = 1; j < N + 1; j++)
                    {
                        dist[i][j] = Math.Floor(CalcDistEuc(x, y, i, j) + 0.5);
                    }
                }
            }

            else if (ewt == "EUD_3D")
            {
                Console.WriteLine("EUC_3D - Nao suportado!");
            }

            else if (ewt == "MAX_2D")
            {
                Console.WriteLine("MAX_2D - Nao suportado!");
            }

            else if (ewt == "MAX_3D")
            {
                Console.WriteLine("MAX_3D - Nao Suportado!");
            }

            else if (ewt == "MAN_2D")
            {
                Console.WriteLine("MAN_2D - Nao suportado!");
            }

            else if (ewt == "MAN_3D")
            {
                Console.WriteLine("MAN_3D - Nao Suportado!");
            }

            else if (ewt == "CEIL_2D")
            {
                while (arquivo != "NODE_COORD_SECTION")
                {
                    arquivo = inp.ReadToken();
                }
                // ler coordenadas
                int tempCity;
                for (int i = 1; i < N + 1; i++)
                {
                    tempCity = inp.ReadInt();
                    x[i] = inp.ReadDouble();
                    y[i] = inp.ReadDouble();
                }

                // Calcular Matriz Distancia (Euclidiana)
                for (int i = 1; i < N + 1; i++)
                {
                    for (int j = 1; j < N + 1; j++)
                    {
                        dist[i][j] = Math.Ceiling(CalcDistEuc(x, y, i, j));
                    }
                }
            }

            else if (ewt == "GEO")
            {
                while (arquivo != "NODE_COORD_SECTION")
                {
                    arquivo = inp.ReadToken();
                }
                // ler coordenadas
                int tempCity;
                for (int i = 1; i < N + 1; i++)
                {
                    tempCity = inp.ReadInt();
                    x[i] = inp.ReadDouble();
                    y[i] = inp.ReadDouble();
                }

                double[] latitude = new double[N + 1];
                double[] longitude = new double[N + 1];

                CalcLatLong(x, y, N, latitude, longitude);

                // Calcular Matriz Distancia
                for (int i = 1; i < N + 1; i++)
                {
                    for (int j = 1; j < N + 1; j++)
                    {
                        dist[i][j] = CalcDistGeo(latitude, longitude, i, j);
                    }
                }
            }

            else if (ewt == "ATT")
            {
                while (arquivo != "NODE_COORD_SECTION")
                {
                    arquivo = inp.ReadToken();
                }

                // ler coordenadas
                int tempCity;
                int[] tempX = new int[N + 1];
                int[] tempY = new int[N + 1];

                for (int i = 1; i < N + 1; i++)
                {
                    tempCity = inp.ReadInt();
                    tempX[i] = inp.ReadInt();
                    tempY[i] = inp.ReadInt();
                    x[i] = tempX[i];
                    y[i] = tempY[i];
                }

                // Calcular Matriz Distancia (Pesudo-Euclidiana)
                for (int i = 1; i < N + 1; i++)
                {
                    for (int j = 1; j < N + 1; j++)
                    {
                        dist[i][j] = CalcDistAtt(x, y, i, j);
                    }
                }
            }

            else if (ewt == "XRAY1")
            {
                Console.WriteLine("XRAY1 - Nao suportado!");
            }

            else if (ewt == "XRAY2")
            {
                Console.WriteLine("XRAY2 - Nao suportado!");
            }

            else if (ewt == "SPECIAL")
            {
                Console.WriteLine("SPECIAL - Nao suportado!");
            }

            file.Close();

            Dimension = N;
            Mdist = dist;
        }

        private static double CalcDistEuc(double[] X, double[] Y, int I, int J)
        {
            return
            Math.Sqrt(Math.Pow(X[I] - X[J], 2) + Math.Pow(Y[I] - Y[J], 2));
        }

        private static double CalcDistAtt(double[] X, double[] Y, int I, int J)
        {
            // Calcula Pseudo Distancia Euclidiana
            double rij, tij, dij;

            rij = Math.Sqrt((Math.Pow(X[I] - X[J], 2) + Math.Pow(Y[I] - Y[J], 2)) / 10);
            tij = Math.Floor(rij + 0.5);

            if (tij < rij)
                dij = tij + 1;
            else
                dij = tij;

            return dij;
        }

        private static void CalcLatLong(double[] X, double[] Y, int n, double[] latit, double[] longit)
        {
            double PI = 3.141592, min;
            int deg;

            for (int i = 1; i < n + 1; i++)
            {
                deg = (int)X[i];
                min = X[i] - deg;
                latit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
            }

            for (int i = 1; i < n + 1; i++)
            {
                deg = (int)Y[i];
                min = Y[i] - deg;
                longit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
            }
        }

        private static double CalcDistGeo(double[] latit, double[] longit, int I, int J)
        {
            double q1, q2, q3, RRR = 6378.388;

            q1 = Math.Cos(longit[I] - longit[J]);
            q2 = Math.Cos(latit[I] - latit[J]);
            q3 = Math.Cos(latit[I] + latit[J]);

            return
            (int)(RRR * Math.Acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
        }
    }
}
