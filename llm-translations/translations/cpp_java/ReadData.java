import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.StreamTokenizer;

public class ReadData {

    public int dimension;
    public double[][] dist;

    // whitespace delimited token reader, mimics C++'s "in >> x"
    private static class In {
        private final StreamTokenizer st;

        In(String file) throws IOException {
            st = new StreamTokenizer(new BufferedReader(new FileReader(file)));
            st.resetSyntax();
            st.wordChars(33, 255);
            st.whitespaceChars(0, 32);
        }

        String nextWord() {
            try {
                if (st.nextToken() == StreamTokenizer.TT_EOF) return "";
                return st.sval;
            } catch (IOException e) {
                return "";
            }
        }

        double nextDouble() {
            return Double.parseDouble(nextWord());
        }

        int nextInt() {
            return (int) Double.parseDouble(nextWord());
        }
    }

    public void readData(String[] argv) {

        if (argv.length < 1) {
            System.out.println("\nFaltando parametros\n");
            System.out.println(" ./exec [Instancia] ");
            System.exit(1);
        }

        if (argv.length > 1) {
            System.out.println("\nMuitos parametros\n");
            System.out.println(" ./exec [Instancia] ");
            System.exit(1);
        }

        int N;
        String arquivo = "", ewt;

        String instancia = argv[0];

        In in = null;
        try {
            in = new In(instancia);
        } catch (IOException e) {
            System.out.println("arquivo nao pode ser aberto");
            System.exit(1);
        }

        while (!arquivo.equals("DIMENSION:") && !arquivo.equals("DIMENSION")) {
            arquivo = in.nextWord();
        }

        if (arquivo.equals("DIMENSION")) arquivo = in.nextWord();

        N = in.nextInt();

        while (!arquivo.equals("EDGE_WEIGHT_TYPE:") && !arquivo.equals("EDGE_WEIGHT_TYPE")) {
            arquivo = in.nextWord();
        }
        if (arquivo.equals("EDGE_WEIGHT_TYPE")) arquivo = in.nextWord();

        ewt = in.nextWord();

        double[] x = new double[N + 1];
        double[] y = new double[N + 1];

        // Alocar matriz 2D
        double[][] dist = new double[N + 1][N + 1];

        if (ewt.equals("EXPLICIT")) {

            while (!arquivo.equals("EDGE_WEIGHT_FORMAT:") && !arquivo.equals("EDGE_WEIGHT_FORMAT")) {
                arquivo = in.nextWord();
            }

            String ewf;
            if (arquivo.equals("EDGE_WEIGHT_FORMAT")) arquivo = in.nextWord();
            ewf = in.nextWord();

            if (ewf.equals("FUNCTION")) {
                System.out.println("FUNCTION - Nao suportado!");
            } else if (ewf.equals("FULL_MATRIX")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int i = 1; i < N + 1; i++) {
                    for (int j = 1; j < N + 1; j++) {
                        dist[i][j] = in.nextDouble();
                    }
                }
            } else if (ewf.equals("UPPER_ROW")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int i = 1; i < N; i++) {
                    for (int j = i + 1; j < N + 1; j++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }

                for (int i = 1; i < N + 1; i++) {
                    dist[i][i] = 0;
                }

            } else if (ewf.equals("LOWER_ROW")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int i = 2; i < N + 1; i++) {
                    for (int j = 1; j < i; j++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }

                for (int i = 1; i < N + 1; i++) {
                    dist[i][i] = 0;
                }
            } else if (ewf.equals("UPPER_DIAG_ROW")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int i = 1; i < N + 1; i++) {
                    for (int j = i; j < N + 1; j++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }
            } else if (ewf.equals("LOWER_DIAG_ROW")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int i = 1; i < N + 1; i++) {
                    for (int j = 1; j <= i; j++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }
            } else if (ewf.equals("UPPER_COL")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int j = 2; j < N + 1; j++) {
                    for (int i = 1; i < j; i++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }

                for (int i = 1; i < N + 1; i++) {
                    dist[i][i] = 0;
                }

            } else if (ewf.equals("LOWER_COL")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int j = 1; j < N; j++) {
                    for (int i = j + 1; i < N + 1; j++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }

                for (int i = 1; i < N + 1; i++) {
                    dist[i][i] = 0;
                }

            } else if (ewf.equals("UPPER_DIAG_COL")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int j = 1; j < N + 1; j++) {
                    for (int i = 1; i <= j; i++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }
            } else if (ewf.equals("LOWER_DIAG_COL")) {

                while (!arquivo.equals("EDGE_WEIGHT_SECTION")) {
                    arquivo = in.nextWord();
                }

                // Preencher Matriz Distancia
                for (int j = 1; j < N + 1; j++) {
                    for (int i = j; i < N + 1; j++) {
                        dist[i][j] = in.nextDouble();
                        dist[j][i] = dist[i][j];
                    }
                }
            }

        } else if (ewt.equals("EUC_2D")) {

            while (!arquivo.equals("NODE_COORD_SECTION")) {
                arquivo = in.nextWord();
            }
            // ler coordenadas
            for (int i = 1; i < N + 1; i++) {
                in.nextInt();
                x[i] = in.nextDouble();
                y[i] = in.nextDouble();
            }

            // Calcular Matriz Distancia (Euclidiana)
            for (int i = 1; i < N + 1; i++) {
                for (int j = 1; j < N + 1; j++) {
                    dist[i][j] = Math.floor(calcDistEuc(x, y, i, j) + 0.5);
                }
            }
        } else if (ewt.equals("EUD_3D")) {
            System.out.println("EUC_3D - Nao suportado!");
        } else if (ewt.equals("MAX_2D")) {
            System.out.println("MAX_2D - Nao suportado!");
        } else if (ewt.equals("MAX_3D")) {
            System.out.println("MAX_3D - Nao Suportado!");
        } else if (ewt.equals("MAN_2D")) {
            System.out.println("MAN_2D - Nao suportado!");
        } else if (ewt.equals("MAN_3D")) {
            System.out.println("MAN_3D - Nao Suportado!");
        } else if (ewt.equals("CEIL_2D")) {

            while (!arquivo.equals("NODE_COORD_SECTION")) {
                arquivo = in.nextWord();
            }
            // ler coordenadas
            for (int i = 1; i < N + 1; i++) {
                in.nextInt();
                x[i] = in.nextDouble();
                y[i] = in.nextDouble();
            }

            // Calcular Matriz Distancia (Euclidiana)
            for (int i = 1; i < N + 1; i++) {
                for (int j = 1; j < N + 1; j++) {
                    dist[i][j] = Math.ceil(calcDistEuc(x, y, i, j));
                }
            }
        } else if (ewt.equals("GEO")) {

            while (!arquivo.equals("NODE_COORD_SECTION")) {
                arquivo = in.nextWord();
            }
            // ler coordenadas
            for (int i = 1; i < N + 1; i++) {
                in.nextInt();
                x[i] = in.nextDouble();
                y[i] = in.nextDouble();
            }

            double[] latitude = new double[N + 1];
            double[] longitude = new double[N + 1];

            calcLatLong(x, y, N, latitude, longitude);

            // Calcular Matriz Distancia
            for (int i = 1; i < N + 1; i++) {
                for (int j = 1; j < N + 1; j++) {
                    dist[i][j] = calcDistGeo(latitude, longitude, i, j);
                }
            }

        } else if (ewt.equals("ATT")) {

            while (!arquivo.equals("NODE_COORD_SECTION")) {
                arquivo = in.nextWord();
            }

            // ler coordenadas
            int[] tempX = new int[N + 1];
            int[] tempY = new int[N + 1];

            for (int i = 1; i < N + 1; i++) {
                in.nextInt();
                tempX[i] = in.nextInt();
                tempY[i] = in.nextInt();
                x[i] = tempX[i];
                y[i] = tempY[i];
            }

            // Calcular Matriz Distancia (Pesudo-Euclidiana)
            for (int i = 1; i < N + 1; i++) {
                for (int j = 1; j < N + 1; j++) {
                    dist[i][j] = calcDistAtt(x, y, i, j);
                }
            }

        } else if (ewt.equals("XRAY1")) {
            System.out.println("XRAY1 - Nao suportado!");
        } else if (ewt.equals("XRAY2")) {
            System.out.println("XRAY2 - Nao suportado!");
        } else if (ewt.equals("SPECIAL")) {
            System.out.println("SPECIAL - Nao suportado!");
        }

        this.dimension = N;
        this.dist = dist;
    }

    static double calcDistEuc(double[] X, double[] Y, int I, int J) {
        return Math.sqrt(Math.pow(X[I] - X[J], 2) + Math.pow(Y[I] - Y[J], 2));
    }

    static double calcDistAtt(double[] X, double[] Y, int I, int J) {
        // Calcula Pseudo Distancia Euclidiana
        double rij, tij, dij;

        rij = Math.sqrt((Math.pow(X[I] - X[J], 2) + Math.pow(Y[I] - Y[J], 2)) / 10);
        tij = Math.floor(rij + 0.5);

        if (tij < rij)
            dij = tij + 1;
        else
            dij = tij;

        return dij;
    }

    static void calcLatLong(double[] X, double[] Y, int n, double[] latit, double[] longit) {
        double PI = 3.141592, min;
        int deg;

        for (int i = 1; i < n + 1; i++) {
            deg = (int) X[i];
            min = X[i] - deg;
            latit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
        }

        for (int i = 1; i < n + 1; i++) {
            deg = (int) Y[i];
            min = Y[i] - deg;
            longit[i] = PI * (deg + 5.0 * min / 3.0) / 180.0;
        }
    }

    static double calcDistGeo(double[] latit, double[] longit, int I, int J) {
        double q1, q2, q3, RRR = 6378.388;

        q1 = Math.cos(longit[I] - longit[J]);
        q2 = Math.cos(latit[I] - latit[J]);
        q3 = Math.cos(latit[I] + latit[J]);

        return (int) (RRR * Math.acos(0.5 * ((1.0 + q1) * q2 - (1.0 - q1) * q3)) + 1.0);
    }
}
