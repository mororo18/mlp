import java.util.Random;


class Main {
    private int line_size;
    private static boolean linear = true;

    private static final double D = 1e9;

    private static int to_1D(int i, int j, int k, int size, int Z) {
        return Z*(size * i + j) + k;
    }

    public static void main(String args[]) {
        int size = 5000;
        int value;
        long time;

        int [] index = new int [size];
        Random rand = new Random();

        for (int i = 0; i < size; i++) {
            index[i] = rand.nextInt(size);
        }


        int [][] nested2D = new int [size][size];
        int [] flat2D = new int [size*size]; 


        System.out.print("Teste com acesso ");
        if (linear) System.out.print("SEQUENCIAL");
        else System.out.print("ALEATORIO");
        System.out.print(" em Matriz (");
        System.out.print(size);
        System.out.print(" x ");
        System.out.print(size);
        System.out.println(")");



        /*====================================================*/

        time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                if (linear) value = nested2D[i][j];
                else value = nested2D[index[i]][index[j]];
            }
        }
        time += System.nanoTime();

        System.out.print("Nested 2D\t");
        System.out.println(time/ D);


        /*====================================================*/

        time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                if (linear) value = flat2D[i * size + j];
                else value = flat2D[index[i] * size + index[j]];
            }
        }
        time += System.nanoTime();

        System.out.print("Flat 2D\t\t");
        System.out.println(time/ D);

        /*====================================================*/



        /*====================================================*/
        /*====================================================*/
        /*====================================================*/



        int Z = 3;
        int [][][] nested3D = new int [size][size][Z];
        int [][] partNested3D = new int [size][size*Z];
        int [] flat3D = new int [size*size*Z]; 


        System.out.print("Teste com acesso ");
        if (linear) System.out.println("SEQUENCIAL");
        else System.out.print("ALEATORIO");
        System.out.print(" em array 3D (");
        System.out.print(size);
        System.out.print(" x ");
        System.out.print(size);
        System.out.print(" x ");
        System.out.print(Z);
        System.out.println(")");

        /*====================================================*/

        time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                for (int k = 0; k < Z; k++) {
                    if (linear) value = nested3D[i][j][k];
                    else value = nested3D[index[i]][index[j]][k];
                }
            }
        }
        time += System.nanoTime();

        System.out.print("Nested 3D\t");
        System.out.println(time/ D);


        /*====================================================*/

        time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                for (int k = 0; k < Z; k++) {
                    if (linear) value = flat3D[to_1D(i, j, k, size, Z)];
                    else value = flat3D[to_1D(index[i], index[j], k, size, Z)];
                }
            }
        }
        time += System.nanoTime();

        System.out.print("Flat 3D\t\t");
        System.out.println(time/ D);

        /*====================================================*/

        time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                for (int k = 0; k < Z; k++) {
                    if (linear) value = partNested3D[i][j*Z + k];
                    else value = partNested3D[index[i]][index[j]*Z + k];
                }
            }
        }
        time += System.nanoTime();

        System.out.print("Part. Nested 3D\t");
        System.out.println(time/ D);



        /*
        for (int i = 0; i < size; i++) {
            test.set(i, i, 0);
            for (int j = i+1; j < size; j++) {
                test.set(i, j, j);
                test.set(j, i, j);

            }
        }

        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                System.out.print(test.get(i,j));
                System.out.print("\t");
            }
            System.out.println();
        }
        */
    }

}
