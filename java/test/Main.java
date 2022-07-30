import java.util.Random;


class Main {
    private int line_size;
    private static boolean linear = false;

    private static int to_1D(int i, int j, int k, int size) {
        return 3*(size * i + j) + k;
    }

    public static void main(String args[]) {
        int size = 5000;
        int a = size;

        int z_size = 3;

        int [][][] nested = new int [size][size][z_size];
        int [][] partNested = new int [size][size*z_size];
        int [] flat = new int [size*size*3]; 
        MatrixTest test = new MatrixTest(size, size);

        int [] index = new int [size];
        Random rand = new Random();

        for (int i = 0; i < size; i++) {
            index[i] = rand.nextInt(size);
        }

        System.out.print("Teste com acesso ");
        if (linear) System.out.println("LINEAR");
        else System.out.println("ALEATORIO");

        int value;

        /*====================================================*/

        long time = -System.nanoTime();
        for (int i = 0; i < a; i++) {
            for (int j = i+1; j < a; j++) {
                for (int k = 0; k < z_size; k++) {
                    if (linear) value = nested[i][j][k];
                    else value = nested[index[i]][index[j]][k];
                }
            }
        }
        time += System.nanoTime();

        System.out.print("nested\t\t");
        System.out.println(time/ 10e8);


        /*====================================================*/


        /*
        time = - System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                //value = test.get(i, j);
                value = test.get(indexInt[i], indexInt[j]);
            }
        }
        time += System.nanoTime();
        System.out.print("MatrixTest\t");
        System.out.println(time);
        */

        /*====================================================*/

        time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                for (int k = 0; k < z_size; k++) {
                    if (linear) value = flat[to_1D(i, j, k, size)];
                    else value = flat[to_1D(index[i], index[j], k, size)];
                }
            }
        }
        time += System.nanoTime();

        System.out.print("Flat\t\t");
        System.out.println(time/ 10e8);

        /*====================================================*/

        time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                for (int k = 0; k < z_size; k++) {
                    if (linear) value = partNested[i][j*z_size + k];
                    else value = partNested[index[i]][index[j]*z_size + k];
                }
            }
        }
        time += System.nanoTime();

        System.out.print("Part. Nested\t");
        System.out.println(time/ 10e8);



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
