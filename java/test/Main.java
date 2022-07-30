import java.util.Random;

class Main {
    public static void main(String args[]) {
        int size = 100;

        int [][] original = new int [size][size];
        MatrixTest test = new MatrixTest(size, size);

        int [] indexInt = new int [size];
        long [] index = new long [size];
        Random rand = new Random();

      //for (int i = 0; i < size; i++) {
      //    indexInt[i] = rand.nextInt(size);
      //    index[i] = indexInt[i];
      //}

        int value;

        long time = -System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                //value = original[i][j];
                value = original[indexInt[i]][indexInt[j]];
            }
        }
        time += System.nanoTime();

        System.out.print("original\t");
        System.out.println(time);

        time = - System.nanoTime();
        for (int i = 0; i < size; i++) {
            for (int j = i+1; j < size; j++) {
                //value = test.get(i, j);
                value = test.get(indexInt[i], indexInt[j]);
            }
        }
        time += System.nanoTime();
        System.out.print("meu\t\t");
        System.out.println(time);
        System.out.println(int[].class);

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
