import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public class Data {

    private int dimen;
    private double[][] cost;
    private List<Integer> rnd = new ArrayList<>();

    static double[][] matrixAllocate(int sz) {
        return new double[sz][sz];
    }

    // Reads "../distance_matrix": first line is the dimension, then the upper
    // triangle (one row per line, values separated by spaces), then the "RND"
    // marker followed by the amount of random numbers and the numbers themselves.
    static int loadData(double[][][] matrix, List<Integer> rnd) {
        int dimension = 0;

        try (BufferedReader file = new BufferedReader(new FileReader("../distance_matrix"))) {
            String line;

            line = file.readLine();
            dimension = Integer.parseInt(line.trim());
            matrix[0] = matrixAllocate(dimension);

            boolean flag = false;
            int rnd_size = -1;
            int i = 0;
            while ((line = file.readLine()) != null) {
                int j = i + 1;
                int pos;
                while ((pos = line.indexOf(' ')) != -1) {
                    double value = Double.parseDouble(line.substring(0, pos));
                    line = line.substring(pos + 1);
                    matrix[0][i][j] = value;
                    matrix[0][j][i] = value;
                    j++;
                }

                i++;

                if (flag) {

                    if (rnd_size < 0) {
                        rnd_size = Integer.parseInt(line);
                    } else {
                        rnd.add(Integer.parseInt(line));
                    }
                }

                if (line.equals("RND")) flag = true;
            }

        } catch (IOException e) {
            // file could not be opened/read
        }

        return dimension;
    }

    public void load() {
        double[][][] holder = new double[1][][];
        dimen = loadData(holder, rnd);
        cost = holder[0];
    }

    public List<Integer> getRndVec() {
        return rnd;
    }

    public int getDimen() {
        return dimen;
    }

    public double[][] getCostPtr() {
        return cost;
    }
}
