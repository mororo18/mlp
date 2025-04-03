import java.io.File;  
import java.io.FileNotFoundException; 
import java.util.Scanner; 

class Data{
    private int dimension;
    private double [][] matrix;
    private int [] rnd;
    private int rnd_size;

    Data(){
        dimension = 0;
    }

    public int      getDimension()              {return dimension;}
    public double   getDistance(int i, int j)   {return matrix[i][j];}
    public int []   getRnd()                    {return rnd;}
    public int      getRndSize()                {return rnd_size;}
    
    public void loadData(){
        String file_name = "../distance_matrix";
        try{
            File file = new File(file_name);
            Scanner f_reader = new Scanner(file);

            String line = f_reader.nextLine();

            dimension = Integer.parseInt(line);
            matrix = new double [dimension][dimension];
            for(int k = 0; k < dimension; k++) matrix[k][k] = 0.0;

            for (int i = 0; i < dimension; i++) {
                line = f_reader.nextLine();

                int j = i + 1;
                while(line.indexOf(' ') != -1){
                    String cost = line.substring(0, line.indexOf(' '));
                    line = line.substring(line.indexOf(' ')+1, line.length());

                    double cost_value = Double.parseDouble(cost);
                    matrix[i][j] = cost_value;
                    matrix[j][i] = matrix[i][j];
                    j++;
                }
            }


            line = f_reader.nextLine();
            line = f_reader.nextLine();
            line = f_reader.nextLine();

            rnd_size = Integer.parseInt(line);
            rnd = new int[rnd_size];

            for (int i = 0; i < rnd_size; i++) {
                line = f_reader.nextLine();
                rnd[i] = Integer.parseInt(line);
            }


            f_reader.close();

        }catch (FileNotFoundException e){
            System.out.println("Error while reading " + file_name + " file");
            e.printStackTrace();
        }

    }
}
