import java.io.File;  
import java.io.FileNotFoundException; 
import java.util.Scanner; 

class Data{
    private int dimension;
    private double [][] matrix;

    Data(){
        dimension = 0;
    }

    public int      getDimension()              {return dimension;}
    public double   getDistance(int i, int j)   {return matrix[i][j];}
    
    //programacao orientada a spaghetti
    public void loadData(){
        String file_name = "../distance_matrix";
        try{
            File file = new File(file_name);
            Scanner f_reader = new Scanner(file);

            int i = 0;
            while(f_reader.hasNextLine()){
                String line = f_reader.nextLine();

                if(dimension == 0){
                    dimension = Integer.parseInt(line);
                    matrix = new double [dimension][dimension];
                    for(int k = 0; k < dimension; k++) matrix[k][k] = 0.0;
                    System.out.println(dimension); 
                    continue;
                }

                int j = i + 1;
                while(line.indexOf(' ') != -1){
                    String cost = line.substring(0, line.indexOf(' '));
                    line = line.substring(line.indexOf(' ')+1, line.length());

                    double cost_value = Double.parseDouble(cost);
                    matrix[i][j] = cost_value;
                    matrix[j][i] = matrix[i][j];
                    j++;
                    System.out.print(cost + " ");
                }
                System.out.println(); 
                i++;
            }
            f_reader.close();

        }catch (FileNotFoundException e){
            System.out.println("Error while reading " + file_name + " file");
            e.printStackTrace();
        }

    }
}
