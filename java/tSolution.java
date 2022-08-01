import java.util.ArrayList;

class tSolution {
    public static final boolean TEST = true;
    private ArrayList<Integer> s;

    private static double [][][] seq3d;
    private static double [][] seq2d;

    private double cost;
    private static int Z = 3;


    tSolution(int dimen, double c) {
if (TEST)   
        seq2d = new double [dimen+1][(dimen+1)*Z];
else        
        seq3d = new double [dimen+1][dimen+1][Z];
        cost = c;
    }

    void storeSolut(ArrayList<Integer> sN) {s = sN;}

    public static double getSeq(int i, int j, int k) {
if (TEST)   
        return seq2d[i][j*Z + k];
else        
        return seq3d[i][j][k];
    }

    public static void setSeq(int i, int j, int k, double value) {
if (TEST)   
        seq2d[i][j*Z + k] = value; 
else        
        seq3d[i][j][k] = value;
    }

    double getCost() {return cost;}
    void setCost(double value) {cost = value;}

    ArrayList<Integer> getSolut() {return s;}
    ArrayList<Integer> getSolutCpy() {return new ArrayList<>(s);}
}
