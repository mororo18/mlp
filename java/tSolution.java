//import java.util.ArrayList;

class tSolution {
    public static final boolean TEST = true;
    private tArray s;

    private double [][][] seq3d;
    private double [][] seq2d;

    private double cost;
    private static int Z = 3;


    tSolution(int dimen, double c) {
if (TEST)   
        seq2d = new double [dimen+1][(dimen+1)*Z];
else        
        seq3d = new double [dimen+1][dimen+1][Z];
        cost = c;
    }

    void storeSolut(tArray sN) {s = sN;}

    public double getSeq(int i, int j, int k) {
if (TEST)   
        return seq2d[i][j*Z + k];
else        
        return seq3d[i][j][k];
    }

    public void setSeq(int i, int j, int k, double value) {
if (TEST)   
        seq2d[i][j*Z + k] = value; 
else        
        seq3d[i][j][k] = value;
    }

    double getCost() {return cost;}
    void setCost(double value) {cost = value;}

    tArray getSolut() {return s;}
    tArray getSolutCpy() {return s.getArrayCopy();}
}
