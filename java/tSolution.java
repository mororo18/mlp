import java.util.ArrayList;

class tSolution {
    private ArrayList<Integer> s;
    private double [][][] seq;
    private double cost;


    tSolution(int dimen, double c) {
        seq = new double [dimen+1][dimen+1][3];
        cost = c;
    }

    void storeSolut(ArrayList<Integer> sN) {s = sN;}
    double getSeq(int i, int j, int k) {return seq[i][j][k];}
    void setSeq(int i, int j, int k, double value) {seq[i][j][k] = value;}

    double getCost() {return cost;}
    void setCost(double value) {cost = value;}

    ArrayList<Integer> getSolut() {return s;}
    ArrayList<Integer> getSolutCpy() {return new ArrayList<>(s);}
}
