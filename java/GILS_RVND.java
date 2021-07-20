import java.util.Vector;

class GILS_RVND {
    private int     dimension;
    private double  [][] c;

    private final int C = 0;
    private final int T = 1;
    private final int W = 2;

    private double [][][] subseq;

    GILS_RVND(){
        Data data = new Data();
        data.loadData();

        dimension = data.getDimension();
        c = new double [dimension][dimension];
        for(int i = 0; i < dimension; i++){
            for(int j = i; j < dimension; j++){
                c[i][j] = data.getDistance(i, j);
                c[j][i] = data.getDistance(i, j);
                System.out.print(Double.toString(c[i][j])+ " ");
            }
            System.out.println();
        }

        subseq = new double [dimension+1][dimension+1][3];
    }

    private void subseq_load(Vector<Integer> s, double [][][] seq){
        double [] arr_branch = {0.0, 1.0};
        for(int i = 0; i < dimension+1; i++){
            int k = 1 - i - (i != 0 ? 0 : 1);

            seq[i][i][T] = 0.0;
            seq[i][i][C] = 0.0;
            seq[i][i][W] = (i != 0 ? 1.0 : 0.0);

            for(int j = i+1; j < dimension+1; j++){
                int a = j-1;

                seq[i][j][T] = c[s.get(a)][s.get(j)] + seq[i][a][T];
                seq[i][j][C] = seq[i][j][T] + seq[i][a][C];
                seq[i][j][W] = j + k;

                seq[j][i][T] = seq[i][j][T];
                seq[j][i][C] = seq[i][j][C];
                seq[j][i][W] = seq[i][j][W];

                System.out.print(seq[i][j][C]);
                System.out.print(" ");
            }
            System.out.println();
        }
    }

    public void solve(){
        Vector<Integer> s = new Vector<>(dimension);

        for(int i = 0; i < dimension; i++){
            s.add(i); 
            System.out.print(s.get(i));
            System.out.print(" ");
        }
        s.add(0); 
        System.out.println();

        subseq_load(s, subseq);
    }
}
