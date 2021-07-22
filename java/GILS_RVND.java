import java.util.ArrayList;
import java.util.Random;
import java.util.Collections;

class GILS_RVND {
    private final double EPSILON = 1e-8;
    private int     dimension;
    private double  [][] c;
    private double  [][][] subseq;

    private static final int C = 0;
    private static final int T = 1;
    private static final int W = 2;

    // Neighborhoods ids
    private static final int SWAP       = 0;
    private static final int REINSERTION= 1;
    private static final int OR_OPT2    = 2;
    private static final int OR_OPT3    = 3;
    private static final int TWO_OPT    = 4;

    private boolean improv_flag;

    // parameters
    private int Iils;
    private final int Imax = 10;
    private static final double [] R = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                        0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    private final int R_size = 26;

    Random rand;

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
        Iils = (dimension < 100 ? dimension : 100);

        rand = new Random();
    }

    private void subseq_load(ArrayList<Integer> s, double [][][] seq){
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

    private ArrayList<Integer> construction(double alpha){
        ArrayList<Integer> s = new ArrayList<>();
        s.add(0);

        ArrayList<Integer> cList = new ArrayList<>();
        for(int i = 1; i < dimension; i++){
            cList.add(i);
        }

        int r = 0;
        while(!cList.isEmpty()){
            cList.sort((Integer i, Integer j) -> Double.compare (c[i][s.get(s.size()-1)], c[j][s.get(s.size()-1)]));
            //erro ao usar variavel 'r' na lambda expression. 'r' n eh final(const)
            //cList.sort((Integer i, Integer j) -> Double.compare (c[i][r], c[j][r]));

            int range = (int)(((double)cList.size()) * alpha)+1;
            System.out.println(range);
            int c = cList.get(rand.nextInt(range));
            s.add(c);
            cList.remove(Integer.valueOf(c));
            r = c;
            System.out.println(cList);
        }

        s.add(0);
        System.out.println(s);
        //System.exit(1);
        return s;
    }

    private void swap(ArrayList<Integer> s, int i, int j){
        Collections.swap(s, i, j);
        //System.out.println(s);
    }

    private void reverse(ArrayList<Integer> s, int i, int j){
        Collections.reverse(s.subList(i,j+1));
        System.out.println(s);
    }

    private void reinsert(ArrayList<Integer> s, int i, int j, int pos){
        if(i > pos){
            int sz = j-i+1;
            s.addAll(pos, s.subList(i, j+1));
            s.subList(i+sz, j+1+sz).clear();
        }else{
            s.addAll(pos, s.subList(i, j+1));
            s.subList(i, j+1).clear();
        }

        System.out.println(s);
    }

    private void search_swap(ArrayList<Integer> s, double [][][] seq){

    } 

    private void search_two_opt(ArrayList<Integer> s, double [][][] seq){

    } 

    private void search_reinsertion(ArrayList<Integer> s, double [][][] seq, int opt){

    } 

    private void RVND(ArrayList<Integer> s, double [][][] subseq){
        ArrayList<Integer> neighbd_list = new ArrayList<>() {{
            add(REINSERTION);
            add(OR_OPT2);
            add(OR_OPT3);
            add(SWAP);
            add(TWO_OPT);
        }};

        while(!neighbd_list.isEmpty()){
            
            System.out.println(s);
            //reverse(s, 2, 5);
            //reinsert(s, 2, 7, 9);
            //reinsert(s, 7, 9, 3);
            //swap(s, 2, 5);
            System.exit(1);
            int i_rand = rand.nextInt(neighbd_list.size());
            int neighbd = neighbd_list.get(i_rand);

            improv_flag = false;

            switch(neighbd){
                case REINSERTION:
                    break;
                case OR_OPT2:
                    break;
                case OR_OPT3:
                    break;
                case SWAP:
                    break;
                case TWO_OPT:
                    break;
            }

            if(improv_flag){
                neighbd_list.clear();
                neighbd_list = new ArrayList<>() {{
                    add(REINSERTION);
                    add(OR_OPT2);
                    add(OR_OPT3);
                    add(SWAP);
                    add(TWO_OPT);
                }};
            }else
                neighbd_list.remove(i_rand);
        }
    }

    private void perturb(ArrayList<Integer> s, ArrayList<Integer> sl){
        
    }

    public void solve(){
        double cost_best = Double.MAX_VALUE;
        ArrayList<Integer> s_best;
        for(int i = 0; i < Imax; i++){
            double alpha = R[rand.nextInt(R_size)];
            System.out.println(alpha);

            ArrayList<Integer> s = construction(alpha);
            ArrayList<Integer> sl = new ArrayList<>(s); 
            subseq_load(s, subseq); 
            double rvnd_cost_best = subseq[0][dimension][C] - EPSILON;
            double rvnd_cost_crnt;

            int iterILS = 0;
            while(iterILS < Iils){
                RVND(s, subseq);
                rvnd_cost_crnt = subseq[0][dimension][C] - EPSILON; 
                if(rvnd_cost_crnt < rvnd_cost_best){
                    rvnd_cost_best = rvnd_cost_crnt;
                    sl.clear();
                    sl = new ArrayList<>(s);
                    iterILS = 0;
                }
                perturb(s, sl);
                subseq_load(s, subseq);
                iterILS++;
            }
            subseq_load(s, subseq);
            double sl_cost = subseq[0][dimension][C] - EPSILON;

            if(sl_cost < cost_best){
                cost_best = sl_cost;
                s_best = sl;
            }
        }

        //ArrayList<Integer> s = new ArrayList<>(dimension);
        //subseq_load(s, subseq);
    }
}
