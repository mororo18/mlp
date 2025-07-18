import java.util.ArrayList;
import java.util.Random;
import java.util.Collections;
import java.lang.Math;

class GILS_RVND {
    private int     []  rnd;
    private int     rnd_index;
    private final double EPSILON = 1e-16;
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
    private int                     Iils;
    private final int               Imax = 10;
    private static final double []  R = {0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12, 
                                        0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.20, 0.21, 0.22, 0.23, 0.24, 0.25};
    private final int               R_size = 26;

    private static long t_construction = 0;
    private static long t_swap = 0;
    private static long t_two_opt = 0;
    private static long t_or_opt2 = 0;
    private static long t_or_opt3 = 0;
    private static long t_reinsertion = 0;
    private static long t_subseq = 0;

    GILS_RVND(){
        Data data = new Data();
        data.loadData();

        dimension = data.getDimension();
        c = new double [dimension][dimension];
        for(int i = 0; i < dimension; i++){
            for(int j = i; j < dimension; j++){
                c[i][j] = data.getDistance(i, j);
                c[j][i] = data.getDistance(i, j);
            }
        }
        System.out.println(Integer.toString(data.getRndSize()));
        rnd = data.getRnd();
        rnd_index = 0;

        subseq = new double [dimension+1][dimension+1][3];
        Iils = (dimension < 100 ? dimension : 100);
    }

    private void subseq_load(ArrayList<Integer> s, double [][][] seq){
        for(int i = 0; i < dimension+1; i++){
            int k = 1 - i - (i != 0 ? 0 : 1);

            seq[i][i][T] = 0.0;
            seq[i][i][C] = 0.0;
            seq[i][i][W] = (i != 0 ? 1.0 : 0.0);

            for(int j = i+1; j < dimension+1; j++){
                int j_prev = j-1;

                seq[i][j][T] = c[s.get(j_prev)][s.get(j)] + seq[i][j_prev][T];
                seq[i][j][C] = seq[i][j][T] + seq[i][j_prev][C];
                seq[i][j][W] = j + k;

            }
        }
    }

    public void sort(ArrayList<Integer> arr, int r) {
        quicksort(arr, 0, arr.size() - 1, r);
    }

    private void quicksort(ArrayList<Integer> arr, int left, int right, int r) {
        if (left < right) {
            int pivot = partition(arr, left, right, r);
            quicksort(arr, left, pivot - 1, r);
            quicksort(arr, pivot + 1, right, r);
        }
    }

    private int partition(ArrayList<Integer> arr, int left, int right, int r) {
        int pivot = arr.get(right);
        int i = left - 1;
        for (int j = left; j < right; j++) {
            if (c[r][arr.get(j)] < c[r][pivot]) {
                i++;
                int temp = arr.get(i);
                arr.set(i, arr.get(j));
                arr.set(j, temp);
            }
        }
        int temp = arr.get(i + 1);
        arr.set(i + 1, arr.get(right));
        arr.set(right, temp);
        return i + 1;
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
            sort(cList, r);

            int range = (int)(((double)cList.size()) * alpha)+1;

            int c = cList.get(rnd[rnd_index++]);

            s.add(c);
            cList.remove(Integer.valueOf(c));
            r = c;
        }

        s.add(0);
        return s;
    }

    private void swap(ArrayList<Integer> s, int i, int j){
        Collections.swap(s, i, j);
    }

    private void reverse(ArrayList<Integer> s, int i, int j){
        Collections.reverse(s.subList(i,j+1));
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
    }

    private void search_swap(ArrayList<Integer> s, double [][][] seq){
        double cost_best = Double.MAX_VALUE;
        double cost_new;
        double cost_concat_1;
        double cost_concat_2;
        double cost_concat_3;
        double cost_concat_4;
        int I = -1;
        int J = -1;
        
        for(int i = 1; i < dimension-1; i++){
            int i_prev = i - 1;
            int i_next = i + 1;
            int s_i_prev = s.get(i_prev);
            int s_i_next = s.get(i_next);
            int s_i = s.get(i);

            // immediate nodes case

            cost_concat_1 =                 seq[0][i_prev][T] + c[s_i_prev][s_i_next];
            cost_concat_2 = cost_concat_1 + seq[i][i_next][T] + c[s_i][s.get(i_next+1)];

            cost_new = seq[0][i_prev][C]                                                            // 1st subseq 
                    + seq[i][i_next][W]           * (cost_concat_1) + c[s_i_next][s_i]              // concat 2nd subseq
                    + seq[i_next+1][dimension][W] * (cost_concat_2) + seq[i_next+1][dimension][C];  // concat 3rd subseq

            if(cost_new < cost_best){
                cost_best = cost_new - EPSILON;
                I = i;
                J = i_next;
            }

            for(int j = i_next+1; j < dimension; j++){
                int j_next = j+1;
                int j_prev = j-1;
                int s_j = s.get(j);
                int s_j_prev = s.get(j_prev);
                int s_j_next = s.get(j_next);

                cost_concat_1 = seq[0][i_prev][T] + c[s_i_prev][s_j];
                cost_concat_2 = cost_concat_1 + c[s_j][s_i_next];
                cost_concat_3 = cost_concat_2 + seq[i_next][j_prev][T] + c[s_j_prev][s_i];
                cost_concat_4 = cost_concat_3 + c[s_i][s_j_next];

                cost_new = seq[0][i_prev][C]                                                        /* first subseq */
                        + cost_concat_1                                                             /* concatenate second subseq (single node) */
                        + seq[i_next][j_prev][W] * cost_concat_2 + seq[i_next][j_prev][C]           /* concatenate third subseq */
                        + cost_concat_3                                                             /* concatenate fourth subseq (single node) */
                        + seq[j_next][dimension][W] * cost_concat_4 + seq[j_next][dimension][C];    /* concatenate fifth subseq */

                if(cost_new < cost_best){
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                }

            }

        }

        if(cost_best < seq[0][dimension][C] - EPSILON){
            swap(s, I, J);
            subseq_load(s, seq);
            improv_flag = true;
        }
    } 

    private void search_two_opt(ArrayList<Integer> s, double [][][] seq){
        int I = -1;
        int J = -1;
        double cost_best = Double.MAX_VALUE;
        double cost_new;
        double cost_concat_1;
        double cost_concat_2;

        for(int i = 1; i < dimension-1; i++){
            int i_prev = i -1;
            int s_i = s.get(i);
            int s_i_prev = s.get(i_prev);
            double rev_seq_cost = seq[i][i+1][T];

            for(int j = i+2; j < dimension; j++){
                int j_next = j+1;
                int s_j_next = s.get(j_next);
                int s_j_prev = s.get(j-1);
                int s_j = s.get(j);

                rev_seq_cost += c[s_j_prev][s_j] * (seq[i][j][W]-1);

                cost_concat_1 =                 seq[0][i_prev][T] + c[s_j][s_i_prev];
                cost_concat_2 = cost_concat_1 + seq[i][j][T] + c[s_j_next][s_i];

                cost_new = seq[0][i_prev][C]                                                        /*        1st subseq */
                        + seq[i][j][W]              * cost_concat_1 + rev_seq_cost                  /* concat 2nd subseq (reversed seq) */
                        + seq[j_next][dimension][W] * cost_concat_2 + seq[j_next][dimension][C];    /* concat 3rd subseq */

                if(cost_new < cost_best){
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                }
            }
        }

        if(cost_best < seq[0][dimension][C] - EPSILON){
            reverse(s, I, J);
            subseq_load(s, seq);
            improv_flag = true;
        }
    } 

    private void search_reinsertion(ArrayList<Integer> s, double [][][] seq, int opt){
        double cost_best = Double.MAX_VALUE;
        double cost_new;
        double cost_concat_1;
        double cost_concat_2;
        double cost_concat_3;
        int I = -1;
        int J = -1;
        int POS = -1;

        for(int i = 1; i < dimension - opt + 1; i++){
            int j = opt + i - 1;
            int i_prev = i-1;
            int j_next = j+1;
            int s_i = s.get(i);
            int s_j = s.get(j);
            int s_i_prev = s.get(i_prev);
            int s_j_next = s.get(j_next);

            // k -> reinsertion places
            for(int k = 0; k < i_prev; k++){
                int k_next = k+1;
                int s_k = s.get(k);
                int s_k_next = s.get(k_next);
                
                cost_concat_1 = seq[0][k][T] + c[s_k][s_i];
                cost_concat_2 = cost_concat_1 + seq[i][j][T] + c[s_j][s_k_next];
                cost_concat_3 = cost_concat_2 + seq[k_next][i_prev][T] + c[s_i_prev][s_j_next];

                cost_new = seq[0][k][C]                                                             /*        1st subseq */
                        + seq[i][j][W]              * cost_concat_1 + seq[i][j][C]                  /* concat 2nd subseq (reinserted seq) */
                        + seq[k_next][i_prev][W]    * cost_concat_2 + seq[k_next][i_prev][C]        /* concat 3rd subseq */
                        + seq[j_next][dimension][W] * cost_concat_3 + seq[j_next][dimension][C];    /* concat 4th subseq */

                if(cost_new < cost_best){
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                    POS = k;
                }
            }

            for(int k = i+opt; k < dimension; k++){
                int k_next = k+1;
                int s_k = s.get(k);
                int s_k_next = s.get(k_next);

                cost_concat_1 = seq[0][i_prev][T] + c[s_i_prev][s_j_next];
                cost_concat_2 = cost_concat_1 + seq[j_next][k][T] + c[s_k][s_i];
                cost_concat_3 = cost_concat_2 + seq[i][j][T] + c[s_j][s_k_next];

                cost_new = seq[0][i_prev][C]                                                        /*      1st subseq */
                        + seq[j_next][k][W]         * cost_concat_1 + seq[j_next][k][C]             /* concat 2nd subseq */
                        + seq[i][j][W]              * cost_concat_2 + seq[i][j][C]                  /* concat 3rd subseq (reinserted seq) */
                        + seq[k_next][dimension][W] * cost_concat_3 + seq[k_next][dimension][C];    /* concat 4th subseq */

                if(cost_new < cost_best){
                    cost_best = cost_new - EPSILON;
                    I = i;
                    J = j;
                    POS = k;
                }
            }
        }

        if(cost_best < seq[0][dimension][C] - EPSILON){
            reinsert(s, I, J, POS+1);
            subseq_load(s, seq);
            improv_flag = true;
        }
    } 

    private void RVND(ArrayList<Integer> s, double [][][] subseq){
        ArrayList<Integer> neighbd_list = new ArrayList<Integer>() {{
            add(SWAP);
            add(TWO_OPT);
            add(REINSERTION);
            add(OR_OPT2);
            add(OR_OPT3);
        }};

        while(!neighbd_list.isEmpty()){
            
            int i_rand = rnd[rnd_index++];
            int neighbd = neighbd_list.get(i_rand);

            improv_flag = false;

            switch(neighbd){
                case REINSERTION:
                    search_reinsertion(s, subseq, REINSERTION);
                    break;
                case OR_OPT2:
                    search_reinsertion(s, subseq, OR_OPT2);
                    break;
                case OR_OPT3:
                    search_reinsertion(s, subseq, OR_OPT3);
                    break;
                case SWAP:
                    search_swap(s, subseq);
                    break;
                case TWO_OPT:
                    search_two_opt(s, subseq);
                    break;
            }

            if(improv_flag){
                neighbd_list.clear();
                neighbd_list = new ArrayList<Integer>() {{
                    add(SWAP);
                    add(TWO_OPT);
                    add(REINSERTION);
                    add(OR_OPT2);
                    add(OR_OPT3);
                }};
            }else
                neighbd_list.remove(i_rand);
        }
    }

    private ArrayList<Integer> perturb(ArrayList<Integer> sl){
        ArrayList<Integer> sl_cpy = new ArrayList<>(sl); 

        int A_start = 1, A_end = 1;
        int B_start = 1, B_end = 1;

        while((A_start <= B_start && B_start <= A_end) || (B_start <= A_start && A_start <= B_end)){

            A_start = rnd[rnd_index++];
            A_end = A_start + rnd[rnd_index++];

            B_start = rnd[rnd_index++];
            B_end = B_start + rnd[rnd_index++];
        }

        if(A_start < B_start){
            reinsert(sl_cpy, B_start, B_end-1, A_end);
            reinsert(sl_cpy, A_start, A_end-1, B_end);
        }else{
            reinsert(sl_cpy, A_start, A_end-1, B_end);
            reinsert(sl_cpy, B_start, B_end-1, A_end);
        }

        return sl_cpy;
    }

    public void solve(){
        double cost_best = Double.MAX_VALUE;
        ArrayList<Integer> s_best = new ArrayList<>();
        for(int i = 0; i < Imax; i++){
            double alpha = R[rnd[rnd_index++]];

            System.out.print("[+] Local Search ");
            System.out.println(i+1);
            System.out.println("\t[+] Constructing Inital Solution..");
            ArrayList<Integer> s = construction(alpha);
            ArrayList<Integer> sl = new ArrayList<>(s); 
            subseq_load(s, subseq); 
            double rvnd_cost_best = subseq[0][dimension][C] - EPSILON;
            double rvnd_cost_crnt;

            int iterILS = 0;
            System.out.println("\t[+] Looking for the best Neighbor..");
            while(iterILS < Iils){
                RVND(s, subseq);
                rvnd_cost_crnt = subseq[0][dimension][C] - EPSILON; 
                if(rvnd_cost_crnt < rvnd_cost_best){
                    rvnd_cost_best = rvnd_cost_crnt;
                    sl.clear();
                    sl = new ArrayList<>(s);
                    iterILS = 0;
                }
                s = perturb(sl);
                subseq_load(s, subseq);
                iterILS++;
            }
            subseq_load(sl, subseq);
            double sl_cost = subseq[0][dimension][C] - EPSILON;
            
            if(sl_cost < cost_best){
                cost_best = sl_cost;
                s_best = new ArrayList<>(sl);
            }

            System.out.print("\tCurrent best solution cost: ");
            System.out.println(cost_best);
        }

        System.out.print("COST: ");
        System.out.println(cost_best);
        System.out.print("SOLUTION: ");
        System.out.println(s_best);
        System.out.print("Construction: ");
        System.out.println(t_construction/10e8);
        System.out.print("swap: ");
        System.out.println(t_swap/10e8);
        System.out.print("2_opt: ");
        System.out.println(t_two_opt/10e8);
        System.out.print("or_opt2: ");
        System.out.println(t_or_opt2/10e8);
        System.out.print("or_opt3: ");
        System.out.println(t_or_opt3/10e8);
        System.out.print("reinsertion: ");
        System.out.println(t_reinsertion/10e8);
        System.out.print("subseq_load: ");
        System.out.println(t_subseq/10e8);
    }
}
