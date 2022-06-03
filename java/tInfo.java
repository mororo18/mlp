class tInfo {
    private int     []  rnd;
    private int     rnd_index;

    final double EPSILON = 1e-16;
    private int     dimension;
    private double  [][] c;

    static final int C = 0;
    static final int T = 1;
    static final int W = 2;

    public static  int SWAP       = 0;
    static  int REINSERTION= 1;
    static  int OR_OPT2    = 2;
    static  int OR_OPT3    = 3;
    static  int TWO_OPT    = 4;

    tInfo( int dimen, double [][] cost, int [] rnd_arr) {
        rnd = rnd_arr;
        dimension = dimen;
        c = cost;
        rnd_index = 0;
    }

    int getDimen()                  {return dimension;}
    double getCost(int i, int j)    {return c[i][j];}
    int rndCrnt()                   {return rnd[rnd_index++];}
}
