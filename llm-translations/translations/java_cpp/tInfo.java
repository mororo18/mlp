class tInfo {
    private int     []  rnd;
    private int     rnd_index;

    static final double EPSILON = 1e-16;
    private int     dimension;
    private double  [][] c;

    static final int T = 0;
    static final int C = 1;
    static final int W = 2;

    public static final int SWAP        = 0;
    public static final int REINSERTION = 1;
    public static final int OR_OPT2     = 2;
    public static final int OR_OPT3     = 3;
    public static final int TWO_OPT     = 4;

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
