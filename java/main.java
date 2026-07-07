class Main {
    public static void main(String args[]){
        boolean verbose = false;
        for (String arg : args) {
            if (arg.equals("-v") || arg.equals("--verbose")) {
                verbose = true;
            }
        }
        read(verbose);
    }

    private static void read(boolean verbose){
        GILS_RVND tsp = new GILS_RVND();
        long time = - System.nanoTime();
        tsp.solve(verbose);
        time += System.nanoTime();

        System.out.print("TIME: ");
        System.out.println(time/10e8);
    }
}
