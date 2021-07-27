class Main {
    public static void main(String args[]){
        read();
    }

    private static void read(){
        GILS_RVND tsp = new GILS_RVND();
        long time = - System.nanoTime();
        tsp.solve();
        time += System.nanoTime();

        System.out.print("TIME: ");
        System.out.println(time/10e8);
    }
}
