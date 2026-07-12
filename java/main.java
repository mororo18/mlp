import java.lang.management.ManagementFactory;
import com.sun.management.OperatingSystemMXBean;

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
        OperatingSystemMXBean osBean = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean();
        long cpuTime = - osBean.getProcessCpuTime();
        long time = - System.nanoTime();
        tsp.solve(verbose);
        time += System.nanoTime();
        cpuTime += osBean.getProcessCpuTime();

        System.out.print("TIME: ");
        System.out.println(cpuTime/10e8);
        System.out.print("wall clock (s): ");
        System.out.println(time/10e8);
    }
}
