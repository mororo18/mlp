import sun.misc.Unsafe;
import java.lang.reflect.Field;

class MatrixTest {

    private         Unsafe unsafe;
    private static final  int INT_SIZE = 4;
    private         int rowSize;
    private         int colSize;
    private         int size;
    private         int rowOffSet;
    private         long addrList;
    private final long ADDR_SIZE = 8;

    private 
    void loadUnsafe() throws SecurityException,
         NoSuchFieldException, IllegalArgumentException,
         IllegalAccessException
         {
             Field theUnsafe = Unsafe.class.getDeclaredField("theUnsafe");
             theUnsafe.setAccessible(true);
             unsafe = (Unsafe) theUnsafe.get(null);
         }

    MatrixTest(int row, int col) {
        this.rowSize = row;                                 // quantidade de elementos por linha
        this.colSize = col;                                 // quantidade de elementos por coluna
        this.size = row*col;                                // tamanho da matriz em bytes
        this.rowOffSet = this.colSize * this.INT_SIZE;      // tamanho da linha me bytes
        try {
        loadUnsafe();
        } catch (Exception e) {
            System.out.println("qebro");
            System.exit(0);
        }

        this.addrList = unsafe.allocateMemory(this.ADDR_SIZE * this.rowSize);

        for (int r = 0; r < this.rowSize; r++) {
            long addr = unsafe.allocateMemory(INT_SIZE * this.colSize);
            unsafe.putLong(this.addrList + ADDR_SIZE * r, addr);
        }
    }


    private long offset(int i, int j) {
        long colOffSet = j * this.INT_SIZE;
        long totalRowOffSet = i * this.rowOffSet;  
        return totalRowOffSet + colOffSet;
    } 

    private long getRowAddr(int row) {
        return unsafe.getLong(this.addrList + this.ADDR_SIZE * row);
    }

    public void free() {
        for (int r = 0; r < this.rowSize; r++) {
            unsafe.freeMemory(getRowAddr(r));
        }
        unsafe.freeMemory(this.addrList);
    }

    public int get(int i, int j) {
      //long colOffSet = j * this.INT_SIZE;
      //long totalRowOffSet = i * this.rowOffSet;  
        return unsafe.getInt(getRowAddr(i) + j * this.INT_SIZE);
    }

    public void set(int i, int j, int value) {
      //long colOffSet = j * this.INT_SIZE;
      //long totalRowOffSet = i * this.rowOffSet;  
        unsafe.putInt(getRowAddr(i) + j * this.INT_SIZE, value);
        //unsafe.putInt(this.addr + j * this.INT_SIZE + i * this.rowOffSet, value);
        //unsafe.putInt(this.addr + totalRowOffSet + colOffSet, value);
    }

   /*
    public static void main(String[] args) {
        System.out.println("opas");

        Unsafe unsafe;
        try {
            unsafe = loadUnsafe();
            byte INT_SIZE = 4;
            final long addr = unsafe.allocateMemory(INT_SIZE * 10);

            int a = 123;
            loadUnsafe().putInt(addr + INT_SIZE *5, (int)123);
            loadUnsafe().putInt(addr + INT_SIZE *6, (int)124);
            loadUnsafe().putInt(addr + INT_SIZE *4, (int)122);
            loadUnsafe().putInt(addr , (int)12);
            loadUnsafe().putInt(addr + INT_SIZE *3, (int)22);
            int meuInt = loadUnsafe().getInt(addr + INT_SIZE *6);

            System.out.println(meuInt);

        } catch (Exception e) {
            System.out.println("unsafe qebro");
            System.exit(0);
        }

    }
    */
}
