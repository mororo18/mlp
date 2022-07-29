import sun.misc.Unsafe;
import java.lang.reflect.Field;

class Test {


    public static Unsafe loadUnsafe() throws SecurityException,
           NoSuchFieldException, IllegalArgumentException,
           IllegalAccessException
   {
       Field theUnsafe = Unsafe.class.getDeclaredField("theUnsafe");
       theUnsafe.setAccessible(true);
       Unsafe unsafe = (Unsafe) theUnsafe.get(null);
       return unsafe;
   }

   //public static Unsafe getUnsafe() throws NoSuchFieldException 
   //{return unsafe;}

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
}
