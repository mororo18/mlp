import  java.lang.System;


class tArray {

    private int [] arr;
    private int [] arrBuffer;
    private int sizeCrnt;
    private int sizeBuffer;

    private int size;

    interface Comparator {
        public int evaluate(int i, int j);
    }

    tArray() {
        this.size = 0;
        this.sizeCrnt = 0;
    }

    tArray(int []arr) {
        this.arr        = arr.clone();
        this.sizeCrnt   = arr.length;
        this.size       = arr.length;

        this.arrBuffer  = new int[3];
        this.sizeBuffer = 0;
    }

    tArray(int size) {
        this.arr        = new int[size];
        this.arrBuffer  = new int[3];
        this.sizeCrnt   = size;
        this.sizeBuffer = 0;
        this.size       = arr.length;
    }

    public int size() {return sizeCrnt;}
    public boolean isEmpty() {
        if (sizeCrnt > 0) return false; 
        else              return true;
    }

    public void add(int element) {
        if (sizeCrnt >= size) {
        } else {
            this.arr[sizeCrnt] = element;
            sizeCrnt++;
        }
    }
    public void reserve(int sz) {
        this.arr = new int[sz];
        this.arrBuffer = new int[3];
        this.size = sz;
        this.sizeCrnt = 0;
        this.sizeBuffer = 0;
    }

    public tArray getArrayCopy() {
        return new tArray(this.arr.clone());
        //return arr.clone();
    }

    public void assign(int []arr) {
        this.arr = arr.clone();
        sizeCrnt = arr.length;
    }

    private void shift(int from, int to, int sz) {
        if (from < to) {
            int dist = to - from;

            for (int i = from+sz-1, j = to+sz-1; i >= from; i--,j--) {
                arr[j] = arr[i];
            }
        } else {
            for (int i = from, j = to; i < from+sz; i++, j++) {
                arr[j] = arr[i];
            }
        }

    }

    public void swap(int i, int j) {
        int tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    }

    public void reverse(int i, int j) {
        int m = (i+j)/2;
        for (int first = i, last = j;
                first <= m;
                first++, last--) {

            swap(first, last);
        }
    }

    public void reinsert(int i, int j, int pos) {
        this.sizeBuffer = j-i+1;

        System.arraycopy(this.arr, i, this.arrBuffer, 0, this.sizeBuffer);

        if (pos < i) {
            int sz = i-pos;
            shift(pos, j+1-sz, sz);
            System.arraycopy(this.arrBuffer, 0, arr, pos, this.sizeBuffer);
        } else {
            int sz = pos-j-1;
            shift(j+1, i, sz);
            System.arraycopy(this.arrBuffer, 0, arr, i+sz, this.sizeBuffer);
        }
        
        if (feasibility() != true) {
            System.out.println("infisibu");
            System.exit(0);
        }
    }

    public void remove(int index) {
        shift(index+1, index, sizeCrnt-index-1);
        sizeCrnt--;
    }

    public void fill() {
        this.sizeCrnt = this.size;
    }
    
    public int set(int i, int element) {
        return arr[i] = element;
    }

    public int get(int i) {
        return arr[i];
    }

    public boolean feasibility() {
        boolean [] check = new boolean[this.sizeCrnt];
        for (int i = 0; i < sizeCrnt; i++) {
            check[this.arr[i]] = true;
        }

        for (int i = 0; i < sizeCrnt; i++) {
            if (!check[i]) return false;
        }
         return true;
    }


    public void sort(Comparator comp) {

        for (int i = 0; i < sizeCrnt; i++) {
            for (int j = 0; j < sizeCrnt-i-1; j++) {
                if (comp.evaluate(arr[j], arr[j+1]) > 0) {
                    int tmp = arr[j];
                    arr[j] = arr[j+1];
                    arr[j+1] = tmp;
                }
            }
        }
    }

    public void printPretty() {
        System.out.print("[");
        for (int i = 0; i < sizeCrnt-1; i++) {
            System.out.print(arr[i]);
            System.out.print(", ");
        }
        System.out.print(arr[sizeCrnt-1]);
        System.out.print("]");
        System.out.println();
    }


    public void print() {
        for (int i = 0; i < sizeCrnt; i++) {
            System.out.print(arr[i]);
            System.out.print(" ");
        }
        System.out.println();
    }

}
