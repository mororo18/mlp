import  java.lang.System;


class tArray {

    private int [] arr;
    private int size;
    private int sizeCrnt;

    interface Comparator {
        public int evaluate(int i, int j);
    }

    tArray() {
        this.size = 0;
        this.sizeCrnt = 0;
    }

    tArray(int []arr) {
        this.arr = arr.clone();
        sizeCrnt = arr.length;
        size = arr.length;
    }

    tArray(int size) {
        arr = new int[size];
        sizeCrnt = size;
        size = arr.length;
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
        this.size = sz;
        this.sizeCrnt = 0;
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
        int size = j-i+1;
        int [] seq = new int[size];

        System.arraycopy(this.arr, i, seq, 0, size);

        if (pos < i) {
            int sz = i-pos;
            shift(pos, j+1-sz, sz);
            System.arraycopy(seq, 0, arr, pos, size);
        } else {
            int sz = pos-j-1;
            shift(j+1, i, sz);
            System.arraycopy(seq, 0, arr, i+sz, size);
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
