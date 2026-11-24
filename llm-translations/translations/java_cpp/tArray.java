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

        this.arrBuffer  = new int[this.size];
        this.sizeBuffer = 0;
    }

    tArray(int size) {
        this.arr        = new int[size];
        this.arrBuffer  = new int[size];
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
        this.arr        = new int[sz];
        this.arrBuffer  = new int[sz];
        this.size       = sz;
        this.sizeCrnt   = 0;
        this.sizeBuffer = 0;
    }

    public tArray getArrayCopy() {
        return new tArray(this.arr.clone());
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
        boolean [] check = new boolean[this.sizeCrnt-1];

        for (int i = 0; i < check.length; i++) {
            check[i] = false;
        }

        for (int i = 0; i < this.arr.length; i++) {
            check[this.arr[i]] = true;
        }

        for (int i = 0; i < check.length; i++) {
            if (!check[i]) return false;
        }
         return true;
    }


    public void sort(int r, tInfo info) {
	quicksort(arr, 0, size() - 1, info, r);
    }

    private static void quicksort(int[] arr, int left, int right, tInfo info, int r) {
        if (left < right) {
            int pivot = partition(arr, left, right, info, r);
            quicksort(arr, left, pivot - 1, info, r);
            quicksort(arr, pivot + 1, right, info, r);
        }
    }

    private static int partition(int[] arr, int left, int right, tInfo info, int r) {
        int pivot = arr[right];
        int i = left - 1;
        for (int j = left; j < right; j++) {
            if (info.getCost(r, arr[j]) < info.getCost(r, pivot)) {
                i++;
                int temp = arr[i];
                arr[i] = arr[j];
                arr[j] = temp;
            }
        }
        int temp = arr[i + 1];
        arr[i + 1] = arr[right];
        arr[right] = temp;
        return i + 1;
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
