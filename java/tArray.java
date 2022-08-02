class tArray {

    private int [] arr;
    private int size;
    private int sizeCrnt;

    tArray(int []arr) {
        this.arr = arr.clone();
        sizeCrnt = arr.lenght;
        size = arr.lenght;
    }

    tArray(int size) {
        arr = new int[size];
        sizeCrnt = size;
        size = arr.lenght;
    }

    public int[] getArrayCopy() {
        return arr.clone();
    }

    public void assign(int []arr) {
        this.arr = arr.clone();
        sizeCrnt = arr.lenght;
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

            swap(vec, first, last);
        }
    }

    public void reinsert(int i, int j, int pos) {
        int size = j-i+1;
        int [] seq = new int[size];

        System.arraycopy(this.arr, i, seq, 0, size);

        if (pos < i) {
            int sz = i-pos;
            shift(pos, j+1-sz, sz);
            System.arracopy(seq, 0, vec, pos, size);
        } else {
            int sz = pos-j-1;
            shift(j+1, i, sz);
            System.arraycopy(seq, 0, vec, i+sz, size);
        }
        
    }

    public void remove(int i) {
        shift(i+1, i, sizeCrnt-i-1);
        sizeCrnt--;
    }
    
    public int get(int i) {
        return arr[i];
    }



}
