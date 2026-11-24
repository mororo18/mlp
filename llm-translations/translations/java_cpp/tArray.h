#ifndef TARRAY_H
#define TARRAY_H

#include "tInfo.h"

#include <algorithm>
#include <iostream>
#include <vector>

class tArray {
private:
    std::vector<int> arr;
    std::vector<int> arrBuffer;
    int sizeCrnt;
    int sizeBuffer;

    int size;

public:
    tArray() {
        this->size      = 0;
        this->sizeCrnt  = 0;
        this->sizeBuffer = 0;
    }

    tArray(const std::vector<int> & arr) {
        this->arr       = arr;
        this->sizeCrnt  = (int) arr.size();
        this->size      = (int) arr.size();

        this->arrBuffer.assign(this->size, 0);
        this->sizeBuffer = 0;
    }

    tArray(int size) {
        this->arr.assign(size, 0);
        this->arrBuffer.assign(size, 0);
        this->sizeCrnt  = size;
        this->sizeBuffer = 0;
        this->size      = (int) this->arr.size();
    }

    int size_() {return sizeCrnt;}
    bool isEmpty() {
        if (sizeCrnt > 0) return false;
        else              return true;
    }

    void add(int element) {
        if (sizeCrnt >= size) {
        } else {
            this->arr[sizeCrnt] = element;
            sizeCrnt++;
        }
    }

    void reserve(int sz) {
        this->arr.assign(sz, 0);
        this->arrBuffer.assign(sz, 0);
        this->size      = sz;
        this->sizeCrnt  = 0;
        this->sizeBuffer = 0;
    }

    tArray getArrayCopy() {
        return tArray(this->arr);
    }

    void assign(const std::vector<int> & arr) {
        this->arr = arr;
        sizeCrnt = (int) arr.size();
    }

private:
    void shift(int from, int to, int sz) {
        if (from < to) {
            for (int i = from+sz-1, j = to+sz-1; i >= from; i--, j--) {
                arr[j] = arr[i];
            }
        } else {
            for (int i = from, j = to; i < from+sz; i++, j++) {
                arr[j] = arr[i];
            }
        }
    }

public:
    void swap(int i, int j) {
        int tmp = arr[i]; arr[i] = arr[j]; arr[j] = tmp;
    }

    void reverse(int i, int j) {
        int m = (i+j)/2;
        for (int first = i, last = j;
                first <= m;
                first++, last--) {

            swap(first, last);
        }
    }

    void reinsert(int i, int j, int pos) {
        this->sizeBuffer = j-i+1;

        std::copy(this->arr.begin()+i, this->arr.begin()+i+this->sizeBuffer, this->arrBuffer.begin());

        if (pos < i) {
            int sz = i-pos;
            shift(pos, j+1-sz, sz);
            std::copy(this->arrBuffer.begin(), this->arrBuffer.begin()+this->sizeBuffer, arr.begin()+pos);
        } else {
            int sz = pos-j-1;
            shift(j+1, i, sz);
            std::copy(this->arrBuffer.begin(), this->arrBuffer.begin()+this->sizeBuffer, arr.begin()+i+sz);
        }
    }

    void remove(int index) {
        shift(index+1, index, sizeCrnt-index-1);
        sizeCrnt--;
    }

    void fill() {
        this->sizeCrnt = this->size;
    }

    int set(int i, int element) {
        return arr[i] = element;
    }

    int get(int i) {
        return arr[i];
    }

    bool feasibility() {
        std::vector<bool> check(this->sizeCrnt-1);

        for (int i = 0; i < (int) check.size(); i++) {
            check[i] = false;
        }

        for (int i = 0; i < (int) this->arr.size(); i++) {
            check[this->arr[i]] = true;
        }

        for (int i = 0; i < (int) check.size(); i++) {
            if (!check[i]) return false;
        }
        return true;
    }

    void sort(int r, tInfo & info) {
        quicksort(arr, 0, size_() - 1, info, r);
    }

private:
    static void quicksort(std::vector<int> & arr, int left, int right, tInfo & info, int r) {
        if (left < right) {
            int pivot = partition(arr, left, right, info, r);
            quicksort(arr, left, pivot - 1, info, r);
            quicksort(arr, pivot + 1, right, info, r);
        }
    }

    static int partition(std::vector<int> & arr, int left, int right, tInfo & info, int r) {
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

public:
    void printPretty() {
        std::cout << "[";
        for (int i = 0; i < sizeCrnt-1; i++) {
            std::cout << arr[i];
            std::cout << ", ";
        }
        std::cout << arr[sizeCrnt-1];
        std::cout << "]";
        std::cout << std::endl;
    }

    void print() {
        for (int i = 0; i < sizeCrnt; i++) {
            std::cout << arr[i];
            std::cout << " ";
        }
        std::cout << std::endl;
    }
};

#endif
