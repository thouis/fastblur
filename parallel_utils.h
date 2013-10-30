#ifndef PARALLEL_UTILS_H
#define PARALLEL_UTILS_H

#include <cilk/cilk.h>

// Reduces two arrays, puts the result in the
//  left, and destroys the right.
template<typename T>
T* reduce_arrays_sub(T* left, T* right, int size) {
  cilk_for (int i = 0; i < size; i++) {
    left[i] -= right[i];
  }
  free(right);
  return left;
}

// Transposes a matrix in row-major form to column-major. Currently creates
// a new array. Should be able to replace this with a cache-efficient
// implementation to get better performance.
template<typename T>
T* transpose (T* row_major, int rows, int cols) {
    T* column_major = (T*) malloc(sizeof(T) * rows * cols);

    cilk_for (int col = 0; col < cols; col++) {
        for (int j = 0; j < rows; j++) {
          column_major[col*rows + j] = row_major[j*cols + col];
        }
    }
    return column_major;
}
#endif
