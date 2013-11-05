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


// Arguments separated to reduce recursion overheads.
template<typename T>
struct transpose_arguments {
  T* row_major;
  T* col_major;
  int rows;
  int cols;
};

// This helper function is separated from the cilk dc_transpose function for compiler's benefit.
template<typename T>
void serial_transpose (T* row_major, T* col_major, int rows, int cols, int r_offset, int c_offset, int r_size, int c_size) {
   for (int r = r_offset; r < r_offset+r_size; r++) {
    int offset2 = r*cols;
   for (int c = c_offset; c < c_offset+c_size; c++) {
        col_major[rows*c + r] = row_major[offset2 + c];
      }
    }
}

template<typename T>
void dc_transpose(transpose_arguments<T>* arguments, int r_offset, int c_offset, int r_size, int c_size) {
  if (r_size*c_size < 512) {
    return serial_transpose(arguments->row_major,arguments->col_major,arguments->rows, arguments->cols, r_offset, c_offset, r_size, c_size); // separate function to use better compiler.
  }
  if (r_size > c_size) {
    // divide on rows.
    int new_r_size1 = r_size / 2;
    int new_r_size2 = r_size - new_r_size1;
    cilk_spawn dc_transpose(arguments, r_offset, c_offset, new_r_size1, c_size);
    dc_transpose(arguments, r_offset + new_r_size1, c_offset, new_r_size2, c_size);
    cilk_sync;
  } else {
    int new_c_size1 = c_size/2;
    int new_c_size2 = c_size - new_c_size1;
    cilk_spawn dc_transpose(arguments, r_offset, c_offset, r_size, new_c_size1);
    dc_transpose(arguments, r_offset, c_offset + new_c_size1, r_size, new_c_size2);
    cilk_sync;
  }
}

// Transposes a matrix in row-major form to column-major. Currently creates
// a new array. Should be able to replace this with a cache-efficient
// implementation to get better performance.
template<typename T>
T* transpose (T* row_major, int rows, int cols) {
    T* column_major = (T*) malloc(sizeof(T) * rows * cols);
    transpose_arguments<T>* arguments = (transpose_arguments<T>*) calloc(1,sizeof(transpose_arguments<T>));
    arguments->row_major = row_major;
    arguments->col_major = column_major;
    arguments->rows = rows;
    arguments->cols = cols;
    dc_transpose<T>(arguments, 0, 0, rows, cols);
   // Below is the old transpose code.
   /*
    cilk_for (int col = 0; col < cols; col++) {
        for (int j = 0; j < rows; j++) {
          column_major[col*rows + j] = row_major[j*cols + col];
        }
    }*/
    return column_major;
}


#endif
