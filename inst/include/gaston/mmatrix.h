#ifndef MMATRIX_H
#define MMATRIX_H

#include <fstream>
#include <iostream>
#include <string>
#include <system_error>
#include <vector>

#include "mio.hpp"

/*
 The goal of this class is to add the possibility for big matrices
 to become a memory mapped file using mio library */
template <typename T>
class MMatrix
{
public:
    // Constructor opening the file containing the matrix
    // if path exists else creating one.
    MMatrix(std::string path, size_t ncol, size_t nrow);
    // Destructor flushing changes to disk before unmapping
    ~MMatrix();

    size_t nrow() const;
    size_t ncol() const;
    std::string path() const;
    T *data() const;

    // IN BASE 0 :
    T &operator[](size_t ind);
    const T &operator[](size_t ind) const;
    T &operator()(size_t i, size_t j);
    const T &operator()(size_t i, size_t j) const;
    T &at(size_t ind) const;
    T &at(size_t i, size_t j) const;
    template <typename U>
    std::vector<U> sum() const;

    //  helper function, creates and loads a descriptor path
    void create_descriptor_file();

protected:
    // Number of columns of the matrix, (base 1).
    size_t ncol_;
    // Number of rows of the matrix (base 1).
    size_t nrow_;
    size_t size_;
    // (Relative ?) path of the file containing the matrix
    std::string path_;
    // Mio object handling the matrix.
    mio::mmap_sink matrix_file_;
    // type T pointer to the first byte of data in matrix_file_
    T *data_ptr_;
};

#endif // MMATRIX_H