#include <cassert>
#include <filesystem>
#include <fstream> // for ofstream, loading the file
#include <iostream> // for std::cout
#include <stdexcept>
#include <system_error> // for std::error_code

#include "mio.hpp"
#include "mmatrix.h"

// Constructor opening the file containing the matrix if path exists, else
// creating one.
template <typename T>
MMatrix<T>::MMatrix(std::string path, size_t ncol, size_t nrow)
    : ncol_(ncol)
    , nrow_(nrow)
    , path_(path)
{
    size_ = ncol * nrow;
    size_t matrix_size = size_ * sizeof(T);

    std::ifstream file(path, std::ifstream::binary);
    if (!file.good())
    {
        std::cout << "This file does not exist, creating one ...  ";
        std::ofstream file(path, std::ios::binary); // to load with \0
        if (!file || !file.is_open())
        {
            throw std::runtime_error("Failed to open the file;");
        }
        /* Writing a null byte to the end of the file
        after using seekp to the before last byte will
        ensure we have a non-empty file of the desired size
        without loading it into memory */
        file.seekp(matrix_size - 1);
        file.put('\0');
        // mio will reopen it
        file.close();
        std::cout << "Done !" << std::endl;
    }
    else {
        std::cout << "Using and overwritting already existing file: " << path <<"...\n";
    }
    // then opening the file
    std::error_code error;

    matrix_file_ = mio::make_mmap_sink(path, 0, mio::map_entire_file, error);
    if (error)
    {
        std::string errMsg = "Error code " + std::to_string(error.value())
            + ", Failed to map the file : " + error.message();
        throw std::runtime_error(errMsg);
    }

    // TODO : see if necessary
    if (matrix_size > matrix_file_.size())
    {
        std::string errMsg = "Error : Trying to access a matrix of "
            + std::to_string(matrix_size) + " byte(s), so bigger than " + path
            + " file.";
        throw std::invalid_argument(errMsg);
    }
    data_ptr_ = reinterpret_cast<T *>(matrix_file_.data());
}

// Destructor flushing changes to disk before unmapping
template <typename T>
MMatrix<T>::~MMatrix()
{
    //TO DEBUG
    std::cout << "unmapping mmatrix " << path_ << "...\n";
    std::error_code error;
    if (matrix_file_.is_mapped())
    {
        matrix_file_.sync(error);
        if (error)
        {
            // here no exception not to disturb the unstacking
            std::cerr << "Failed to unsync the file: " << error.message()
                      << '\n';
        }
        matrix_file_.unmap();
    }
}

// Getters :
template <typename T>
size_t MMatrix<T>::nrow() const
{
    return nrow_;
}
template <typename T>
size_t MMatrix<T>::ncol() const
{
    return ncol_;
}
template <typename T>
std::string MMatrix<T>::path() const
{
    std::cout << "THIS IS THE " << path_ << "\n";
    return path_;
}
template <typename T>
T *MMatrix<T>::data() const
{
    return data_ptr_;
}

// Operator [] gives back the data at index, UNSAFE.
template <typename T>
T &MMatrix<T>::operator[](size_t ind)
{
    return data_ptr_[ind];
}
template <typename T>
const T &MMatrix<T>::operator[](size_t ind) const
{
    return data_ptr_[ind];
}

// Operator () gives back data at row i, col j
template <typename T>
T &MMatrix<T>::operator()(size_t i, size_t j)
{
    return data_ptr_[(j * nrow_) + i];
}
template <typename T>
const T &MMatrix<T>::operator()(size_t i, size_t j) const
{
    return data_ptr_[(j * nrow_) + i];
}

// Same as operators but safe, with bound checking.
// Can be used with one (like []) or two parameters.
template <typename T>
T &MMatrix<T>::at(size_t ind) const
{
    if (ind >= size_)
    {
        throw std::out_of_range("Index out of range");
    }
    return data_ptr_[ind];
}
template <typename T>
T &MMatrix<T>::at(size_t i, size_t j) const
{
    if (i > nrow_ || j > ncol_)
    {
        std::string errMsg = "Matrix indices out of range, only goes up to "
            + std::to_string(ncol_) + " columns and " + std::to_string(nrow_)
            + " rows";
        throw std::out_of_range(errMsg);
    }
    return data_ptr_[(j * nrow_) + i];
}

// UNSAFE, calling ()
template <typename T>
template <typename U>
std::vector<U> MMatrix<T>::sum() const
{
    std::vector<U> results(ncol_); // Allocates AND initialises w/ zero

    for (size_t i = 0; i < ncol_; ++i)
    {
        for (size_t j = 0; j < nrow_; ++j)
        {
            // Add the element in column i, row j.
            results[i] += static_cast<U>((*this)(j, i));
        }
    }
    return results;
}

// HELPER FUNCTIONS :

// get_type_name() is to get the template type by comparing it to known types
// it's necessary to complete the descriptor file
template <typename T>
inline std::string get_type_name()
{
    if (std::is_same<T, int>::value)
    {
        return "integer"; // written in full cos need for descfile
    }
    else if (std::is_same<T, float>::value)
    {
        return "float";
    }
    else if (std::is_same<T, double>::value)
    {
        return "double";
    }
    else if (std::is_same<T, char>::value)
    {
        return "char";
    }
    else
    {
        return "unknown"; // to expand later ?
    }
}

#if defined(_WIN32) || defined(_WIN64)
#    include <windows.h>
#else // For Posix and Posix compliant (e.g. macOS)
#    include <limits.h>
#    include <unistd.h>
#endif

inline std::string get_path()
{
    char buffer[PATH_MAX]; // PATH_MAX defined in limits.h

#if defined(_WIN32) || defined(_WIN64)
    // For Windows :
    if (GetCurrentDirectoryA(sizeof(buffer), buffer))
    {
        return std::string(buffer);
    }
    else
    {
        std::string errMsg =
            "Failed to retrieve the working directory on Windows";
        throw std::runtime_error(errMsg);
        return "";
    }
#else
    // For POSIX :
    if (getcwd(buffer, sizeof(buffer)))
    {
        return std::string(buffer);
    }
    else
    {
        std::string errMsg =
            "Failed to retrieve the working directory on POSIX";
        throw std::runtime_error(errMsg);
        return "";
    }
#endif
}

/* Method to create a descriptor file with metadata
making the matrix compatible with bigmemory for example
*/

// Create the descriptor file
// for now, will only deal with ncols nrows size and paths.
// will be fed to attach.big.matrix in bigmemory
// like so : y <- attach.big.matrix(matrix.desc)
// TODO : for a better readability, check if possible to add \n
template <typename T>
void MMatrix<T>::create_descriptor_file()
{
    // Generate the descriptor file name
    std::string descriptor_file_ = path_ + ".desc";

    // TODO : print to debug
    // std::cout << "path of descriptor_file " << descriptor_file_ << " \n";

    // Create and write to the descriptor file
    std::ofstream desc_file(descriptor_file_);
    if (!desc_file || !desc_file.is_open())
    {
        throw std::runtime_error("Failed to open the descriptor file.");
    }

    desc_file << "new(\"big.matrix.descriptor\", description = list(";
    desc_file << "sharedType = \"FileBacked\", "; // always TRUE if the
                                                  // big.matrix is file-backed
    desc_file << "filename = \"" << path_ << "\", ";

    std::string working_dir = get_path();
    // TODO : print to debug
    std::cout << "current directory : " << working_dir << "\n";

    desc_file << "dirname = \"" << working_dir << "\", ";
    desc_file << "totalRows = " << nrow_ << "L, ";
    desc_file << "totalCols = " << ncol_ << "L, ";
    // TODO : same here
    desc_file << "rowOffset = c(0, " << nrow_ << "), ";
    desc_file << "colOffset = c(0, " << ncol_ << "), ";
    desc_file << "nrow = " << nrow_ << ", ";
    desc_file << "ncol = " << ncol_ << ", ";

    // TODO : Add col names here if necessary
    desc_file << "rowNames = NULL, colNames = NULL, ";

    std::string type_ = get_type_name<T>();
    // TODO : print to debug
    // std::cout << "type in matrix: " << type_ << " \n";
    desc_file << "type = \"" << type_ << "\", ";

    desc_file << "separated = FALSE"; // from cran bigmemory, use separated
                                      // column organization of the data;
    desc_file << "))";

    desc_file.close();

    // TODO : add error codes
    // TODO : ask myself, keep void as return type ?
}
