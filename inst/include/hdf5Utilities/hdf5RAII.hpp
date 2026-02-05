#ifndef BIGDATASTATMETH_HDF5_RAII_HPP
#define BIGDATASTATMETH_HDF5_RAII_HPP

namespace BigDataStatMeth {

template<typename T>
class HDF5Handle {
private:
    T* ptr_;
public:
    explicit HDF5Handle(T* p = nullptr) : ptr_(p) {}
    
    ~HDF5Handle() {
        if (ptr_) {
            delete ptr_;
            ptr_ = nullptr;
        }
    }
    
    T* get() const { return ptr_; }
    T* operator->() const { 
        if (!ptr_) {
            Rcpp::stop("Attempting to access null HDF5Handle");
        }
        return ptr_; 
    }
    
    T& operator*() const { 
        if (!ptr_) {
            Rcpp::stop("Attempting to dereference null HDF5Handle");
        }
        return *ptr_; 
    }
    
    explicit operator bool() const { return ptr_ != nullptr; } 
    
    void reset(T* p = nullptr) {
        if (ptr_) delete ptr_;
        ptr_ = p;
    }
    
    // Non-copyable
    HDF5Handle(const HDF5Handle&) = delete;
    HDF5Handle& operator=(const HDF5Handle&) = delete;
};

// Aliases - declared here, but only usable after the classes are included
// BigDataStatMeth.hpp incluye este archivo después de las clases
using hdf5FileHandle           = HDF5Handle<hdf5File>;
using hdf5GroupHandle          = HDF5Handle<hdf5Group>;
using hdf5DatasetHandle        = HDF5Handle<hdf5Dataset>;
using hdf5DatasetInternalHandle = HDF5Handle<hdf5DatasetInternal>;
using hdf5DatasetDimsHandle     = HDF5Handle<hdf5Dims>;

} // namespace BigDataStatMeth

#endif