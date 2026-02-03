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
    T* operator->() const { return ptr_; }
    T& operator*() const { return *ptr_; }
    
    // Non-copyable
    HDF5Handle(const HDF5Handle&) = delete;
    HDF5Handle& operator=(const HDF5Handle&) = delete;
};

} // namespace BigDataStatMeth

#endif