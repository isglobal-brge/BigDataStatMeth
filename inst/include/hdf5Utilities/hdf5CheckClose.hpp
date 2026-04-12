/**
 * @file hdf5CheckClose.hpp
 * @brief Safe closing helpers for HDF5 HLA datasets.
 *
 * @details
 * Header-only utilities to safely close one or more HDF5 dataset wrappers
 * used in BigDataStatMeth. The API is null-safe, idempotent (closing an
 * already-closed handle is a no-op), and does not throw. Exceptions raised
 * by the underlying HDF5 calls are caught and ignored to keep cleanup paths
 * simple and robust.
 *
 * Key features:
 *  - Close N datasets in a single call (variadic C++17).
 *  - Null-pointer checks for every argument.
 *  - No-throw guarantees for error paths.
 *  - Suitable for use in destructors and failure cleanups.
 *
 * @note Requires C++17. The dataset type must expose:
 *       `void close_file()` and `auto getDatasetptr() -> void*` (or
 *       a pointer-like type convertible to bool/nullable).
 *
 * @see H5::H5File, H5::DataSet
 * @since 0.99.0
 */

#ifndef BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP
#define BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP

#include <RcppEigen.h>
#include "H5Cpp.h"

namespace BigDataStatMeth {


    /**
     * @brief Safely close any number of dataset handles (C++17).
     *
     * @tparam Ds Parameter pack of dataset wrapper types. Each @p Ds must
     *         provide `getDatasetptr()` and `close_file()`.
     * @param ds  Zero or more pointers to dataset wrappers (nullable).
     *
     * @pre Each pointer either is `nullptr` or points to a valid dataset
     *      wrapper instance.
     * @post All non-null pointers whose `getDatasetptr()` is non-null are
     *       asked to `close_file()`. Pointers themselves are not modified.
     *
     * @exception None. The function is `noexcept`; any exception thrown by
     *            `close_file()` is caught and ignored.
     *
     * @warning Do not pass the *same* pointer more than once.
     * @warning This function does not set pointers to `nullptr`. Do that
     *          yourself if needed to avoid double closes elsewhere.
     *
     * @par Thread-safety
     * The function performs only local operations on the passed objects.
     * It is thread-safe as long as you do not pass the *same* object
     * concurrently from multiple threads.
     *
     * @par Example
     * @code
     * hdf5Dataset *a = nullptr; // dataset A
     * hdf5Dataset *b = nullptr; // dataset B
     * hdf5Dataset *c = nullptr;
     * BigDataStatMeth::checkClose_file(a, b, c);
     * a = b = c = nullptr; // optional: clear caller-owned pointers
     * @endcode
     */
    template <class... Ds>
    inline void checkClose_file(Ds*... ds) noexcept {
        ([&](auto* p) {
            try {
                if (p && p->getDatasetptr() != nullptr) p->close_file();
            } catch (...) {}
        }(ds), ...);
    }
    
    /**
     * @brief Safely close a list of dataset handles (initializer-list form).
     *
     * @param list Initializer list of dataset wrapper pointers (nullable).
     *
     * @details Convenience overload to allow brace-initializer syntax:
     *          `checkClose_file({ds1, ds2, ds3});`
     *
     * @exception None. Exceptions from `close_file()` are caught and ignored.
     *
     * @par Example
     * @code
     * BigDataStatMeth::checkClose_file({a, b, c});
     * @endcode
     */
    inline void checkClose_file(
            std::initializer_list<BigDataStatMeth::hdf5Dataset*> list) noexcept {
        for (auto* p : list) {
            try {
                if (p && p->getDatasetptr() != nullptr) p->close_file();
            } catch (...) {}
        }
    }
    
    /**
     * @brief Close all open HDF5 handles at the C library level
     * @details Uses HDF5 C API to find and close all open objects
     * (files, datasets, groups, datatypes, attributes) regardless
     * of how they were opened. Equivalent to rhdf5::h5closeAll().
     * Safe to call at any time — invalid IDs are skipped silently.
     */
    inline void closeAllHDF5Handles() {
        try {
            ssize_t n = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_ALL);
            if (n <= 0) return;
            std::vector<hid_t> ids(static_cast<size_t>(n));
            H5Fget_obj_ids(H5F_OBJ_ALL, H5F_OBJ_ALL,
                           static_cast<size_t>(n), ids.data());
            for (hid_t id : ids)
                if (H5Iis_valid(id) > 0) H5Oclose(id);
        } catch (...) {}
    }


}

#endif // BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP