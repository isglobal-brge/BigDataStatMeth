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
            // Exclude H5F_OBJ_DATATYPE: the HDF5 library keeps pre-defined types
            // (H5T_NATIVE_DOUBLE, etc.) permanently open as library-level objects.
            // Closing them corrupts the library state for subsequent operations.
            const unsigned types = H5F_OBJ_FILE | H5F_OBJ_DATASET |
                H5F_OBJ_GROUP | H5F_OBJ_ATTR;
            ssize_t n = H5Fget_obj_count(H5F_OBJ_ALL, types);
            if (n <= 0) return;
            std::vector<hid_t> ids(static_cast<size_t>(n));
            H5Fget_obj_ids(H5F_OBJ_ALL, types,
                           static_cast<size_t>(n), ids.data());
            for (hid_t id : ids)
                if (H5Iis_valid(id) > 0) H5Oclose(id);
        } catch (...) {}
    }
    
    
    /**
     * @brief Close all HDF5 handles associated with a specific file.
     *
     * Finds all open HDF5 file handles matching \p filename (by absolute path),
     * closes their associated datasets/groups/attributes, then closes the file
     * handles themselves. Pre-defined HDF5 library types are never touched.
     *
     * @param filename Absolute path to the HDF5 file.
     */
    inline void closeHDF5HandlesForFile(const std::string& filename) {
        try {
            // Get all open file handles
            ssize_t n_files = H5Fget_obj_count(H5F_OBJ_ALL, H5F_OBJ_FILE);
            if (n_files <= 0) return;
            
            std::vector<hid_t> file_ids(static_cast<size_t>(n_files));
            H5Fget_obj_ids(H5F_OBJ_ALL, H5F_OBJ_FILE,
                           static_cast<size_t>(n_files), file_ids.data());
            
            for (hid_t fid : file_ids) {
                if (H5Iis_valid(fid) <= 0) continue;
                
                // Get the absolute path HDF5 stored for this handle
                ssize_t len = H5Fget_name(fid, nullptr, 0);
                if (len <= 0) continue;
                std::string hdf5_name(static_cast<size_t>(len + 1), '\0');
                H5Fget_name(fid, &hdf5_name[0], static_cast<size_t>(len + 1));
                hdf5_name.resize(static_cast<size_t>(len));
                
                if (hdf5_name != filename) continue;
                
                // Close datasets, groups, attributes for this file only
                // Never touch H5F_OBJ_DATATYPE — pre-defined types are library-internal
                const unsigned types = H5F_OBJ_DATASET | H5F_OBJ_GROUP | H5F_OBJ_ATTR;
                ssize_t n = H5Fget_obj_count(fid, types);
                if (n > 0) {
                    std::vector<hid_t> ids(static_cast<size_t>(n));
                    H5Fget_obj_ids(fid, types, static_cast<size_t>(n), ids.data());
                    for (hid_t id : ids)
                        if (H5Iis_valid(id) > 0) H5Oclose(id);
                }
                H5Fclose(fid);
            }
        } catch (...) {}
    }


}

#endif // BIGDATASTATMETH_HDF5_CHECK_CLOSE_FILE_HPP