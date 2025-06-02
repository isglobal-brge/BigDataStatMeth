/**
 * @file hdf5Utilities.hpp
 * @brief Header file for HDF5 utility functions and operations
 * @details This header file provides a comprehensive set of utilities for working
 * with HDF5 files and datasets. The implementation includes:
 * 
 * Key features:
 * - Dataset manipulation and transformation
 * - File import and export operations
 * - Data imputation and cleaning
 * - Dataset binding and splitting
 * - Dimension handling
 * 
 * Components:
 * - hdf5BindDatasets: Functions for combining multiple datasets
 * - hdf5Dims: Dimension handling utilities
 * - hdf5ImportFiles: File import and conversion utilities
 * - hdf5ImputeData: Missing data imputation functions
 * - hdf5Methods: General HDF5 operation methods
 * - hdf5ReduceDataset: Dataset reduction operations
 * - hdf5RemoveElements: Element removal utilities
 * - hdf5RemoveLowData: Low-quality data filtering
 * - hdf5SortDataset: Dataset sorting operations
 * - hdf5SplitDataset: Dataset splitting utilities
 * - hdf5ApplytoDatasets: Batch operations on datasets
 * 
 * The module provides essential utilities for managing and manipulating
 * large-scale data stored in HDF5 format.
 */
#ifndef BIGDATASTATMETH_HDF5UTILITIES_HPP
#define BIGDATASTATMETH_HDF5UTILITIES_HPP

    #include "hdf5Utilities/hdf5BindDatasets.hpp"
    #include "hdf5Utilities/hdf5Dims.hpp"
    #include "hdf5Utilities/hdf5ImportFiles.hpp"
    #include "hdf5Utilities/hdf5ImputeData.hpp"
    #include "hdf5Utilities/hdf5Methods.hpp"
    #include "hdf5Utilities/hdf5ReduceDataset.hpp"
    #include "hdf5Utilities/hdf5RemoveElements.hpp"
    #include "hdf5Utilities/hdf5RemoveLowData.hpp"
    #include "hdf5Utilities/hdf5SortDataset.hpp"
    #include "hdf5Utilities/hdf5SplitDataset.hpp"
    #include "hdf5Utilities/hdf5ApplytoDatasets.hpp"

#endif // BIGDATASTATMETH_HD5ALGEBRA_HPP