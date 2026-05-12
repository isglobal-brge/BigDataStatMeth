/**
 * @file hdf5_r6_matvec.cpp
 * @brief R6 wrappers for vector-matrix broadcast (sweep) and diagonal operations.
 *
 * ARCHITECTURAL RULE (roadmap section 1.4): calls header functions directly,
 * never Rcpp-exported bd* symbols.
 *
 * Functions provided:
 *   rcpp_hdf5dataset_sweep      -- broadcast vector over matrix rows/cols
 *                                 (hdf5_matrixVector_calculus)
 *   rcpp_hdf5dataset_diag_get   -- extract diagonal -> numeric vector in memory
 *                                 (getDiagonalfromMatrix)
 *   rcpp_hdf5dataset_diag_set   -- write diagonal in-place
 *                                 (setDiagonalMatrix)
 *   rcpp_hdf5dataset_diag_op    -- operate on diagonals of two HDF5 datasets
 *                                 (DiagonalOps::addDiagonals / subtractDiagonals /
 *                                  multiplyDiagonals / divideDiagonals)
 *   rcpp_hdf5dataset_diag_scale -- apply scalar to diagonal elements
 *                                 (DiagonalOps::scalarOperation)
 *
 * @note Uses Rf_error() in catch blocks to avoid macOS ARM64 double-free.
 *
 * NOTE on byrows convention (R vs HDF5):
 *   R stores matrices in column-major; HDF5 stores them transposed.
 *   The user-facing 'byrows' parameter follows base-R sweep() MARGIN:
 *     byrows=FALSE (MARGIN=2): apply vector to every R column,
 *                              vector length == ncol(matrix).
 *     byrows=TRUE  (MARGIN=1): apply vector to every R row,
 *                              vector length == nrow(matrix).
 *   The header hdf5_matrixVector_calculus() uses HDF5 coordinates, so
 *   the boolean must be FLIPPED before the call:
 *     R byrows=FALSE -> header bbyrows=TRUE
 *     R byrows=TRUE  -> header bbyrows=FALSE
 */

#include <BigDataStatMeth.hpp>
#include <sstream>
#include <ctime>
#include <cstdlib>


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_sweep
// ---------------------------------------------------------------------------

/**
 * @brief Broadcast a vector over a matrix (base-R sweep equivalent).
 *
 * Applies an element-wise binary operation between an HDF5 matrix and a
 * 1-row HDF5 vector stored as a 1xN R matrix.
 *
 *   byrows=FALSE (MARGIN=2): vector length == ncol(matrix), applied column-wise.
 *   byrows=TRUE  (MARGIN=1): vector length == nrow(matrix), applied row-wise.
 *
 * Delegates to BigDataStatMeth::hdf5_matrixVector_calculus() with FLIPPED
 * byrows (R/HDF5 transposition).
 *
 * @param ptr_mat     External pointer (SEXP) for the matrix.
 * @param ptr_vec     External pointer (SEXP) for the vector (stored as 1-row R matrix).
 * @param func        Operation: "+", "-", "*" (default), "/", "pow".
 * @param byrows      FALSE (default, MARGIN=2): vector length == ncol(matrix).
 *                    TRUE (MARGIN=1): vector length == nrow(matrix).
 * @param paral       Logical or NULL; enable OpenMP parallelisation.
 * @param threads     Integer or NULL; thread count.
 * @param compression Integer or NULL; gzip level for output (NULL = inherit from A).
 * @return Named list with elements "filename" and "path" of the result dataset.
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_sweep(SEXP ptr_mat,
                                   SEXP ptr_vec,
                                   std::string func   = "*",
                                   bool        byrows = false,
                                   Rcpp::Nullable<bool> paral       = R_NilValue,
                                   Rcpp::Nullable<int>  threads     = R_NilValue,
                                   Rcpp::Nullable<int>  compression = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "",
        Rcpp::Named("path")     = "");

    try {
        H5::Exception::dontPrint();

        auto* rawA = static_cast<BigDataStatMeth::hdf5Dataset*>( R_ExternalPtrAddr(ptr_mat));
        auto* rawB = static_cast<BigDataStatMeth::hdf5Dataset*>( R_ExternalPtrAddr(ptr_vec));

        if (rawA == nullptr || rawB == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!rawA->isOpen() || !rawB->isOpen())
            throw std::runtime_error("Dataset is closed");

        const std::string filename = rawA->getFullPath();
        const std::string groupA   = rawA->getGroup();
        const std::string nameA    = rawA->getDatasetName();
        const std::string groupB   = rawB->getGroup();
        const std::string nameB    = rawB->getDatasetName();

        // Operation map: "+":0, "-":1, "*":2, "/":3, "pow":4
        Rcpp::NumericVector oper_map = {0, 1, 2, 3, 4};
        oper_map.names() = Rcpp::CharacterVector({"+", "-", "*", "/", "pow"});
        int ioper = static_cast<int>(oper_map[oper_map.findName(func)]);

        bool bparal = (!paral.isNull()) && Rcpp::as<bool>(paral);

        // Open fresh handles (Rule 1.4: independent of the XPtr handles)
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA( new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB( new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        // Validate vector dimension using R-native coordinates.
        // Vectors are stored as 1xN R matrices, so their length is ncols_r().
        if (!byrows && dsB->ncols_r() != dsA->ncols_r())
            throw std::runtime_error( "sweep: vector length != ncol(matrix) for byrows=FALSE");
        if (byrows && dsB->ncols_r() != dsA->nrows_r())
            throw std::runtime_error( "sweep: vector length != nrow(matrix) for byrows=TRUE");
        
        
        
        int comp_level = compression.isNotNull()
                          ? Rcpp::as<int>(compression)
                          : static_cast<int>(dsA->getCompressionLevel());

        const std::string out_name = "sweep_" + nameA;
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, "OUTPUT", out_name, true));
        dsC->setCompressionLevel(comp_level);

        // FLIP byrows: R MARGIN convention is opposite to HDF5 row convention.
        // byrows=FALSE (MARGIN=2) -> header bbyrows=TRUE
        // byrows=TRUE  (MARGIN=1) -> header bbyrows=FALSE
        BigDataStatMeth::hdf5_matrixVector_calculus( dsA.get(), dsB.get(), dsC.get(), ioper, !byrows, bparal, threads);
        // BigDataStatMeth::hdf5_matrixVector_calculus( dsA.get(), dsB.get(), dsC.get(), ioper, byrows, bparal, threads);
        
        lst["filename"] = filename;
        lst["path"] = "OUTPUT/" + out_name;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_sweep (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_sweep (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_sweep: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_diag_get
// ---------------------------------------------------------------------------

/**
 * @brief Extract diagonal elements from an HDF5 matrix into memory.
 *
 * Returns the diagonal of an HDF5 dataset as an in-memory numeric vector.
 * Delegates to BigDataStatMeth::getDiagonalfromMatrix().
 *
 * @param ptr_mat External pointer (SEXP) to the HDF5 dataset.
 * @return Numeric vector of diagonal elements (length = min(nrow, ncol)).
 */
// [[Rcpp::export]]
Rcpp::NumericVector rcpp_hdf5dataset_diag_get(SEXP ptr_mat)
{
    try {
        H5::Exception::dontPrint();

        auto* raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_mat));
        if (raw == nullptr)
            Rf_error("rcpp_hdf5dataset_diag_get: invalid external pointer");
        if (!raw->isOpen())
            Rf_error("rcpp_hdf5dataset_diag_get: dataset is closed");

        // const std::string filename = raw->getFullPath();
        // const std::string group    = raw->getGroup();
        // const std::string name     = raw->getDatasetName();
        // 
        // std::unique_ptr<BigDataStatMeth::hdf5Dataset> ds(
        //     new BigDataStatMeth::hdf5Dataset(filename, group, name, false));
        // ds->openDataset();
        // 
        // if (ds->getDatasetptr() == nullptr)
        //     Rf_error("rcpp_hdf5dataset_diag_get: cannot open dataset");
        // 
        // return BigDataStatMeth::getDiagonalfromMatrix(ds.get());

        return BigDataStatMeth::getDiagonalfromMatrix(raw);
        
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_get: %s", e.what());
    }
    return Rcpp::NumericVector(0);  // unreachable
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_diag_set
// ---------------------------------------------------------------------------

/**
 * @brief Write diagonal elements into an HDF5 matrix (in-place).
 *
 * Replaces the diagonal elements of an HDF5 dataset with the supplied values.
 * The matrix must be square. Delegates to BigDataStatMeth::setDiagonalMatrix().
 *
 * @param ptr_mat External pointer (SEXP) to the HDF5 dataset.
 * @param values  Numeric vector of replacement diagonal values.
 * @return TRUE on success (invisibly).
 */
// [[Rcpp::export]]
bool rcpp_hdf5dataset_diag_set(SEXP ptr_mat, Rcpp::NumericVector values)
{
    
    try {
        H5::Exception::dontPrint();

        auto* ds = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_mat));
        if (ds == nullptr)
            Rf_error("rcpp_hdf5dataset_diag_set: invalid external pointer");
        if (!ds->isOpen())
            Rf_error("rcpp_hdf5dataset_diag_set: dataset is closed");

        BigDataStatMeth::setDiagonalMatrix(ds, values);
        return true;

    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_set: %s", e.what());
    }
    return false;

}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_diag_op
// ---------------------------------------------------------------------------

/**
 * @brief Element-wise operation on the diagonals of two HDF5 datasets.
 *
 * Operates on the diagonal elements (or diagonal vectors) of two HDF5
 * datasets and writes the result to a new dataset. Inputs may be square
 * matrices (diagonal is extracted) or 1xN / Nx1 vectors (used directly).
 * Delegates to BigDataStatMeth::DiagonalOps::addDiagonals() etc.
 *
 * @param ptr_a       External pointer (SEXP) for dataset A.
 * @param ptr_b       External pointer (SEXP) for dataset B.
 * @param op          Operation: "+" (default), "-", "*", "/".
 * @param paral       Logical or NULL; enable parallelisation.
 * @param threads     Integer or NULL; thread count.
 * @param compression Integer or NULL; gzip level (NULL = inherit from A).
 * @param outgroup   Character or NULL. Output group. Default \code{"OUTPUT"}.
 * @param outdataset Character or NULL. Output dataset name.
 *   Default \code{"diag_OP_A_B"} where OP is the operation and A, B the input names.
 * @return Named list with elements "filename" and "path" of the result dataset.
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_diag_op(SEXP ptr_a,
                                     SEXP ptr_b,
                                     std::string op        = "+",
                                     Rcpp::Nullable<bool> paral       = R_NilValue,
                                     Rcpp::Nullable<int>  threads     = R_NilValue,
                                     Rcpp::Nullable<int>  compression = R_NilValue,
                                     Rcpp::Nullable<std::string> outgroup    = R_NilValue,
                                     Rcpp::Nullable<std::string> outdataset  = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "",
        Rcpp::Named("path")     = "");

    try {
        H5::Exception::dontPrint();

        auto* rawA = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_a));
        auto* rawB = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_b));

        if (rawA == nullptr || rawB == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!rawA->isOpen() || !rawB->isOpen())
            throw std::runtime_error("Dataset is closed");
        if (rawA->getFullPath() != rawB->getFullPath())
            throw std::runtime_error(
                "diag_op: both datasets must be in the same HDF5 file");

        const std::string filename = rawA->getFullPath();
        const std::string groupA   = rawA->getGroup();
        const std::string nameA    = rawA->getDatasetName();
        const std::string groupB   = rawB->getGroup();
        const std::string nameB    = rawB->getDatasetName();

        bool bparal = (!paral.isNull()) && Rcpp::as<bool>(paral);
        int comp_level = compression.isNotNull()
                          ? Rcpp::as<int>(compression)
                          : static_cast<int>(rawA->getCompressionLevel());

        // Fresh handles
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsA(
            new BigDataStatMeth::hdf5Dataset(filename, groupA, nameA, false));
        dsA->openDataset();
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsB(
            new BigDataStatMeth::hdf5Dataset(filename, groupB, nameB, false));
        dsB->openDataset();

        if (dsA->getDatasetptr() == nullptr || dsB->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open datasets");

        const std::string out_grp  = outgroup.isNull()
            ? std::string("OUTPUT")
                : Rcpp::as<std::string>(outgroup);
        const std::string out_name = outdataset.isNull()
            ? ("diag_" + op + "_" + nameA + "_" + nameB)
            : Rcpp::as<std::string>(outdataset);
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsC(
                new BigDataStatMeth::hdf5Dataset(filename, out_grp, out_name, true));
        dsC->setCompressionLevel(comp_level);

        if      (op == "+") BigDataStatMeth::DiagonalOps::addDiagonals     (dsA.get(), dsB.get(), dsC.get(), "new", bparal, threads);
        else if (op == "-") BigDataStatMeth::DiagonalOps::subtractDiagonals (dsA.get(), dsB.get(), dsC.get(), "new", bparal, threads);
        else if (op == "*") BigDataStatMeth::DiagonalOps::multiplyDiagonals (dsA.get(), dsB.get(), dsC.get(), "new", bparal, threads);
        else if (op == "/") BigDataStatMeth::DiagonalOps::divideDiagonals   (dsA.get(), dsB.get(), dsC.get(), "new", bparal, threads);
        else throw std::runtime_error("diag_op: op must be one of '+', '-', '*', '/'");

        lst["filename"] = filename;
        lst["path"] = out_grp + "/" + out_name;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_op (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_op (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_op: %s", e.what());
    }

    return lst;
}


// ---------------------------------------------------------------------------
// rcpp_hdf5dataset_diag_scale
// ---------------------------------------------------------------------------

/**
 * @brief Apply a scalar operation to the diagonal elements of an HDF5 matrix.
 *
 * Applies a scalar operation (add / subtract / multiply / divide) to the
 * diagonal elements of a square HDF5 matrix. Off-diagonal elements are
 * unchanged. Result is written to a new dataset.
 * Delegates to BigDataStatMeth::DiagonalOps::scalarOperation().
 * Operation codes: 0=add, 1=subtract, 2=multiply, 3=divide.
 *
 * @param ptr_mat     External pointer (SEXP) to the input matrix.
 * @param scalar      Numeric scalar value.
 * @param op_code     Integer operation code (0=add, 1=subtract, 2=multiply, 3=divide).
 * @param paral       Logical or NULL; enable parallelisation.
 * @param threads     Integer or NULL; thread count.
 * @param compression Integer or NULL; gzip level (NULL = inherit from input).
 * @param outgroup   Character or NULL. Output group. Default \code{"OUTPUT"}.
 * @param outdataset Character or NULL. Output dataset name.
 *  Default \code{"diagscale_OP_A"} where OP is the operation and A the input name.
 * @return Named list with elements "filename" and "path" of the result dataset.
 */
// [[Rcpp::export]]
Rcpp::List rcpp_hdf5dataset_diag_scale(SEXP ptr_mat,
                                        double      scalar,
                                        int         op_code     = 2,
                                        Rcpp::Nullable<bool> paral       = R_NilValue,
                                        Rcpp::Nullable<int>  threads     = R_NilValue,
                                        Rcpp::Nullable<int>  compression = R_NilValue,
                                        Rcpp::Nullable<std::string> outgroup   = R_NilValue,
                                        Rcpp::Nullable<std::string> outdataset = R_NilValue)
{
    Rcpp::List lst = Rcpp::List::create(
        Rcpp::Named("filename") = "",
        Rcpp::Named("path")     = "");

    try {
        H5::Exception::dontPrint();

        auto* raw = static_cast<BigDataStatMeth::hdf5Dataset*>(
            R_ExternalPtrAddr(ptr_mat));
        if (raw == nullptr)
            throw std::runtime_error("Invalid external pointer");
        if (!raw->isOpen())
            throw std::runtime_error("Dataset is closed");

        if (op_code < 0 || op_code > 3)
            throw std::runtime_error(
                "diag_scale: op_code must be 0=add, 1=subtract, 2=multiply, 3=divide");

        const std::string filename = raw->getFullPath();
        const std::string group    = raw->getGroup();
        const std::string name     = raw->getDatasetName();

        bool bparal = (!paral.isNull()) && Rcpp::as<bool>(paral);
        int comp_level = compression.isNotNull()
                          ? Rcpp::as<int>(compression)
                          : static_cast<int>(raw->getCompressionLevel());

        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsIn(
            new BigDataStatMeth::hdf5Dataset(filename, group, name, false));
        dsIn->openDataset();

        if (dsIn->getDatasetptr() == nullptr)
            throw std::runtime_error("Failed to open dataset");

        static const char* op_names[] = {"add","sub","mul","div"};
        const std::string out_grp  = outgroup.isNull()
            ? std::string("OUTPUT")
                : Rcpp::as<std::string>(outgroup);
        const std::string out_name = outdataset.isNull()
            ? ("diagscale_" + std::string(op_names[op_code]) + "_" + name)
            : Rcpp::as<std::string>(outdataset);
        std::unique_ptr<BigDataStatMeth::hdf5Dataset> dsOut(
                new BigDataStatMeth::hdf5Dataset(filename, out_grp, out_name, true));
        
        dsOut->setCompressionLevel(comp_level);

        
        const hsize_t n     = dsIn->nrows_r();   // R: filas lógicas
        const hsize_t m     = dsIn->ncols_r();   // R: columnas lógicas
        const bool is_vec   = (n == 1 || m == 1);
        
        if (is_vec) {
            const hsize_t sz = std::max(n, m);
            std::vector<double> buf(sz);
            std::vector<hsize_t> st = {1,1}, blk = {1,1};
            
            // HDF5 almacena transpuesto: 1×200 en R → 200×1 en HDF5
            const hsize_t hdf5_rows = dsIn->nrows();  // HDF5-nativo
            const hsize_t hdf5_cols = dsIn->ncols();  // HDF5-nativo
            dsIn->readDatasetBlock({0,0}, {hdf5_rows, hdf5_cols}, st, blk, buf.data());
            
            for (hsize_t i = 0; i < sz; ++i) {
                switch (op_code) {
                case 0: buf[i] += scalar; break;
                case 1: buf[i] -= scalar; break;
                case 2: buf[i] *= scalar; break;
                case 3: buf[i] /= scalar; break;
                }
            }
            
            // Crear output con mismas dimensiones R que input
            dsOut->createDataset(n, m, "real");  // createDataset usa coordenadas R
            dsOut->writeDatasetBlock(buf, {0,0}, {hdf5_rows, hdf5_cols}, st, blk);
        } else {
            BigDataStatMeth::DiagonalOps::scalarOperation(
                dsIn.get(), dsOut.get(), scalar, op_code, "new", bparal, threads);
        }
        
        //.. 17/03/2026  ..// BigDataStatMeth::DiagonalOps::scalarOperation(dsIn.get(), dsOut.get(), scalar, op_code, "new", bparal, threads);

        // // Create output with same R dimensions as input (full matrix copy)
        // const hsize_t src_nr = dsIn->nrows();
        // const hsize_t src_nc = dsIn->ncols();
        // dsOut->createDataset(src_nc, src_nr, "real");  // createDataset(R_rows, R_cols)
        // 
        // // Block-copy input matrix to output
        // {
        //     const hsize_t COPY_BLOCK = 500;
        //     std::vector<hsize_t> st = {1, 1}, blk = {1, 1};
        //     for (hsize_t row = 0; row < src_nr; row += COPY_BLOCK) {
        //         hsize_t toread = std::min(COPY_BLOCK, src_nr - row);
        //         std::vector<double> buf(toread * src_nc);
        //         dsIn->readDatasetBlock({row, 0}, {toread, src_nc}, st, blk, buf.data());
        //         dsOut->writeDatasetBlock(buf, {row, 0}, {toread, src_nc}, st, blk);
        //     }
        // }
        // 
        // // Get diagonal, apply scalar in memory, write back
        // Rcpp::NumericVector diag_vals = BigDataStatMeth::getDiagonalfromMatrix(dsIn.get());
        // for (int i = 0; i < diag_vals.size(); ++i) {
        //     switch (op_code) {
        //     case 0: diag_vals[i] += scalar; break;
        //     case 1: diag_vals[i] -= scalar; break;
        //     case 2: diag_vals[i] *= scalar; break;
        //     case 3: diag_vals[i] /= scalar; break;
        //     }
        // }
        // BigDataStatMeth::setDiagonalMatrix(dsOut.get(), diag_vals);
        
        lst["filename"] = filename;
        lst["path"] = out_grp + "/" + out_name;

    } catch (H5::FileIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_scale (File IException): %s",
                 e.getCDetailMsg());
    } catch (H5::DataSetIException& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_scale (DataSet IException): %s",
                 e.getCDetailMsg());
    } catch (std::exception& e) {
        Rf_error("c++ exception rcpp_hdf5dataset_diag_scale: %s", e.what());
    }

    return lst;
}
