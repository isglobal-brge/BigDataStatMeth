// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppEigen.h>
#include <Rcpp.h>

using namespace Rcpp;

// bdRemovelowdata
Rcpp::RObject bdRemovelowdata(std::string filename, std::string group, std::string dataset, std::string outgroup, std::string outdataset, Rcpp::Nullable<double> pcent, Rcpp::Nullable<bool> SNPincols);
RcppExport SEXP _BigDataStatMeth_bdRemovelowdata(SEXP filenameSEXP, SEXP groupSEXP, SEXP datasetSEXP, SEXP outgroupSEXP, SEXP outdatasetSEXP, SEXP pcentSEXP, SEXP SNPincolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< std::string >::type outgroup(outgroupSEXP);
    Rcpp::traits::input_parameter< std::string >::type outdataset(outdatasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type pcent(pcentSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type SNPincols(SNPincolsSEXP);
    rcpp_result_gen = Rcpp::wrap(bdRemovelowdata(filename, group, dataset, outgroup, outdataset, pcent, SNPincols));
    return rcpp_result_gen;
END_RCPP
}
// Crossprod_hdf5
Rcpp::RObject Crossprod_hdf5(std::string filename, const std::string group, std::string A, Rcpp::Nullable<std::string> groupB, Rcpp::Nullable<std::string> B, Rcpp::Nullable<int> block_size, Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads, Rcpp::Nullable<double> mixblock_size, Rcpp::Nullable<std::string> outgroup);
RcppExport SEXP _BigDataStatMeth_Crossprod_hdf5(SEXP filenameSEXP, SEXP groupSEXP, SEXP ASEXP, SEXP groupBSEXP, SEXP BSEXP, SEXP block_sizeSEXP, SEXP paralSEXP, SEXP threadsSEXP, SEXP mixblock_sizeSEXP, SEXP outgroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type groupB(groupBSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type paral(paralSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type mixblock_size(mixblock_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outgroup(outgroupSEXP);
    rcpp_result_gen = Rcpp::wrap(Crossprod_hdf5(filename, group, A, groupB, B, block_size, paral, threads, mixblock_size, outgroup));
    return rcpp_result_gen;
END_RCPP
}
// blockmult_hdf5
Rcpp::RObject blockmult_hdf5(std::string filename, const std::string group, std::string A, std::string B, Rcpp::Nullable<int> block_size, Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads, Rcpp::Nullable<double> mixblock_size, Rcpp::Nullable<std::string> outgroup);
RcppExport SEXP _BigDataStatMeth_blockmult_hdf5(SEXP filenameSEXP, SEXP groupSEXP, SEXP ASEXP, SEXP BSEXP, SEXP block_sizeSEXP, SEXP paralSEXP, SEXP threadsSEXP, SEXP mixblock_sizeSEXP, SEXP outgroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::string >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type paral(paralSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type mixblock_size(mixblock_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outgroup(outgroupSEXP);
    rcpp_result_gen = Rcpp::wrap(blockmult_hdf5(filename, group, A, B, block_size, paral, threads, mixblock_size, outgroup));
    return rcpp_result_gen;
END_RCPP
}
// blockmult_sparse_hdf5
Rcpp::RObject blockmult_sparse_hdf5(std::string filename, const std::string group, std::string A, std::string B, Rcpp::Nullable<std::string> outgroup);
RcppExport SEXP _BigDataStatMeth_blockmult_sparse_hdf5(SEXP filenameSEXP, SEXP groupSEXP, SEXP ASEXP, SEXP BSEXP, SEXP outgroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type A(ASEXP);
    Rcpp::traits::input_parameter< std::string >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outgroup(outgroupSEXP);
    rcpp_result_gen = Rcpp::wrap(blockmult_sparse_hdf5(filename, group, A, B, outgroup));
    return rcpp_result_gen;
END_RCPP
}
// tCrossprod_hdf5
Rcpp::RObject tCrossprod_hdf5(std::string filename, const std::string group, std::string A, Rcpp::Nullable<std::string> groupB, Rcpp::Nullable<std::string> B, Rcpp::Nullable<int> block_size, Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads, Rcpp::Nullable<double> mixblock_size, Rcpp::Nullable<std::string> outgroup);
RcppExport SEXP _BigDataStatMeth_tCrossprod_hdf5(SEXP filenameSEXP, SEXP groupSEXP, SEXP ASEXP, SEXP groupBSEXP, SEXP BSEXP, SEXP block_sizeSEXP, SEXP paralSEXP, SEXP threadsSEXP, SEXP mixblock_sizeSEXP, SEXP outgroupSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type groupB(groupBSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type paral(paralSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type mixblock_size(mixblock_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outgroup(outgroupSEXP);
    rcpp_result_gen = Rcpp::wrap(tCrossprod_hdf5(filename, group, A, groupB, B, block_size, paral, threads, mixblock_size, outgroup));
    return rcpp_result_gen;
END_RCPP
}
// Import_text_to_HDF5
int Import_text_to_HDF5(Rcpp::CharacterVector filename, std::string outputfile, std::string outGroup, std::string outDataset, Rcpp::Nullable<std::string> sep, Rcpp::Nullable<bool> header, Rcpp::Nullable<bool> rownames, Rcpp::Nullable<bool> overwrite);
RcppExport SEXP _BigDataStatMeth_Import_text_to_HDF5(SEXP filenameSEXP, SEXP outputfileSEXP, SEXP outGroupSEXP, SEXP outDatasetSEXP, SEXP sepSEXP, SEXP headerSEXP, SEXP rownamesSEXP, SEXP overwriteSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::CharacterVector >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type outputfile(outputfileSEXP);
    Rcpp::traits::input_parameter< std::string >::type outGroup(outGroupSEXP);
    Rcpp::traits::input_parameter< std::string >::type outDataset(outDatasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type sep(sepSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type header(headerSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type rownames(rownamesSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type overwrite(overwriteSEXP);
    rcpp_result_gen = Rcpp::wrap(Import_text_to_HDF5(filename, outputfile, outGroup, outDataset, sep, header, rownames, overwrite));
    return rcpp_result_gen;
END_RCPP
}
// bdImpute_snps_hdf5
Rcpp::RObject bdImpute_snps_hdf5(std::string filename, std::string group, std::string dataset, Rcpp::Nullable<std::string> outgroup, Rcpp::Nullable<std::string> outdataset, Rcpp::Nullable<bool> bycols);
RcppExport SEXP _BigDataStatMeth_bdImpute_snps_hdf5(SEXP filenameSEXP, SEXP groupSEXP, SEXP datasetSEXP, SEXP outgroupSEXP, SEXP outdatasetSEXP, SEXP bycolsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outgroup(outgroupSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outdataset(outdatasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bycols(bycolsSEXP);
    rcpp_result_gen = Rcpp::wrap(bdImpute_snps_hdf5(filename, group, dataset, outgroup, outdataset, bycols));
    return rcpp_result_gen;
END_RCPP
}
// Normalize_hdf5
Rcpp::RObject Normalize_hdf5(std::string filename, const std::string group, std::string dataset, Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale, Rcpp::Nullable<int> wsize);
RcppExport SEXP _BigDataStatMeth_Normalize_hdf5(SEXP filenameSEXP, SEXP groupSEXP, SEXP datasetSEXP, SEXP bcenterSEXP, SEXP bscaleSEXP, SEXP wsizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< const std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bcenter(bcenterSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bscale(bscaleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type wsize(wsizeSEXP);
    rcpp_result_gen = Rcpp::wrap(Normalize_hdf5(filename, group, dataset, bcenter, bscale, wsize));
    return rcpp_result_gen;
END_RCPP
}
// bdPCA_hdf5
Rcpp::RObject bdPCA_hdf5(std::string filename, std::string group, std::string dataset, Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale, Rcpp::Nullable<int> k, Rcpp::Nullable<int> q, Rcpp::Nullable<bool> force, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_bdPCA_hdf5(SEXP filenameSEXP, SEXP groupSEXP, SEXP datasetSEXP, SEXP bcenterSEXP, SEXP bscaleSEXP, SEXP kSEXP, SEXP qSEXP, SEXP forceSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bcenter(bcenterSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bscale(bscaleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type force(forceSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bdPCA_hdf5(filename, group, dataset, bcenter, bscale, k, q, force, threads));
    return rcpp_result_gen;
END_RCPP
}
// Normalize_Data
Rcpp::RObject Normalize_Data(Rcpp::RObject& X, Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale);
RcppExport SEXP _BigDataStatMeth_Normalize_Data(SEXP XSEXP, SEXP bcenterSEXP, SEXP bscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bcenter(bcenterSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bscale(bscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(Normalize_Data(X, bcenter, bscale));
    return rcpp_result_gen;
END_RCPP
}
// bdMLR_MR
Rcpp::RObject bdMLR_MR(Rcpp::RObject X, Rcpp::RObject y, int blocks, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_bdMLR_MR(SEXP XSEXP, SEXP ySEXP, SEXP blocksSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type blocks(blocksSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bdMLR_MR(X, y, blocks, threads));
    return rcpp_result_gen;
END_RCPP
}
// bdCrossprod_generic
Eigen::MatrixXd bdCrossprod_generic(Rcpp::RObject A, Rcpp::Nullable<Rcpp::RObject> B, Rcpp::Nullable<bool> transposed, Rcpp::Nullable<int> block_size, Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_bdCrossprod_generic(SEXP ASEXP, SEXP BSEXP, SEXP transposedSEXP, SEXP block_sizeSEXP, SEXP paralSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::RObject> >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type transposed(transposedSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type paral(paralSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bdCrossprod_generic(A, B, transposed, block_size, paral, threads));
    return rcpp_result_gen;
END_RCPP
}
// bdwproduct
Eigen::MatrixXd bdwproduct(Rcpp::RObject X, Rcpp::RObject w, std::string op);
RcppExport SEXP _BigDataStatMeth_bdwproduct(SEXP XSEXP, SEXP wSEXP, SEXP opSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type w(wSEXP);
    Rcpp::traits::input_parameter< std::string >::type op(opSEXP);
    rcpp_result_gen = Rcpp::wrap(bdwproduct(X, w, op));
    return rcpp_result_gen;
END_RCPP
}
// bdScalarwproduct
Eigen::MatrixXd bdScalarwproduct(Rcpp::RObject A, double w, std::string op);
RcppExport SEXP _BigDataStatMeth_bdScalarwproduct(SEXP ASEXP, SEXP wSEXP, SEXP opSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type A(ASEXP);
    Rcpp::traits::input_parameter< double >::type w(wSEXP);
    Rcpp::traits::input_parameter< std::string >::type op(opSEXP);
    rcpp_result_gen = Rcpp::wrap(bdScalarwproduct(A, w, op));
    return rcpp_result_gen;
END_RCPP
}
// blockmult_sparse
Rcpp::RObject blockmult_sparse(Rcpp::RObject A, Rcpp::RObject B, Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_blockmult_sparse(SEXP ASEXP, SEXP BSEXP, SEXP paralSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type B(BSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type paral(paralSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(blockmult_sparse(A, B, paral, threads));
    return rcpp_result_gen;
END_RCPP
}
// blockmult
Rcpp::List blockmult(Rcpp::RObject a, Rcpp::RObject b, Rcpp::Nullable<int> block_size, Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads, Rcpp::Nullable<double> bigmatrix, Rcpp::Nullable<double> mixblock_size, Rcpp::Nullable<std::string> outfile, Rcpp::Nullable<bool> onmemory);
RcppExport SEXP _BigDataStatMeth_blockmult(SEXP aSEXP, SEXP bSEXP, SEXP block_sizeSEXP, SEXP paralSEXP, SEXP threadsSEXP, SEXP bigmatrixSEXP, SEXP mixblock_sizeSEXP, SEXP outfileSEXP, SEXP onmemorySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type a(aSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type b(bSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type paral(paralSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type bigmatrix(bigmatrixSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type mixblock_size(mixblock_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type outfile(outfileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type onmemory(onmemorySEXP);
    rcpp_result_gen = Rcpp::wrap(blockmult(a, b, block_size, paral, threads, bigmatrix, mixblock_size, outfile, onmemory));
    return rcpp_result_gen;
END_RCPP
}
// CholFactor
Eigen::MatrixXd CholFactor(Rcpp::RObject a);
RcppExport SEXP _BigDataStatMeth_CholFactor(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(CholFactor(a));
    return rcpp_result_gen;
END_RCPP
}
// CholSolve
Eigen::MatrixXd CholSolve(Rcpp::RObject a, Rcpp::RObject b);
RcppExport SEXP _BigDataStatMeth_CholSolve(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type a(aSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(CholSolve(a, b));
    return rcpp_result_gen;
END_RCPP
}
// inversechol_par
Eigen::MatrixXd inversechol_par(Rcpp::RObject a, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_inversechol_par(SEXP aSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type a(aSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(inversechol_par(a, threads));
    return rcpp_result_gen;
END_RCPP
}
// partCrossProd
Rcpp::RObject partCrossProd(Rcpp::RObject X);
RcppExport SEXP _BigDataStatMeth_partCrossProd(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(partCrossProd(X));
    return rcpp_result_gen;
END_RCPP
}
// partCrossProd_block
Rcpp::RObject partCrossProd_block(Rcpp::RObject X);
RcppExport SEXP _BigDataStatMeth_partCrossProd_block(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(partCrossProd_block(X));
    return rcpp_result_gen;
END_RCPP
}
// parCrossProd
Rcpp::RObject parCrossProd(Rcpp::RObject X);
RcppExport SEXP _BigDataStatMeth_parCrossProd(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(parCrossProd(X));
    return rcpp_result_gen;
END_RCPP
}
// parCrossProd_block
Rcpp::RObject parCrossProd_block(Rcpp::RObject X);
RcppExport SEXP _BigDataStatMeth_parCrossProd_block(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(parCrossProd_block(X));
    return rcpp_result_gen;
END_RCPP
}
// partCrossProdEigen
Rcpp::RObject partCrossProdEigen(Rcpp::RObject X);
RcppExport SEXP _BigDataStatMeth_partCrossProdEigen(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(partCrossProdEigen(X));
    return rcpp_result_gen;
END_RCPP
}
// parxwxt
Rcpp::RObject parxwxt(Rcpp::RObject X, Rcpp::RObject W);
RcppExport SEXP _BigDataStatMeth_parxwxt(SEXP XSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(parxwxt(X, W));
    return rcpp_result_gen;
END_RCPP
}
// parxtwx
Rcpp::RObject parxtwx(Rcpp::RObject X, Rcpp::RObject W);
RcppExport SEXP _BigDataStatMeth_parxtwx(SEXP XSEXP, SEXP WSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type W(WSEXP);
    rcpp_result_gen = Rcpp::wrap(parxtwx(X, W));
    return rcpp_result_gen;
END_RCPP
}
// parXy
Rcpp::RObject parXy(Rcpp::RObject X, Rcpp::RObject Y);
RcppExport SEXP _BigDataStatMeth_parXy(SEXP XSEXP, SEXP YSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type Y(YSEXP);
    rcpp_result_gen = Rcpp::wrap(parXy(X, Y));
    return rcpp_result_gen;
END_RCPP
}
// bdpseudoinv
Rcpp::RObject bdpseudoinv(const Rcpp::RObject& X);
RcppExport SEXP _BigDataStatMeth_bdpseudoinv(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(bdpseudoinv(X));
    return rcpp_result_gen;
END_RCPP
}
// bdQR
Rcpp::RObject bdQR(const Rcpp::RObject& X, Rcpp::Nullable<bool> thin);
RcppExport SEXP _BigDataStatMeth_bdQR(SEXP XSEXP, SEXP thinSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type thin(thinSEXP);
    rcpp_result_gen = Rcpp::wrap(bdQR(X, thin));
    return rcpp_result_gen;
END_RCPP
}
// bddtrsm
Rcpp::RObject bddtrsm(Rcpp::RObject R, Rcpp::RObject Z, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_bddtrsm(SEXP RSEXP, SEXP ZSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type R(RSEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bddtrsm(R, Z, threads));
    return rcpp_result_gen;
END_RCPP
}
// Create_HDF5_matrix_file
Rcpp::RObject Create_HDF5_matrix_file(std::string filename, RObject object, Rcpp::Nullable<std::string> group, Rcpp::Nullable<std::string> dataset, Rcpp::Nullable<bool> transp);
RcppExport SEXP _BigDataStatMeth_Create_HDF5_matrix_file(SEXP filenameSEXP, SEXP objectSEXP, SEXP groupSEXP, SEXP datasetSEXP, SEXP transpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< RObject >::type object(objectSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type group(groupSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<std::string> >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type transp(transpSEXP);
    rcpp_result_gen = Rcpp::wrap(Create_HDF5_matrix_file(filename, object, group, dataset, transp));
    return rcpp_result_gen;
END_RCPP
}
// Add_HDF5_matrix
Rcpp::RObject Add_HDF5_matrix(RObject object, std::string filename, std::string group, std::string dataset, Rcpp::Nullable<bool> transp);
RcppExport SEXP _BigDataStatMeth_Add_HDF5_matrix(SEXP objectSEXP, SEXP filenameSEXP, SEXP groupSEXP, SEXP datasetSEXP, SEXP transpSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< RObject >::type object(objectSEXP);
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type group(groupSEXP);
    Rcpp::traits::input_parameter< std::string >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type transp(transpSEXP);
    rcpp_result_gen = Rcpp::wrap(Add_HDF5_matrix(object, filename, group, dataset, transp));
    return rcpp_result_gen;
END_RCPP
}
// Remove_HDF5_element
Rcpp::RObject Remove_HDF5_element(std::string filename, std::string element);
RcppExport SEXP _BigDataStatMeth_Remove_HDF5_element(SEXP filenameSEXP, SEXP elementSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type filename(filenameSEXP);
    Rcpp::traits::input_parameter< std::string >::type element(elementSEXP);
    rcpp_result_gen = Rcpp::wrap(Remove_HDF5_element(filename, element));
    return rcpp_result_gen;
END_RCPP
}
// bdSolve
Rcpp::RObject bdSolve(const Rcpp::RObject A, const Rcpp::RObject B);
RcppExport SEXP _BigDataStatMeth_bdSolve(SEXP ASEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::RObject >::type A(ASEXP);
    Rcpp::traits::input_parameter< const Rcpp::RObject >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(bdSolve(A, B));
    return rcpp_result_gen;
END_RCPP
}
// bdInvCholesky
Eigen::MatrixXd bdInvCholesky(const Rcpp::RObject& X);
RcppExport SEXP _BigDataStatMeth_bdInvCholesky(SEXP XSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type X(XSEXP);
    rcpp_result_gen = Rcpp::wrap(bdInvCholesky(X));
    return rcpp_result_gen;
END_RCPP
}
// bdSVD
Rcpp::RObject bdSVD(const Rcpp::RObject& X, Rcpp::Nullable<int> k, Rcpp::Nullable<int> nev, Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale);
RcppExport SEXP _BigDataStatMeth_bdSVD(SEXP XSEXP, SEXP kSEXP, SEXP nevSEXP, SEXP bcenterSEXP, SEXP bscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type nev(nevSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bcenter(bcenterSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bscale(bscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(bdSVD(X, k, nev, bcenter, bscale));
    return rcpp_result_gen;
END_RCPP
}
// bdSVD_hdf5
Rcpp::RObject bdSVD_hdf5(const Rcpp::RObject& file, Rcpp::Nullable<CharacterVector> group, Rcpp::Nullable<CharacterVector> dataset, Rcpp::Nullable<int> k, Rcpp::Nullable<int> q, Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_bdSVD_hdf5(SEXP fileSEXP, SEXP groupSEXP, SEXP datasetSEXP, SEXP kSEXP, SEXP qSEXP, SEXP bcenterSEXP, SEXP bscaleSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::RObject& >::type file(fileSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<CharacterVector> >::type group(groupSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<CharacterVector> >::type dataset(datasetSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type k(kSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type q(qSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bcenter(bcenterSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bscale(bscaleSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(bdSVD_hdf5(file, group, dataset, k, q, bcenter, bscale, threads));
    return rcpp_result_gen;
END_RCPP
}
// bdSVD_lapack
Rcpp::RObject bdSVD_lapack(Rcpp::RObject X, Rcpp::Nullable<bool> bcenter, Rcpp::Nullable<bool> bscale);
RcppExport SEXP _BigDataStatMeth_bdSVD_lapack(SEXP XSEXP, SEXP bcenterSEXP, SEXP bscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bcenter(bcenterSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type bscale(bscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(bdSVD_lapack(X, bcenter, bscale));
    return rcpp_result_gen;
END_RCPP
}
// tCrossprod_Weighted
Rcpp::RObject tCrossprod_Weighted(Rcpp::RObject A, Rcpp::RObject W, Rcpp::Nullable<int> block_size, Rcpp::Nullable<bool> paral, Rcpp::Nullable<int> threads);
RcppExport SEXP _BigDataStatMeth_tCrossprod_Weighted(SEXP ASEXP, SEXP WSEXP, SEXP block_sizeSEXP, SEXP paralSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::RObject >::type A(ASEXP);
    Rcpp::traits::input_parameter< Rcpp::RObject >::type W(WSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type block_size(block_sizeSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<bool> >::type paral(paralSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<int> >::type threads(threadsSEXP);
    rcpp_result_gen = Rcpp::wrap(tCrossprod_Weighted(A, W, block_size, paral, threads));
    return rcpp_result_gen;
END_RCPP
}
// parallelVectorSum
double parallelVectorSum(Rcpp::NumericVector x);
RcppExport SEXP _BigDataStatMeth_parallelVectorSum(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelVectorSum(x));
    return rcpp_result_gen;
END_RCPP
}
// parallelpow2
Rcpp::NumericVector parallelpow2(Rcpp::NumericVector x);
RcppExport SEXP _BigDataStatMeth_parallelpow2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(parallelpow2(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BigDataStatMeth_bdRemovelowdata", (DL_FUNC) &_BigDataStatMeth_bdRemovelowdata, 7},
    {"_BigDataStatMeth_Crossprod_hdf5", (DL_FUNC) &_BigDataStatMeth_Crossprod_hdf5, 10},
    {"_BigDataStatMeth_blockmult_hdf5", (DL_FUNC) &_BigDataStatMeth_blockmult_hdf5, 9},
    {"_BigDataStatMeth_blockmult_sparse_hdf5", (DL_FUNC) &_BigDataStatMeth_blockmult_sparse_hdf5, 5},
    {"_BigDataStatMeth_tCrossprod_hdf5", (DL_FUNC) &_BigDataStatMeth_tCrossprod_hdf5, 10},
    {"_BigDataStatMeth_Import_text_to_HDF5", (DL_FUNC) &_BigDataStatMeth_Import_text_to_HDF5, 8},
    {"_BigDataStatMeth_bdImpute_snps_hdf5", (DL_FUNC) &_BigDataStatMeth_bdImpute_snps_hdf5, 6},
    {"_BigDataStatMeth_Normalize_hdf5", (DL_FUNC) &_BigDataStatMeth_Normalize_hdf5, 6},
    {"_BigDataStatMeth_bdPCA_hdf5", (DL_FUNC) &_BigDataStatMeth_bdPCA_hdf5, 9},
    {"_BigDataStatMeth_Normalize_Data", (DL_FUNC) &_BigDataStatMeth_Normalize_Data, 3},
    {"_BigDataStatMeth_bdMLR_MR", (DL_FUNC) &_BigDataStatMeth_bdMLR_MR, 4},
    {"_BigDataStatMeth_bdCrossprod_generic", (DL_FUNC) &_BigDataStatMeth_bdCrossprod_generic, 6},
    {"_BigDataStatMeth_bdwproduct", (DL_FUNC) &_BigDataStatMeth_bdwproduct, 3},
    {"_BigDataStatMeth_bdScalarwproduct", (DL_FUNC) &_BigDataStatMeth_bdScalarwproduct, 3},
    {"_BigDataStatMeth_blockmult_sparse", (DL_FUNC) &_BigDataStatMeth_blockmult_sparse, 4},
    {"_BigDataStatMeth_blockmult", (DL_FUNC) &_BigDataStatMeth_blockmult, 9},
    {"_BigDataStatMeth_CholFactor", (DL_FUNC) &_BigDataStatMeth_CholFactor, 1},
    {"_BigDataStatMeth_CholSolve", (DL_FUNC) &_BigDataStatMeth_CholSolve, 2},
    {"_BigDataStatMeth_inversechol_par", (DL_FUNC) &_BigDataStatMeth_inversechol_par, 2},
    {"_BigDataStatMeth_partCrossProd", (DL_FUNC) &_BigDataStatMeth_partCrossProd, 1},
    {"_BigDataStatMeth_partCrossProd_block", (DL_FUNC) &_BigDataStatMeth_partCrossProd_block, 1},
    {"_BigDataStatMeth_parCrossProd", (DL_FUNC) &_BigDataStatMeth_parCrossProd, 1},
    {"_BigDataStatMeth_parCrossProd_block", (DL_FUNC) &_BigDataStatMeth_parCrossProd_block, 1},
    {"_BigDataStatMeth_partCrossProdEigen", (DL_FUNC) &_BigDataStatMeth_partCrossProdEigen, 1},
    {"_BigDataStatMeth_parxwxt", (DL_FUNC) &_BigDataStatMeth_parxwxt, 2},
    {"_BigDataStatMeth_parxtwx", (DL_FUNC) &_BigDataStatMeth_parxtwx, 2},
    {"_BigDataStatMeth_parXy", (DL_FUNC) &_BigDataStatMeth_parXy, 2},
    {"_BigDataStatMeth_bdpseudoinv", (DL_FUNC) &_BigDataStatMeth_bdpseudoinv, 1},
    {"_BigDataStatMeth_bdQR", (DL_FUNC) &_BigDataStatMeth_bdQR, 2},
    {"_BigDataStatMeth_bddtrsm", (DL_FUNC) &_BigDataStatMeth_bddtrsm, 3},
    {"_BigDataStatMeth_Create_HDF5_matrix_file", (DL_FUNC) &_BigDataStatMeth_Create_HDF5_matrix_file, 5},
    {"_BigDataStatMeth_Add_HDF5_matrix", (DL_FUNC) &_BigDataStatMeth_Add_HDF5_matrix, 5},
    {"_BigDataStatMeth_Remove_HDF5_element", (DL_FUNC) &_BigDataStatMeth_Remove_HDF5_element, 2},
    {"_BigDataStatMeth_bdSolve", (DL_FUNC) &_BigDataStatMeth_bdSolve, 2},
    {"_BigDataStatMeth_bdInvCholesky", (DL_FUNC) &_BigDataStatMeth_bdInvCholesky, 1},
    {"_BigDataStatMeth_bdSVD", (DL_FUNC) &_BigDataStatMeth_bdSVD, 5},
    {"_BigDataStatMeth_bdSVD_hdf5", (DL_FUNC) &_BigDataStatMeth_bdSVD_hdf5, 8},
    {"_BigDataStatMeth_bdSVD_lapack", (DL_FUNC) &_BigDataStatMeth_bdSVD_lapack, 3},
    {"_BigDataStatMeth_tCrossprod_Weighted", (DL_FUNC) &_BigDataStatMeth_tCrossprod_Weighted, 5},
    {"_BigDataStatMeth_parallelVectorSum", (DL_FUNC) &_BigDataStatMeth_parallelVectorSum, 1},
    {"_BigDataStatMeth_parallelpow2", (DL_FUNC) &_BigDataStatMeth_parallelpow2, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_BigDataStatMeth(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
