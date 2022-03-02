#include "include/BigDataStatMeth.h"
#ifdef _OPENMP
#include <pthread.h>
#endif
#include <errno.h>     // errno
#include <ctype.h>     // isspace

static int  DTthreads = -1;   // Never read directly hence static; use getDTthreads(n, /*throttle=*/0|1). -1 so we know for sure initDTthreads() ran and set it >= 1.
static int  DTthrottle = -1;  // Thread 1 is assigned DTthrottle iterations before a 2nd thread is utilized; #4484.
static bool RestoreAfterFork = true;  // see #2885 in v1.12.0
static int getIntEnv(const char *name, int def);

static int getIntEnv(const char *name, int def)
{
    const char *val = getenv(name);
    if (val==NULL) return def;
    size_t nchar = strlen(val);
    if (nchar==0) return def;
    char *end;
    errno = 0;
    long int ans = strtol(val, &end, 10);  // ignores leading whitespace. If it fully consumed the string, *end=='\0' and isspace('\0')==false
    while (isspace(*end)) end++;  // ignore trailing whitespace
    if (errno || (size_t)(end-val)!=nchar || ans<1 || ans>INT_MAX) {
        Rcpp::warning(("Ignoring invalid %s==\"%s\". Not an integer >= 1. Please remove any characters that are not a digit [0-9]."), name, val);
        return def;
        
    }
    return (int)ans;
}

static inline int imin(int a, int b) { return a < b ? a : b; }
static inline int imax(int a, int b) { return a > b ? a : b; }

void initDTthreads() {
    // called at package startup from init.c
    // also called by setDTthreads(threads=NULL) (default) to reread environment variables; see setDTthreads below
    // No verbosity here in this setter. Verbosity is in getDTthreads(verbose=TRUE)
    int ans = getIntEnv("R_DATATABLE_NUM_THREADS", INT_MIN);
    if (ans>=1) {
        ans = imin(ans, omp_get_num_procs());  // num_procs is a hard limit; user cannot achieve more. ifndef _OPENMP then myomp.h defines this to be 1
    } else {
        // Only when R_DATATABLE_NUM_THREADS is unset (or <=0) do we use PROCS_PERCENT; #4514
        int perc = getIntEnv("R_DATATABLE_NUM_PROCS_PERCENT", 50); // use "NUM_PROCS" to use the same name as the OpenMP function this uses
        // 50% of logical CPUs by default; half of 8 is 4 on laptop with 4 cores. Leaves plenty of room for other processes: #3395 & #3298
        if (perc<=1 || perc>100) {
            Rcpp::warning(("Ignoring invalid R_DATATABLE_NUM_PROCS_PERCENT==%d. If used it must be an integer between 2 and 100. Default is 50. See ?setDTtheads."), perc);
            // not allowing 1 is to catch attempts to use 1 or 1.0 to represent 100%.
            perc = 50;
        }
        ans = imax(omp_get_num_procs()*perc/100, 1); // imax for when formula would result in 0.
    }
    ans = imin(ans, omp_get_thread_limit());  // honors OMP_THREAD_LIMIT when OpenMP started; e.g. CRAN sets this to 2. Often INT_MAX meaning unlimited/unset
    ans = imin(ans, omp_get_max_threads());   // honors OMP_NUM_THREADS when OpenMP started, plus reflects any omp_set_* calls made since
    // max_threads() -vs- num_procs(): https://software.intel.com/en-us/forums/intel-visual-fortran-compiler-for-windows/topic/302866
    ans = imin(ans, getIntEnv("OMP_THREAD_LIMIT", INT_MAX));  // user might expect `Sys.setenv(OMP_THREAD_LIMIT=2);setDTthreads()` to work. Satisfy this
    ans = imin(ans, getIntEnv("OMP_NUM_THREADS", INT_MAX));   //   expectation by reading them again now. OpenMP just reads them on startup (quite reasonably)
    ans = imax(ans, 1);  // just in case omp_get_* returned <=0 for any reason, or the env variables above are set <=0
    DTthreads = ans;
    DTthrottle = imax(1, getIntEnv("R_DATATABLE_THROTTLE", 1024)); // 2nd thread is used only when n>1024, 3rd thread when n>2048, etc
}

int getDTthreads(const int64_t n, const bool throttle) {
    
    initDTthreads();
    
    // throttle==true  : a number of iterations per thread (DTthrottle) is applied before a second thread is utilized
    // throttle==false : parallel region is already pre-chunked such as in fread; e.g. two batches intended for two threads
    if (n<1) return 1; // 0 or negative could be deliberate in calling code for edge cases where loop is not intended to run at all
    int64_t ans = throttle ? 1+(n-1)/DTthrottle :  // 1 thread for n<=1024, 2 thread for n<=2048, etc
        n;                    // don't use 20 threads for just one or two batches
    return ans>=DTthreads ? DTthreads : (int)ans;  // apply limit in static local DTthreads saved there by initDTthreads() and setDTthreads()
}

static const char *mygetenv(const char *name, const char *unset) {
    const char *ans = getenv(name);
    return (ans==NULL || ans[0]=='\0') ? unset : ans;
}

SEXP getDTthreads_R(SEXP verbose) {
    if(!IS_TRUE_OR_FALSE(verbose))
        Rf_error(("%s must be TRUE or FALSE"), "verbose");
        if (LOGICAL(verbose)[0]) {
#ifndef _OPENMP
        Rprintf(("This installation of BigDataStatMeth has not been compiled with OpenMP support.\n"));
#else
        Rprintf(("  OpenMP version (_OPENMP)       %d\n"), _OPENMP); // user can use Google to map 201511 to 4.5; it's odd that OpenMP API does not provide 4.5
#endif
        Rprintf(("  omp_get_num_procs()            %d\n"), omp_get_num_procs());
        Rprintf(("  R_DATATABLE_NUM_PROCS_PERCENT  %s\n"), mygetenv("R_DATATABLE_NUM_PROCS_PERCENT", "unset (default 50)"));
        Rprintf(("  R_DATATABLE_NUM_THREADS        %s\n"), mygetenv("R_DATATABLE_NUM_THREADS", "unset"));
        Rprintf(("  R_DATATABLE_THROTTLE           %s\n"), mygetenv("R_DATATABLE_THROTTLE", "unset (default 1024)"));
        Rprintf(("  omp_get_thread_limit()         %d\n"), omp_get_thread_limit());
        Rprintf(("  omp_get_max_threads()          %d\n"), omp_get_max_threads());
        Rprintf(("  OMP_THREAD_LIMIT               %s\n"), mygetenv("OMP_THREAD_LIMIT", "unset"));  // CRAN sets to 2
        Rprintf(("  OMP_NUM_THREADS                %s\n"), mygetenv("OMP_NUM_THREADS", "unset"));
        Rprintf(("  RestoreAfterFork               %s\n"), RestoreAfterFork ? "true" : "false");
        Rprintf(("  BigDataStatMeth is using %d threads with throttle==%d.\n"), getDTthreads(INT_MAX, false), DTthrottle);
    }
    return Rf_ScalarInteger(getDTthreads(INT_MAX, false));
}

static int pre_fork_DTthreads = 0;

void when_fork() {
    pre_fork_DTthreads = DTthreads;
    DTthreads = 1;
}

void after_fork() {
    if (RestoreAfterFork) DTthreads = pre_fork_DTthreads;
}

void avoid_openmp_hang_within_fork() {
    // Called once on loading BigDataStatMeth from init.c
#ifdef _OPENMP
    pthread_atfork(&when_fork, &after_fork, NULL);
#endif
}

