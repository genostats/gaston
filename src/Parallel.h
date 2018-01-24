
#define GASTON_use_RcppParallel false

#if GASTON_use_RcppParallel
#include <RcppParallel.h>
#define Parallel RcppParallel
#else
#include "PoorMansParallel.h"
#define Parallel PoorMansParallel
#endif
