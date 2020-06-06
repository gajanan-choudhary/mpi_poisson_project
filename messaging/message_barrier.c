#include "global_header.h"

/*!
 *    \brief Places a barrier in the communications 
 *     */
void message_barrier(MPI_Comm COMM){
#ifdef _MESSG
    int ierr_code = MPI_ERR_UNKNOWN;  /* the error code from an mpi call */

    fflush(stdout); /* Forcing all output to be dumped.  */

    ierr_code = MPI_Barrier(COMM);
    if (ierr_code != MPI_SUCCESS){
        message_error(ierr_code);
    }

#endif
    return;
}
