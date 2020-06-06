#include "global_header.h"

static int DEBUG = OFF;

/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/
int unit_tests(
#ifdef _MPI
    int **neighbors,
    int *common,
    int npes,
    int myid
#endif
)
{
    int nfailures; /* Number of tests that failed. */

#ifdef _MPI
    /* Test if vectors are updated correctly in parallel. Suppose node 1 is shared with PE0 and PE1. */
    /* If in parallel, PE0 has v={v1, v2}, and PE1 has v={v3, v4}, */
    /* then in serial, the corresponding correct vector would be v={v1, v2+v3, v4}^T, */
    /* and the vectors after updates should be PE0: v={v1, v2+v3}, and PE1: v={v2+v3, v4}^T. */
    
    
#endif
    return (nfailures);
}
/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/




/***********************************************************************************************************************************************************************/
void test_update_double(MODEL_STRUCT *model, int myid, int numprocs)
{
    int i;
    //double *v=model->Fs;
    if (myid == 0){
        printf("\n*****************************************************\n");
        printf("Main initializing\n");
    }

    // Allocating memory and initializing stuff
    if (myid == 0){
        printf("Reading input\n");
    }
    //update(v /*Vector to be updated */, numprocs, model->buf, model->nodes, model->neighbors, model->common);
}

/***********************************************************************************************************************************************************************/
void main_run(MODEL_STRUCT *model, int nmodels, int myid, int numprocs)
{
    int i, j, ierr;
    if (myid == 0){
        printf("\n*****************************************************\n");
        printf("Main running\n");
    }

    for (i=0; i<nmodels; i++){
        model_alloc_matrices(&(model[i]));
        model_build_system(&(model[i]));
        // Printing initial data
        model_print(&(model[i]), 1, myid);

        if (myid == 0){
            printf("\n*******************************************\n");
            printf("Solving the system using conjugate gradient iterations\n");
            printf("*******************************************\n\n");
        }
        ierr = parallel_conjugate_gradient(&(model[i]), myid, numprocs, 1e-4, 20);

#ifdef _MESSG
        message_barrier(MPI_COMM_WORLD);
        int proc_count;
        for(proc_count=0;proc_count<4;proc_count++){
            if(proc_count==myid){
                printf("\n" COLOR_RED "*********** Processor %d printing below:" COLOR_RESET, myid);
#endif
//                for (j=0; j<1000000; j++){ /* Do Nothing! */}
                printf("\nSolution:\n");
                for (j=0; j<model[i].nnodes; j++){
                    if (model[i].nodes[j].type == 1){
                        printf(COLOR_RED "Proc %i" COLOR_RESET ": %f %f %f\n", myid, model[i].nodes[j].xy[0],
                                              model[i].nodes[j].xy[1],
                                              model[i].interior_solution[model[i].nodes[j].local]);
                    }
#ifdef _MPI
                    else if (model[i].nodes[j].type > 1){
                        printf(COLOR_GREEN "Proc %i" COLOR_RESET ": %f %f %f\n", myid, model[i].nodes[j].xy[0],
                                              model[i].nodes[j].xy[1],
                                              model[i].constrained_solution[model[i].nodes[j].local]);
                    }
#endif
                }
#ifdef _MESSG
            }
            message_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }

#ifdef _MESSG
        message_barrier(MPI_COMM_WORLD);
        int proc_count;
        for(proc_count=0;proc_count<4;proc_count++){
            if(proc_count==myid){
                printf("\n" COLOR_RED "*********** Processor %d printing below:" COLOR_RESET, myid);
#endif
                if (ierr == 0){
                    printf("\n    Warning! CG iterations failed to converge!", myid);
                }
                else {
                    printf("\n    CG iterations converged.");
                }
#ifdef _MESSG
            }
            message_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
}

/***********************************************************************************************************************************************************************/
void main_finalize(MODEL_STRUCT *model, int nmodels, int myid, int numprocs)
{
//////////////////////////////////////////
#ifdef _MPI
    MPI_Barrier(MPI_COMM_WORLD);
#endif
//////////////////////////////////////////

    if (myid == 0){
        printf("\n\n*****************************************************\n");
        printf("Main Finalizing\n");
    }

    // Freeing all memory
    model_free(model, nmodels, myid, numprocs);

    if (myid == 0){
        printf("All data freed.\n");
    }

}
