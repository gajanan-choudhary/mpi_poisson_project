#include "global_header.h"

static time_t time1, time2;    /* For total simulation time calculation */

/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/
int main(int argc, char *argv[])
{
    time(&time1);

    int nmodels = 1;
    int rank, world_size;
#ifdef _MPI
    MPI_Status status;
    MPI_Init (&argc, &argv);    /* !!!!!!!!!!!! CHECK THIS !!!!!!!!!!!!!! */
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#else
    rank       = 0;
    world_size = 1;
#endif

    if (argc < 2) {
        if (rank == 0 ) {
            printf("\nError: No command line arguments supplied!\n");
            printf("Please enter the name of the model as a command line argument...\n");
            printf("Terminating program run...\n");
        }
#ifdef _MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        return(1);
    }
    else if (argc > 2) {
        if (rank == 0 ) {
            printf("\nError: Too many command line arguments supplied!\n");
            printf("Please enter only name of the model as a command line argument...\n");
            printf("Terminating program run...\n");
        }
#ifdef _MPI
        MPI_Abort(MPI_COMM_WORLD, -1);
#endif
        return(1);
    }
    else {
        if (rank == 0) {
            printf("\n*****************************************************\n");
            printf("Starting program run\n");
            printf("*****************************************************\n");
        }
    }

    // model data
    MODEL_STRUCT *model = (MODEL_STRUCT *) malloc(sizeof(MODEL_STRUCT) * nmodels);

    main_initialize(model, nmodels, argv[1], rank, world_size);

    main_run(model, nmodels, rank, world_size);

    main_finalize(model, nmodels, rank, world_size);

    time(&time2);
    if (rank == 0) {
        printf("\n*****************************************************\n");
        printf("Program ended successfully. Run time = %.3fs\n", difftime(time2,time1));
        printf("*****************************************************\n");
    }

#ifdef _MPI
    MPI_Finalize();
#endif

    return (0);
}
/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/
/***********************************************************************************************************************************************************************/




/***********************************************************************************************************************************************************************/
void main_initialize(MODEL_STRUCT *model, int nmodels, char *model_name, int myid, int numprocs)
{
    int i;
    if (myid == 0){
        printf("\n*****************************************************\n");
        printf("Main initializing\n");
    }

    // Allocating memory and initializing stuff
    if (myid == 0){
        printf("Reading input\n");
    }
    for (i=0; i<nmodels; i++){
        read_mesh(&(model), myid, numprocs, model_name);
    }
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
