#include "global_header.h"

static int DEBUG = OFF;

void print_vector_double(
    const char *print_string,
    double *interior_v,
    int interior_nnodes,
#ifdef _MPI
    double *constrained_v,
    int interior_boundary_nnodes,
#endif
    int myid,
    int type /* can be 1 or 2 */
){
    int j, irow=0;
    if (DEBUG==ON && type==1){
        printf(COLOR_RED "Processor %i: interior %s:" COLOR_RESET "\n", myid, print_string);
        for (j=0;j<interior_nnodes; j++){
            printf(COLOR_RED "Row % 10i:\t% .4e" COLOR_RESET "\n", j, interior_v[j]);
        }
#ifdef _MPI
        printf(COLOR_GREEN "Processor %i: constrained %s:" COLOR_RESET "\n", myid, print_string);
        for (j=0;j<interior_boundary_nnodes; j++){
            printf(COLOR_GREEN "Row % 10i:\t% .4e" COLOR_RESET "\n", j, constrained_v[j]);
        }
#endif
    }
    else if (DEBUG==OFF && type==1){
        printf("Processor %i: %s:\n", myid, print_string);
        for (j=0;j<interior_nnodes; j++){
            printf(COLOR_RED "Row % 10i:\t% .4e" COLOR_RESET "\n", irow++, interior_v[j]);
        }
#ifdef _MPI
        for (j=0;j<interior_boundary_nnodes; j++){
            printf(COLOR_GREEN "Row % 10i:\t% .4e" COLOR_RESET "\n", irow++, constrained_v[j]);
        }
#endif
    }
    else if (type==2){
#ifdef _MESSG
        message_barrier(MPI_COMM_WORLD);
        int proc_count;
        for(proc_count=0;proc_count<4;proc_count++){
            if(proc_count==myid){
                printf("\n" COLOR_RED "*********** Processor %d printing %s below:" COLOR_RESET "\n", myid, print_string);
#else
                printf("\t%s:\n", myid, print_string);
#endif
                for (j=0;j<interior_nnodes; j++){
                    printf("Row % 10i:\t% .4e\n", irow++, interior_v[j]);
                }
#ifdef _MPI
                for (j=0;j<interior_boundary_nnodes; j++){
                    printf("Row % 10i:\t% .4e\n", irow++, constrained_v[j]);
                }
#endif
#ifdef _MESSG
            }
            message_barrier(MPI_COMM_WORLD);
        }
        //exit(-1);
#endif
    }
    printf("\n");
}

void print_stiffness_matrix(
    double **interior_stiffness,
    int interior_nnodes,
#ifdef _MPI
    double **interaction_stiffness,
    double **constrained_stiffness,
    int interior_boundary_nnodes,
#endif
    int myid
){
    int j, k;
    if (DEBUG==ON){
        printf(COLOR_RED "Processor %i: interior_stiffness (matrix):" COLOR_RESET "\n", myid);
        for (j=0;j<interior_nnodes;j++){
            printf(COLOR_RED "Row % 10i:", j);
            for (k=0;k<interior_nnodes;k++){
                printf("\t% .4e", interior_stiffness[j][k]);
            }
            printf(COLOR_RESET "\n");
        }
#ifdef _MPI
        printf(COLOR_BLUE "Processor %i: interaction_stiffness (matrix):" COLOR_RESET "\n", myid);
        for (j=0;j<interior_nnodes;j++){
            printf(COLOR_BLUE "Row % 10i:", j);
            for (k=0;k<interior_boundary_nnodes;k++){
                printf("\t% .4e", interaction_stiffness[j][k]);
            }
            printf(COLOR_RESET "\n");
        }
    
        printf(COLOR_GREEN "Processor %i: constrained_stiffness (matrix):" COLOR_RESET "\n", myid);
        for (j=0;j<interior_boundary_nnodes;j++){
            printf(COLOR_GREEN "Row % 10i:", j);
            for (k=0;k<interior_boundary_nnodes;k++){
                printf("\t% .4e", constrained_stiffness[j][k]);
            }
            printf(COLOR_RESET "\n");
        }
#endif
    }
    else{
        printf("Processor %i: Full stiffness matrix:\n", myid);
        int irow=0;
        for (j=0;j<interior_nnodes;j++){
            printf(COLOR_RED "Row % 10i:" COLOR_RESET, irow++);
            printf(COLOR_RED);
            for (k=0;k<interior_nnodes;k++){
                printf("\t% .4e", interior_stiffness[j][k]);
            }
            printf(COLOR_RESET);
#ifdef _MPI
            printf(COLOR_BLUE);
            for (k=0;k<interior_boundary_nnodes;k++){
                printf("\t% .4e", interaction_stiffness[j][k]);
            }
            printf(COLOR_RESET);
#endif
            printf("\n");
        }
#ifdef _MPI
        for (j=0;j<interior_boundary_nnodes;j++){
            printf(COLOR_GREEN "Row % 10i:" COLOR_RESET, irow++);
            printf(COLOR_BLUE);
            for (k=0;k<interior_nnodes;k++){
                printf("\t% .4e", interaction_stiffness[k][j]); /* Note the transpose, [k][j] here. */
            }
            printf(COLOR_RESET COLOR_GREEN);
            for (k=0;k<interior_boundary_nnodes;k++){
                printf("\t% .4e", constrained_stiffness[j][k]);
            }
            printf(COLOR_RESET "\n");
        }
#endif
    }
    printf("\n");
}
