#include "global_header.h"

/******************************************************************************************************************/
double inner_product(
    double *vect1_interior,
    double *vect2_interior,
    int interiorlength
#ifdef _MPI
    ,  double *vect1_interiorboundary,
       double *vect2_interiorboundary,
       int interiorboundarylength,
       int *shared
#endif
)
{
    int i;
    double inner_product_local, inner_product_final;
    for (i=0; i<interiorlength; i++){
        inner_product_local += vect1_interior[i] * vect2_interior[i];
    }

#ifdef _MPI
    for (i=0; i<interiorboundarylength; i++){
        inner_product_local += (vect1_interiorboundary[i] * vect2_interiorboundary[i]) / shared[i];
    }
    MPI_Allreduce(&inner_product_local, &inner_product_final, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return (inner_product_final);
#else
    return (inner_product_local);
#endif
}

/******************************************************************************************************************/
#ifdef _MPI
void update(double *a /*Vector to be updated */, int numprocs, double *buf, NODE_STRUCT *nodes,
            int **neighbors, int *common)
{
    int i, j;
    MPI_Status status;
    for (i = 0; i <numprocs; i++){
        for (j=0; j< common[i]; j++){
            buf[j] = a[nodes[neighbors[i][j]].local];
            if (common[i]>0){
                // Send vector of length common[i]
                MPI_Send(&buf[0], common[i], MPI_DOUBLE, i, 10, MPI_COMM_WORLD);
            }
        }
    }

    for (i = 0; i <numprocs; i++){
        if (common[i]>0){
            // Receive vector of length common[i]
            MPI_Recv(&buf[0], common[i], MPI_DOUBLE, i, 10, MPI_COMM_WORLD, &status);
        }
        // Update local copy of the variable a
        for (j=0; j< common[i]; j++){
            a[nodes[neighbors[i][j]].local] += buf[j];
        }
    }
}
#endif

/******************************************************************************************************************/
double apply_dirichlet_bc(double x1, double x2)
{
//    return(1);
    return (x1+x2);
}

/******************************************************************************************************************/
