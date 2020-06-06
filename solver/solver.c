#include "global_header.h"

int parallel_conjugate_gradient(MODEL_STRUCT *model, int myid, int numprocs, double tol, int maxit)
{
    tol *= tol; /* Actual tolerance, since we will compare with the square of the residual norm. */

    int i, j, k;
    double **Ap, *Fp, *Vp;
    double *Rp, *Pp,*Qp;
    double GammaOld, GammaNew, tau, alpha, beta;
    int IntNodes = model->Interior_nodes;
    Ap = model->interior_stiffness;
    Fp = model->interior_forcing;
    Vp = model->interior_solution;
    NODE_STRUCT *nodes = model->nodes;
    ELEMENT_STRUCT *element = model->elems;


#ifdef _MPI
    double **As, **Bp, *Fs, *Vs;
    double *Rs, *Ps,*Qs;
    double *buf = model->buf;
    int *common = model->common;
    int **neighbors = model->neighbors;
    int *shared = model->shared;
    int IBNodes = model->Interior_boundary_nodes;
    Bp = model->interaction_stiffness;
    As = model->constrained_stiffness;
    Fs = model->constrained_forcing;
    Vs = model->constrained_solution;
#endif

    /* Allocate temporary variables. */
    Pp = calloc(sizeof(double), IntNodes);
    Qp = calloc(sizeof(double), IntNodes);
    Rp = calloc(sizeof(double), IntNodes);

#ifdef _MPI
    Ps = calloc(sizeof(double), IBNodes);
    Qs = calloc(sizeof(double), IBNodes);
    Rs = calloc(sizeof(double), IBNodes);
#endif



/********************************************************************************************************************/
#ifdef _MPI
    update(Fs, numprocs, buf, nodes, neighbors, common);
#endif
    for (i=0; i<IntNodes; i++){
        Pp[i] = Rp[i] = Fp[i];
        Vp[i] = 0.0; /* gkc 2017 */
    }
#ifdef _MPI
    for (i=0; i<IBNodes; i++){
        Ps[i] = Rs[i] = Fs[i];
        Vs[i] = 0.0; /* gkc 2017 */
    }
#endif
/********************************************************************************************************************/

/*********************************************************************/
    GammaNew = inner_product(Rp, Rp, IntNodes
#ifdef _MPI
                             , Rs, Rs, IBNodes, shared
#endif
                            );
/*********************************************************************/

    print_vector_double("Forcing before CG", Fp, IntNodes,
#ifdef _MPI
                                    Fs, IBNodes,
#endif
                                    myid, 2);
    print_vector_double("Solution before CG", Vp, IntNodes,
#ifdef _MPI
                                    Vs, IBNodes,
#endif
                                    myid, 2);

    if (GammaNew < tol){
        return(1);
    }
    if (myid == 0){
        printf("\nBefore starting CG iterations: Gamma = % 14.6e\n", GammaNew);
    }
/********************************************************************************************************************/
    for (k=0; k<maxit; k++){
        for (i=0; i<IntNodes; i++){
            Qp[i] = 0.0;
            for (j=0; j<IntNodes; j++){
                Qp[i] += Ap[i][j] * Pp[j];
            }
#ifdef _MPI
            for (j=0; j<IBNodes; j++){
                Qp[i] += Bp[i][j] * Ps[j];
            }
        }
        for (i=0; i<IBNodes; i++){
            Qs[i] = 0.0;
            for (j=0; j<IntNodes; j++){
                Qs[i] += Bp[j][i] * Pp[j]; // Not very good cache use here!
            }
            for (j=0; j<IBNodes; j++){
                Qs[i] += As[i][j] * Ps[j];
            }
#endif
        }


#ifdef _MPI
        update(Qs, numprocs, buf, nodes, neighbors, common);
#endif



/*********************************************************************/
        tau = inner_product(Pp, Qp, IntNodes

#ifdef _MPI
                           , Ps, Qs, IBNodes, shared
#endif
                           );
/*********************************************************************/


        alpha = GammaNew/tau;       // This is ( (R^T)*R ) / ( (P^T)*K*P )
        for (i=0; i<IntNodes; i++){
            Vp[i] += alpha * Pp[i];  // x(k+1) = x(k) + alpha * P(k)
            Rp[i] -= alpha * Qp[i];  // R(k+1) = R(k) - alpha * K*P(k)
        }
#ifdef _MPI
        for (i=0; i<IBNodes; i++){
            Vs[i] += alpha * Ps[i];
            Rs[i] -= alpha * Qs[i];
        }
#endif
        GammaOld = GammaNew;


        if (myid == 0){
            printf("\nIteration % 4i, Gamma = % 14.6e", k, GammaNew);
        }
        if (GammaNew < tol){
            print_vector_double("Solution", Vp, IntNodes,
#ifdef _MPI
                                Vs, IBNodes,
#endif
                                myid, 2);
            return(1);
        }

/*********************************************************************/
        GammaNew = inner_product(Rp, Rp, IntNodes
#ifdef _MPI
                                 , Rs, Rs, IBNodes, shared
#endif
                                );
/*********************************************************************/



        beta = GammaNew/GammaOld;
        for (i=0; i<IntNodes; i++){
            Pp[i] = Rp[i]+beta*Pp[i];
        }
#ifdef _MPI
        for (i=0; i<IBNodes; i++){
            Ps[i] = Rs[i]+beta*Ps[i];
        }
#endif
    }
    
    free(Pp); free(Qp); free(Rp);
#ifdef _MPI
    free(Ps); free(Qs); free(Rs);
#endif

    return(0);  // This means conjugate iterations have failed to converge
}
