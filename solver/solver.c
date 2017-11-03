#include "global_header.h"

int parallel_conjugate_gradient(MODEL_STRUCT *model, int myid, int numprocs, double tol, int maxit)
{
    double **Ap, *Fp, *Vp;
    double *Rp, *Pp,*Qp;
    int IntNodes = model->Interior_nodes;
    Ap = model->interior_stiffness;
    Fp = model->interior_forcing;
    Vp = model->interior_solution;
    double GammaOld, GammaNew, tau, alpha, beta;
    int i, j, k;
    NODE_STRUCT *nodes = model->nodes;
    ELEMENT_STRUCT *element = model->elems;
    Pp = (double *) calloc(sizeof(double) , IntNodes);
    Qp = (double *) calloc(sizeof(double) , IntNodes);
    Rp = (double *) calloc(sizeof(double) , IntNodes);



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
    Ps = (double *) calloc(sizeof(double) , IBNodes);
    Qs = (double *) calloc(sizeof(double) , IBNodes);
    Rs = (double *) calloc(sizeof(double) , IBNodes);
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



    if (GammaNew < tol*tol){
        return(1);
    }
    if (myid == 0){
        printf("\n\nGamma = % 23.15e\n\n", GammaNew);
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
            printf("\n\nGamma = % 23.15e (k = %d)\n\n", GammaNew, k);
        }
        if (GammaNew < tol*tol){
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
    return(0);  // This means conjugate iterations have failed to converge
}
