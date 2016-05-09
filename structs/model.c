#include "global_header.h"

/***********************************************************************************************************************************************************************/
void model_alloc_matrices(MODEL_STRUCT *model)
{
    int i, j;
    model->interior_stiffness = (double **) calloc(sizeof(double *) , model->Interior_nodes);
    for(i=0; i<model->Interior_nodes; i++){
        model->interior_stiffness[i] = (double *) calloc(sizeof(double) , model->Interior_nodes);
    }
    model->interior_forcing  = (double *) calloc(sizeof(double) , model->Interior_nodes);
    model->interior_solution = (double *) calloc(sizeof(double) , model->Interior_nodes);

#ifdef _MPI
    model->interaction_stiffness = (double **) calloc(sizeof(double *) , model->Interior_nodes);
    for(i=0; i<model->Interior_nodes; i++){
        model->interaction_stiffness[i] = (double *) calloc(sizeof(double) , model->Interior_boundary_nodes);
    }

    model->constrained_stiffness = (double **) calloc(sizeof(double *) , model->Interior_boundary_nodes);
    for(i=0; i<model->Interior_boundary_nodes; i++){
        model->constrained_stiffness[i] = (double *) calloc(sizeof(double) , model->Interior_boundary_nodes);
    }

    model->constrained_forcing  = (double *) calloc(sizeof(double) , model->Interior_boundary_nodes);
    model->constrained_solution = (double *) calloc(sizeof(double) , model->Interior_boundary_nodes);
#endif
}

/***********************************************************************************************************************************************************************/
void model_build_system(MODEL_STRUCT *model)
{
    int i, j, k, ii, jj;
    double stiffji;
    double grad0[NDONTRI], grad1[NDONTRI];
    double locxy[NDONTRI][ndim];

    // Assigning pointers for convenience
    NODE_STRUCT *nodes;              nodes = model->nodes;
    ELEMENT_STRUCT *elems;           elems = model->elems;
    double **IS;                     IS    = model->interior_stiffness;
    double *IF;                      IF    = model->interior_forcing;
#ifdef _MPI
    double **IAS;                    IAS   = model->interaction_stiffness;
    double **CS;                     CS    = model->constrained_stiffness;
    double *CF;                      CF    = model->constrained_forcing;
#endif

    for (k = 0; k<model->nelems; k++){
        double area = elems[k].area;
        double twoA = 2*area;
        for (i=0; i<NDONTRI; i++){
            for (j=0; j<ndim; j++){
                locxy[i][j] = nodes[elems[k].vertex[i]].xy[j];
            }
        }
        
        grad0[0] = (locxy[1][1] - locxy[2][1])/twoA;       // (y1 - y2) / 2    // Factor " / Area " redundant
        grad0[1] = (locxy[2][1] - locxy[0][1])/twoA;       // (y2 - y0) / 2
        grad0[2] = (locxy[0][1] - locxy[1][1])/twoA;       // (y0 - y1) / 2

        grad1[0] = (locxy[1][0] - locxy[2][0])/twoA;       // (x1 - x2) / 2
        grad1[1] = (locxy[2][0] - locxy[0][0])/twoA;       // (x2 - x0) / 2
        grad1[2] = (locxy[0][0] - locxy[1][0])/twoA;       // (x0 - x1) / 2
        
        for (j=0; j<NDONTRI; j++){
            jj = elems[k].vertex[j];
//#ifdef _MPI
            if (nodes[jj].type == 1){
//#endif
                for (i=0; i<NDONTRI; i++){
                    ii = elems[k].vertex[i];
                    stiffji = area*(grad0[i]*grad0[j] + grad1[i]*grad1[j]);
                    if (nodes[ii].type == 1){
                        IS[nodes[jj].local][nodes[ii].local] += stiffji;
                    }
#ifdef _MPI
                    else if (nodes[ii].type > 1){
                        IAS[nodes[jj].local][nodes[ii].local] += stiffji;
                    }
#endif
                    else if (nodes[ii].type == 0){
                        IF[nodes[jj].local] -= apply_dirichlet_bc(nodes[ii].xy[0], nodes[ii].xy[1])*stiffji;
                    }
                }
            }
#ifdef _MPI
            else if (nodes[jj].type > 1){
                for (i=0; i<NDONTRI; i++){
                    ii = elems[k].vertex[i];
                    stiffji = area*(grad0[i]*grad0[j] + grad1[i]*grad1[j]);
                    if (nodes[ii].type > 1){
                        CS[nodes[jj].local][nodes[ii].local] += stiffji;
                    }
                    else if (nodes[ii].type == 0){
                        CF[nodes[jj].local] -= apply_dirichlet_bc(nodes[ii].xy[0], nodes[ii].xy[1])*stiffji;
                    }
                }
            }
#endif
        }
    }
}

/***********************************************************************************************************************************************************************/
void model_print(MODEL_STRUCT *model, int nmodels, int myid){
    int i, j, k;
    for (i=0;i<nmodels;i++){
        printf("\nProcessor %i: Model[%i] data below:\n", myid, i);
        printf("\tNumber of nodes    : %i\n\tNumber of elements : %i\n", model[i].nnodes, model[i].nelems);

        printf("\tProcessor %i: Node Data\n", myid);
        printf("\t\t|   Node    |   Type    |     X coordinate     |     Y coordinate     |\n");
        for (j=0;j<model[i].nnodes;j++){
            printf("\t\t| %10i| %10i| %.15e| %.15e|\n", j, model[i].nodes[j].type,
                                    model[i].nodes[j].xy[0], model[i].nodes[j].xy[1]);
        }

        printf("\tProcesssor %i:Element Data\n", myid);
        printf("\t\t|  Element  | Vertex 0  | Vertex 1  | Vertex 2  |         Area         | Processor |\n");
        for (j=0;j<model[i].nelems;j++){
            printf("\t\t| %10i| %10i| %10i| %10i| %.15e| %10i|\n", j,
                model[i].elems[j].vertex[0], model[i].elems[j].vertex[1], model[i].elems[j].vertex[2],
                                                  model[i].elems[j].area, model[i].elems[j].procnum);
        }

        printf("\tProcessor %i: interior_stiffness (matrix):\n", myid);
        for (j=0;j<model[i].Interior_nodes;j++){
            printf("Row %10i:", j);
            for (k=0;k<model[i].Interior_nodes;k++){
                printf("\t%.4e", model[i].interior_stiffness[j][k]);
            }
            printf("\n");
        }
#ifdef _MPI
        printf("\tProcessor %i: interaction_stiffness (matrix):\n", myid);
        for (j=0;j<model[i].Interior_nodes;j++){
            printf("Row %10i:", j);
            for (k=0;k<model[i].Interior_boundary_nodes;k++){
                printf("\t%.4e", model[i].interaction_stiffness[j][k]);
            }
            printf("\n");
        }

        printf("\tProcessor %i: constrained_stiffness (matrix):\n", myid);
        for (j=0;j<model[i].Interior_boundary_nodes;j++){
            printf("Row %10i:", j);
            for (k=0;k<model[i].Interior_boundary_nodes;k++){
                printf("\t%.4e", model[i].constrained_stiffness[j][k]);
            }
            printf("\n");
        }
#endif
        printf("\tProcessor %i: interior_forcing (vector):\n", myid);
        for (j=0;j<model[i].Interior_nodes; j++){
            printf("Row %10i:\t%.4e\n", j, model[i].interior_forcing[j]);
        }
#ifdef _MPI
        printf("\tProcessor %i: constrained_forcing (vector):\n", myid);
        for (j=0;j<model[i].Interior_boundary_nodes; j++){
            printf("Row %10i:\t%.4e\n", j, model[i].constrained_forcing[j]);
        }
#endif
    }
}

/***********************************************************************************************************************************************************************/
void model_free(MODEL_STRUCT *model, int nmodels, int myid, int numprocs){
    int i, j;
    for (i=0;i<nmodels;i++){
        free(model[i].nodes);
        free(model[i].elems);
        for (j=0;j<model[i].Interior_nodes;j++){
            free(model[i].interior_stiffness[j]);
#ifdef _MPI
            free(model[i].interaction_stiffness[j]);
        }
        for (j=0;j<model[i].Interior_boundary_nodes;j++){
            free(model[i].constrained_stiffness[j]);
#endif
        }
        free(model[i].interior_stiffness);
        free(model[i].interior_forcing);
        free(model[i].interior_solution);

#ifdef _MPI
        for (j=0;j<numprocs;j++){
            if (model[i].common[j]>0){
                free(model[i].neighbors[j]);
            }
        }
        free(model[i].neighbors);
        free(model[i].common);
        free(model[i].shared);
        free(model[i].buf);
        free(model[i].interaction_stiffness);
        free(model[i].constrained_stiffness);
        free(model[i].constrained_forcing);
        free(model[i].constrained_solution);
#endif
    }
    free(model);
}

/***********************************************************************************************************************************************************************/
