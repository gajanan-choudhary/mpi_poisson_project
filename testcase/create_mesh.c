#include <stdio.h>
#include <stdlib.h>
#define lengthX 1.0
#define lengthY 1.0

#define NDONTRI 3
#define NDIM 2
#define MAXLINE 300

/******************************************************************************************************************/
typedef struct{
    double xy[NDIM];
    int global_dof;  // Global node numbers, starting from 0 each
    int local_node;
    int type;       // Interior = 0, Exterior Boundary < 0, Interior Boundary > 0
} NODE_STRUCT;
/******************************************************************************************************************/
typedef struct {
    int vertex[NDONTRI];
    int procnum;
} ELEMENT_STRUCT;

/******************************************************************************************************************/
void write_mesh(NODE_STRUCT *nodes, ELEMENT_STRUCT *elems, int nnodes, int nelems, char *output_mesh_filename);

/******************************************************************************************************************/
/******************************************************************************************************************/
/******************************************************************************************************************/
int main(int argc, char *argv[])
{
    int nelemsX, nelemsY, nnodesX, nnodesY, numprocs;
    double deltaX, deltaY;
    if (argc != 4){
        printf("\nError. Need exactly 3 command line arguments\n");
        printf("\t1. Mesh file name\n\t2. #Elems in x direction\n\t3. #Elems in y direction\n\n");
        return (1);
    }
    else{
        nelemsX  = atoi(argv[2]);
        nelemsY  = atoi(argv[3]);
        if (nelemsX <1 || nelemsY<1){
            printf("\nError. Need non-negative command line arguments.\n\n");
            return(1);
        }
        nnodesX = nelemsX/2+1;
        nnodesY = nelemsY/2+1;
        deltaX = lengthX/(nelemsX/2);
        deltaY = lengthY/(nelemsY/2);
    }

    /************************************************************************************************/
    int i, j, k;
    int nnodes = nnodesX*nnodesY;
    int nelems = (nnodesX-1)*(nnodesY-1)*2;
    NODE_STRUCT *nodes;
    ELEMENT_STRUCT *elems;
    /************************************************************************************************/
    // Reading initial data
    nodes = (NODE_STRUCT *)    malloc(nnodes*sizeof(NODE_STRUCT));     // Allocating space for nodes
    elems = (ELEMENT_STRUCT *) malloc(nelems*sizeof(ELEMENT_STRUCT));  // Allocating space for elements

    for (i=0;i<nnodesY;i++){
        int type = 0;
        if (i==0 || i==nnodesY-1)
            type = -1;
        for (j=0;j<nnodesX;j++){
            nodes[i*nnodesX+j].xy[0] = j*deltaX;
            nodes[i*nnodesX+j].xy[1] = i*deltaY;
            nodes[i*nnodesX+j].type = type;
        }
    }

    for (i=0;i<nnodesY-1;i++){
        for (j=0;j<nnodesX-1;j++){

            int tempnodes[4];
            tempnodes[0] = i*nnodesX + j;
            tempnodes[1] = tempnodes[0]+1;
            tempnodes[2] = tempnodes[0]+nnodesX;
            tempnodes[3] = tempnodes[2]+1;

            int e = i*nelemsX+2*j;
            elems[e].vertex[0] = tempnodes[0];
            elems[e].vertex[1] = tempnodes[1];
            elems[e].vertex[2] = tempnodes[2];
//            elems[e].procnum   = i;
            elems[e].procnum   = i*nelemsX/2+j;

            e++;
            elems[e].vertex[0] = tempnodes[2];
            elems[e].vertex[1] = tempnodes[1];
            elems[e].vertex[2] = tempnodes[3];
//            elems[e].procnum   = i;
            elems[e].procnum   = i*nelemsX/2+j;
        }
    }
    /************************************************************************************************/
    // Writing partitioned mesh data
    printf("\nWriting partitioned mesh data\n");
    write_mesh(nodes, elems, nnodes, nelems, argv[1]);

    /************************************************************************************************/
    // Freeing all data
    printf("\nFreeing all data\n");
    free(nodes);
    free(elems);
    return(0);
}

void write_mesh(NODE_STRUCT *nodes, ELEMENT_STRUCT *elems, int nnodes, int nelems, char *model_name)
{
    int i;
    char output_mesh_filename[MAXLINE];
    FILE *outfile;

    sprintf(output_mesh_filename,"%s.2dm", model_name);
    outfile = fopen(output_mesh_filename,"w");

    fprintf(outfile, "%10i %10i\n", nnodes, nelems);
    for(i=0;i<nnodes;i++){
        fprintf(outfile, "%.15e %.15e %10i\n", nodes[i].xy[0], nodes[i].xy[1], nodes[i].type);
    }
    for(i=0;i<nelems;i++){
        fprintf(outfile, "%10i %10i %10i %10i\n", elems[i].vertex[0], elems[i].vertex[1], elems[i].vertex[2], elems[i].procnum);
    }
    fclose(outfile);
}
