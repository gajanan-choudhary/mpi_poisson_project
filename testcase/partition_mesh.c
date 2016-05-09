#include <stdio.h>
#include <stdlib.h>

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
void read_mesh(NODE_STRUCT *nodes, ELEMENT_STRUCT *elems, int nnodes, int nelems, int numprocs,
               FILE *inputmeshfile, int *counter, int *PartSize, int **PartN, int **PartE);
void partition(NODE_STRUCT *nodes, ELEMENT_STRUCT *elems, int nnodes, int nelems, int numprocs,
               char *input_mesh_filename, int *counter, int *common, int *PartSize, int **PartN, int **PartE);

/******************************************************************************************************************/
/******************************************************************************************************************/
/******************************************************************************************************************/
int main(int argc, char *argv[])
{

    if (argc != 3){
        printf("\nError. Need exactly 2 command line arguments\n");
        printf("\t1. The name of the mesh file to be partitioned\n\t2. Number of processors\n\n");
        return (1);
    }
    /************************************************************************************************/
    char input_mesh_filename[MAXLINE];
    FILE *inputmeshfile;                                 // Single input mesh file
    /************************************************************************************************/
    int i, j, k;
    int nnodes, nelems;
    int *PartSize;
    int **PartE;
    int **PartN;
    int *counter;
    int *common;
    NODE_STRUCT *nodes;
    ELEMENT_STRUCT *elems;

    int numprocs = atoi(argv[2]);

    sprintf(input_mesh_filename,"%s.2dm", argv[1]);
    inputmeshfile = fopen(input_mesh_filename, "r");

    fscanf(inputmeshfile,"%d%d", &nnodes, &nelems);        // Number of nodes and elements
    // Sanity check
    if(nelems < numprocs){
        fprintf(stderr, "Error! Number of elements is less than the number of processors\n"
                        "Error encountered in file %s, at line %d\n\n", __FILE__, __LINE__);
        return (1);
    }
    
    counter  = (int *)  malloc(sizeof(int  ) * numprocs); // Recalculated for each processor
    common   = (int *)  malloc(sizeof(int  ) * numprocs); // Recalculated for each processor
    PartSize = (int *)  malloc(sizeof(int  ) * numprocs); // Calculated just once
    PartE    = (int **) malloc(sizeof(int *) * numprocs); // Calculated just once
    PartN    = (int **) malloc(sizeof(int *) * numprocs); // Calculated just once
    
    /************************************************************************************************/
    // Reading initial data
    nodes = (NODE_STRUCT *)    malloc(nnodes*sizeof(NODE_STRUCT));     // Allocating space for nodes
    elems = (ELEMENT_STRUCT *) malloc(nelems*sizeof(ELEMENT_STRUCT));  // Allocating space for elements

    for(i=0;i<numprocs;i++){
        PartN[i] = (int *) malloc(sizeof(int) * nnodes);
    }

    for(i=0;i<numprocs;i++){
        PartE[i] = (int *) malloc(sizeof(int) * nelems);
        PartSize[i] = 0;
    }



    /************************************************************************************************/
    // Reading and interpreting mesh data from input file <argv[1].2dm>
    printf("\nReading mesh data from file %s.2dm\n",argv[1]);
    read_mesh(nodes, elems, nnodes, nelems, numprocs, inputmeshfile, counter, PartSize, PartN, PartE);

    printf("PartSize is given below \n");
    for (i=0; i<numprocs; i++){
        printf("\tProcessor %i : %5i\n", i, PartSize[i]);
    }
    printf("PartN is given below \n");
    for (i=0; i<numprocs; i++){
        for (j=0; j<nnodes; j++){
            printf("%5i", PartN[i][j]);
        }
        printf("\n");
    }

    printf("PartE is given below \n");
    for (i=0; i<numprocs; i++){
        for (j=0; j<nelems; j++){
            printf("%5i", PartE[i][j]);
        }
        printf("\n");
    }

    /************************************************************************************************/
    // Writing partitioned mesh data to <numprocs> files
    printf("\nWriting partitioned mesh data\n");
    partition(nodes, elems, nnodes, nelems, numprocs, argv[1], counter, common, PartSize, PartN, PartE);

    /************************************************************************************************/
    // Freeing all data
    printf("\nFreeing all data\n");
    for(i=0;i<numprocs;i++){
        free(PartN[i]);
        free(PartE[i]);
    }
    free(counter);
    free(common);
    free(PartSize);
    free(PartE);
    free(PartN);
    free(nodes);
    free(elems);
    return(0);
}


void read_mesh(NODE_STRUCT *nodes, ELEMENT_STRUCT *elems, int nnodes, int nelems, int numprocs,
               FILE *inputmeshfile, int *counter, int *PartSize, int **PartN, int **PartE)
{
    int i, j, k, currentvertex;                                  // Counters
    int procnum;
    /************************************************************************************************/
    // Reading node data
    for(i=0;i<nnodes;i++){
        fscanf(inputmeshfile, "%lf%lf%d",&(nodes[i].xy[0]), &(nodes[i].xy[1]), &(nodes[i].global_dof));
        nodes[i].type = 0;
    }
    printf("Finished reading nodes\n");
    /************************************************************************************************/
    // Reading element data
    for(i=0;i<nelems;i++){
        fscanf(inputmeshfile, "%d%d%d%d", &(elems[i].vertex[0]), &(elems[i].vertex[1]),
                                          &(elems[i].vertex[2]), &(elems[i].procnum));
        if (elems[i].procnum>=numprocs){
            fprintf(stderr, "Error! Processor number exceeds supplied number of processors\n"
                            "Error encountered in file %s, at line %d\n\n", __FILE__, __LINE__);
        }
        PartE[elems[i].procnum][PartSize[elems[i].procnum]++] = i; // "This element belongs to that processor"
    }
    fclose(inputmeshfile);
    printf("Finished reading elementss\n");

    printf("\n\nFollowing data was read:");
    printf("%d %d\n", nnodes, nelems);
    for (i=0; i<nnodes; i++){
        printf("%f %f %d\n", nodes[i].xy[0], nodes[i].xy[1], nodes[i].global_dof);
    }
    for (i=0; i<nelems; i++){
        printf("%d %d %d %d\n", elems[i].vertex[0], elems[i].vertex[1], elems[i].vertex[2], elems[i].procnum);
    }

    for(i=0;i<numprocs;i++){
        for(j=0;j<nnodes;j++){
            nodes[j].local_node = -1;
            PartN[i][j] = 0;
        }
        for(j=0;j<PartSize[i];j++){
            for(k=0;k<NDONTRI;k++){
                currentvertex = elems[PartE[i][j]].vertex[k];
                if (nodes[currentvertex].local_node == -1){
                    nodes[currentvertex].local_node = counter[i]++;
                }
                if (nodes[currentvertex].global_dof > -1 && PartN[i][currentvertex] == 0){
                    nodes[currentvertex].type++;
                    PartN[i][currentvertex] = 1;
                }
            }
        }
    }
}
    
void partition(NODE_STRUCT *nodes, ELEMENT_STRUCT *elems, int nnodes, int nelems, int numprocs,
               char *input_mesh_filename, int *counter, int *common, int *PartSize, int **PartN, int **PartE)
{
    int i, j, k, currentvertex;  // Counters
    char output_mesh_filename[MAXLINE];
    FILE *outfile;

    for(i=0;i<numprocs;i++){
        sprintf(output_mesh_filename,"%s_part_%i.2dm", input_mesh_filename, i);
        outfile = fopen(output_mesh_filename,"w");
        fprintf(outfile," %d %d\n", counter[i],PartSize[i]);
        counter[i] = 0;
        for (j=0;j<nnodes;j++){
            nodes[j].local_node = -1;
        }
        for (j=0;j<PartSize[i];j++){
            for(k=0;k<NDONTRI;k++){
                currentvertex = elems[PartE[i][j]].vertex[k];
                if (nodes[currentvertex].local_node == -1){
                    nodes[currentvertex].local_node = counter[i]++;
                    fprintf(outfile, " %.15e %.15e %10i\n", nodes[currentvertex].xy[0], nodes[currentvertex].xy[1],
                                                 nodes[currentvertex].type);
                }
            }
        }
        for (j=0;j<PartSize[i];j++){
            for(k=0;k<NDONTRI;k++){
                currentvertex = elems[PartE[i][j]].vertex[k];
                fprintf(outfile,"%10i",nodes[currentvertex].local_node);
            }
            fprintf(outfile,"\n");
        }
        for(j=0;j<numprocs;j++){
            common[j] = 0;
            if(i!=j){
                for(k=0; k<nnodes; k++){
                    if (PartN[i][k] > 0 && PartN[j][k] > 0){
                        common[j]++;
                    }
                }
            }
        }

        for(j=0;j<numprocs;j++){
            fprintf(outfile, " %i\n", common[j]);
            if (common[j] > 0){
                for(k=0; k<nnodes; k++){
                    if(PartN[i][k] > 0 && PartN[j][k] > 0){
                        fprintf(outfile, "%10i", nodes[k].local_node);
                    }
                }
                fprintf(outfile,"\n");
            }
        }
        fclose(outfile);
    }
}
