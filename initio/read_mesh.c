#include "global_header.h"

void read_mesh(MODEL_STRUCT **mod, int myid, int numprocs, char *model_name)
{
    FILE *partfile;                                       // Input partitioned mesh file
    char part_filename[50];                               // Partition File name
    int i, j;                                             // Counters
    int nnodes_loc, nelems_loc;                           // Number of nodes and elements
    MODEL_STRUCT *model;
    model = *mod;
    model->Interior_nodes = 0;

#ifdef _MPI
    model->Interior_boundary_nodes = 0;
    sprintf(part_filename,"%s_part_%d.2dm",model_name,myid);
#else
    sprintf(part_filename,"%s.2dm",model_name);
#endif
    partfile = fopen(part_filename, "r");
    /************************************************************************************************/
    // Reading initial data
    fscanf(partfile,"%i %i", &nnodes_loc, &nelems_loc); // Local number of nodes and elements
    model->nnodes = nnodes_loc;
    model->nelems = nelems_loc;

    model->nodes = (NODE_STRUCT *)    malloc(nnodes_loc*sizeof(NODE_STRUCT));      // Allocating space for nodes
    model->elems = (ELEMENT_STRUCT *) malloc(nelems_loc*sizeof(ELEMENT_STRUCT));   // Allocating space for elements

    NODE_STRUCT *nodes = model->nodes;
    ELEMENT_STRUCT *elems = model->elems;
    /************************************************************************************************/
    // Reading node data
//    printf("\nReading Node data...\n");
    for(i=0;i<nnodes_loc;i++){
        fscanf(partfile, "%lf %lf %i", &(nodes[i].xy[0]), &(nodes[i].xy[1]), &(nodes[i].type));
        if (nodes[i].type == 1){
            nodes[i].local = model->Interior_nodes++;
        }
#ifdef _MPI
        else if (nodes[i].type > 1){
            nodes[i].local = model->Interior_boundary_nodes++;
        }
#endif
        else{
            nodes[i].local = -1;
        }
    }
#ifdef _MPI
    model->shared = (int *) malloc(sizeof(int) * model->Interior_boundary_nodes);
    for (i=0, model->Interior_boundary_nodes=0 ; i<nnodes_loc; i++){
        if (nodes[i].type > 1){
            model->shared[model->Interior_boundary_nodes++] = nodes[i].type;
        }
    }
#endif
    
    /************************************************************************************************/
    // Reading element data
//    printf("\nReading Element data...\n");
    for(i=0;i<nelems_loc;i++){
        fscanf(partfile, "%i %i %i", &(elems[i].vertex[0]), &(elems[i].vertex[1]), &(elems[i].vertex[2]));
        elems[i].procnum = myid;
        elems[i].area    = 0;
    }
    elem_calc_areas(elems, nodes, nelems_loc);

#ifdef _MPI
    int maxcommon    = 0;
    model->common    = (int *)  malloc(sizeof(int)   * numprocs);
    model->neighbors = (int **) malloc(sizeof(int *) * numprocs);
//    printf("\nCalculating model->common...\n");
    for (i=0; i<numprocs; i++){
        fscanf(partfile, "%d", &(model->common[i]));
        maxcommon = (model->common[i] > maxcommon ? model->common[i] : maxcommon);
        if (model->common[i]>0){
            model->neighbors[i] = (int *)malloc(sizeof(int) * model->common[i]);
        }
        for (j=0; j<model->common[i]; j++){
            fscanf(partfile, "%d", &(model->neighbors[i][j]));
        }
    }
#endif
    fclose(partfile);

#ifdef _MESSG
    message_barrier(MPI_COMM_WORLD);
    int proc_count;
    for(proc_count=0;proc_count<4;proc_count++){
        if(proc_count==myid){
            printf("\n" COLOR_RED "*********** Processor %d printing below:" COLOR_RESET, myid);
            printf("\nProcessor %i:\n\tInterior nodes = %i\n", myid, model->Interior_nodes);
            printf("\tInterior boundary nodes = %i\n", model->Interior_boundary_nodes);
        }
        message_barrier(MPI_COMM_WORLD);
    }
    model->buf = (double *) malloc(sizeof(double) * maxcommon);
    //exit(-1);
#else
    printf("\tInterior nodes = %i\n",  model->Interior_nodes);
#endif

}
