typedef struct{
    int nnodes;
    int nelems;

    int Interior_nodes;
#ifdef _MPI
    int Interior_boundary_nodes;
    int *common;
    int *shared;                     // Length = Interior_boundary_nodes
    int **neighbors;                 // Irregular, #Rows = numprocs   #Columns/Row = common[Row#]
#endif

    // Stiffness matrix in 3 parts, 2 square and one rectangular
    double **interior_stiffness;     // Square       Size = Interior_nodes            x   Interior_nodes
#ifdef _MPI
    double **interaction_stiffness;  // Rectangular  Size = Interior_nodes            x   Interior_boundary_nodes
    double **constrained_stiffness;  // Square       Size = Interior_boundary_nodes   x   Interior_boundary_nodes
#endif





    // Force vector in 2 parts
    double *interior_forcing;        // Length = Interior_nodes
#ifdef _MPI
    double *constrained_forcing;     // Length = Interior_boundary_nodes
#endif




    // Solution vector in 2 parts
    double *interior_solution;       // Length = Interior_nodes
#ifdef _MPI
    double *constrained_solution;    // Length = Interior_boundary_nodes
#endif







#ifdef _MPI
    double *buf;                     // Length = max(common)
#endif



    NODE_STRUCT *nodes;
    ELEMENT_STRUCT *elems;



} MODEL_STRUCT;


void model_free(MODEL_STRUCT *, int, int, int);
void model_print(MODEL_STRUCT *, int, int);
void model_alloc_matrices(MODEL_STRUCT *);
void model_build_system(MODEL_STRUCT *);
