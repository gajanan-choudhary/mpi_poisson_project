/* main */
void main_initialize(MODEL_STRUCT *, int, char *, int, int);
void main_run       (MODEL_STRUCT *, int, int, int);
void main_finalize  (MODEL_STRUCT *, int, int, int);

/* structs*/
void elem_calc_areas(ELEMENT_STRUCT *, NODE_STRUCT *, int);

/*initio */
void read_mesh(MODEL_STRUCT **, int, int, char *);
void partition_mesh(NODE_STRUCT *, ELEMENT_STRUCT *, char *, int);

/* solver */
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
);

#ifdef _MPI
void update(double *a /*Vector to be updated */, int numprocs, double *buf, NODE_STRUCT *nodes,
            int **neighbors, int *common);
#endif

double apply_dirichlet_bc(double, double);

int parallel_conjugate_gradient(MODEL_STRUCT *model, int myid, int numprocs, double tol, int maxit);

