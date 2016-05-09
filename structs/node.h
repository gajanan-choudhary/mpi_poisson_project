typedef struct{
//    int global_id;  // Global node numbers, starting from 0 each
    double xy[2];
    int type;       // Interior = 0, Exterior Boundary < 0, Interior Boundary > 0
    int local;
} NODE_STRUCT;
