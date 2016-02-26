typedef struct {
    int node[2];
    double area;
    double length;
    double E;         /* Material Young's Modulus */
    double **K;   /* Element Stiffness Matrix, 4x4 for 2D */
    double *F;
} ELEMENT;
