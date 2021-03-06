/* Standard Header Files*/
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include <time.h>
#include <stdlib.h>
//#include <stddef.h>
//#include <limits.h>

/* Coloring text! */
#ifdef _COLOR_OUTPUT
#define COLOR_RED     "\x1b[31m"
#define COLOR_GREEN   "\x1b[32m"
#define COLOR_YELLOW  "\x1b[33m"
#define COLOR_BLUE    "\x1b[34m"
#define COLOR_MAGENTA "\x1b[35m"
#define COLOR_CYAN    "\x1b[36m"
#define COLOR_RESET   "\x1b[0m"
#else
#define COLOR_RED ""
#define COLOR_GREEN ""
#define COLOR_YELLOW ""
#define COLOR_BLUE ""
#define COLOR_MAGENTA ""
#define COLOR_CYAN ""
#define COLOR_RESET ""
#endif

/*Dimension of the problem*/
#define ndim 2
#define NDONTRI 3

#define MAXLINE 300

/*Bools*/
#define ON 1
#define OFF 0
#define YES 1
#define NO -3
#define TRUE 1
#define FALSE 0
#define USED 1
#define UNUSED 0
#define TINY DBL_EPSILON    /* Smallest double */
#define SMALL FLT_EPSILON   /* Smallest float */
#define SMALL1 1.0E-1
#define SMALL2 1.0E-2
#define SMALL3 1.0E-3
#define SMALL4 1.0E-4
#define SMALL5 1.0E-5
#define SMALL6 1.0E-6
#define SMALL7 1.0E-7
#define SMALL8 1.0E-8


