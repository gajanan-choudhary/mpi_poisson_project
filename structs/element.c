#include "global_header.h"

void elem_calc_areas(ELEMENT_STRUCT *elems, NODE_STRUCT *nodes, int nelems){
    int i;
    for (i=0;i<nelems;i++){
        elems[i].area = (  nodes[elems[i].vertex[1]].xy[0] * nodes[elems[i].vertex[2]].xy[1]
                         + nodes[elems[i].vertex[2]].xy[0] * nodes[elems[i].vertex[0]].xy[1]
                         + nodes[elems[i].vertex[0]].xy[0] * nodes[elems[i].vertex[1]].xy[1]  )
                      - (  nodes[elems[i].vertex[2]].xy[0] * nodes[elems[i].vertex[1]].xy[1]
                         + nodes[elems[i].vertex[0]].xy[0] * nodes[elems[i].vertex[2]].xy[1]
                         + nodes[elems[i].vertex[1]].xy[0] * nodes[elems[i].vertex[0]].xy[1]  );
        elems[i].area /= 2;
    }
}
