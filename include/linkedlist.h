#ifndef LINKEDLIST_H
#define LINKEDLIST_H

#include <stdlib.h>
#include <stdio.h>

/*
 *   Data structure of a node
 */
typedef struct nodeList_ {
    int index;
    float q;
    struct nodeList_ *sig, *prev;
} nodeList;


// Prototypes


nodeList *new_node(int, float, nodeList *);
int insertIni(int, float, nodeList **, nodeList **, nodeList *);
int insertIni_(int, nodeList **, nodeList *);
int insertEnd(int, float, nodeList **, nodeList **, nodeList *);
int removeFirst(int *, nodeList **, nodeList **);
int insertSort(int, float, nodeList **, nodeList **, nodeList *);

// linkedlist with pivot
int insertSortPivot(int, float, nodeList **, nodeList **, nodeList **, nodeList *, int *);
int removeFirstPivot(int *, nodeList **, nodeList **, nodeList **, int *);


#endif // LINKEDLIST_H
