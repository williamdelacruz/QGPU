#include "linkedlist.h"


/*
 *   Create a new node
 */
 
nodeList *new_node(int index, float q, nodeList *list_nodes)
{
    list_nodes[index].index = index;
    list_nodes[index].q     = q;
    list_nodes[index].sig   = NULL;
    return &list_nodes[index];
}


/*
 *   Insert at the beginning of the list
 */
int insertIni(int index, float q, nodeList **list, nodeList **last, nodeList *list_nodes)
{
    nodeList *nuevo = new_node(index, q, list_nodes);
    
    if (*list==NULL) {
        *list = nuevo;
        *last = *list;
    }
    else {
        nuevo->sig = *list;
        *list = nuevo;
    }
    
    return 0;
}



/*
 *   Insert at the beginning of the list
 */
int insertIni_(int index, nodeList **list, nodeList *list_nodes)
{
    nodeList *nuevo = new_node(index, 0, list_nodes);
    
    if (*list==NULL)
        *list = nuevo;
    else {
        nuevo->sig = *list;
        *list = nuevo;
    }
    
    return 0;
}



/*
 *   Insert at the end of the list
 */
int insertEnd(int index, float q, nodeList **list, nodeList **last, nodeList *list_nodes)
{
    nodeList *nuevo = new_node(index, q, list_nodes);
    
    if (*list==NULL) {
        *list = nuevo;
        *last = *list;
    }
    else {
        (*last)->sig = nuevo;
        *last = nuevo;
    }
    
    return 0;
}



/*
 *    Remove the first element of the list
 */
int removeFirst(int *index, nodeList **list, nodeList **last)
{
    if (*list==NULL)
        return -1;
    else {
        *index = (*list)->index;
        *list = (*list)->sig;
        if (*list==NULL)
            *last = NULL;
    }
    
    return 0;
}


/*
 *    Insert items and sort
 */
int insertSort(int index, float q, nodeList **list, nodeList **last, nodeList *list_nodes)
{
    nodeList *nuevo;
    nodeList *ptr1, *ptr2;
    int flag;
    
    // list is empty
    if (*list==NULL)
    {
        nuevo = new_node(index, q, list_nodes);
        *list = nuevo;
        *last = *list;
    }
    else {
        // insert at the beginning
        if (q >= (*list)->q)
        {
            insertIni(index, q, list, last, list_nodes);
        }
        else
        {
            // insert at the end
            if (q <= (*last)->q)
            {
                insertEnd(index, q, list, last, list_nodes);
            }
            else
            {   // insert in the meedle
                ptr1 = *list;
                ptr2 = *list;
                
                while (ptr1->sig!=NULL && q <= ptr1->q) {
                    ptr2 = ptr1;
                    ptr1 = ptr1->sig;
                }
                
                nuevo = new_node(index, q, list_nodes);
                ptr2->sig = nuevo;
                nuevo->sig = ptr1;
            }
        }
    }
    
    return 0;
}


/*
 *   Insertion in a linked list using a pivot (sort)
 */
int insertSortPivot(int index, float q, nodeList **list, nodeList **last, nodeList **pivot, nodeList *list_nodes, int *num_nodes)
{
    nodeList *ptr, *ptr2;
    nodeList *nuevo = new_node(index, q, list_nodes);
    
    
    // appears the pivot
    
    if (*num_nodes==3 && *list!=NULL)
        *pivot = (*list)->sig;
    
    
    if (*list==NULL) {
        *list = nuevo;
        *last = *list;
        *pivot = NULL;
        (*num_nodes)++;
    }
    else
    {
        // insert at the beggining
        if (q >= (*list)->q)
        {
            //printf("beggining\n");
            nuevo->sig = *list;
            (*list)->prev = nuevo;
            *list = nuevo;
            (*num_nodes)++;
            
            if (*num_nodes==3)
                *pivot = (*list)->sig;
            else if (*num_nodes>3 && *num_nodes%2==1) {
                if (*pivot!=NULL)
                    *pivot = (*pivot)->prev;
            }
        }
        else
        {
            // insert at the end
            if (q <= (*last)->q)
            {
                //printf("end\n");
                (*last)->sig = nuevo;
                nuevo->prev = *last;
                *last = nuevo;
                (*num_nodes)++;
                
                if (*num_nodes==3)
                    *pivot = (*list)->sig;
                else if (*num_nodes>3 && *num_nodes%2==1) {
                    if (*pivot!=NULL)
                        *pivot = (*pivot)->sig;
                }
            }
            else // insert sort
            {
                // insert in the middle
                if (*num_nodes<=2)
                {
                    ptr = (*list)->sig;
                    (*list)->sig = nuevo;
                    nuevo->prev = *list;
                    nuevo->sig = ptr;
                    ptr->prev = nuevo;
                    (*num_nodes)++;
                }
                else // There exists a pivot
                {
                    if (*pivot==NULL)
                        *pivot = (*list)->sig;
                    
                    // search from the pivot
                    if ((*pivot)->q >= q)
                    {
                        for (ptr=*pivot; ptr->sig!=NULL && ptr->q>=q; ptr=ptr->sig);
                        
                        if (ptr->q < q) ptr = ptr->prev;
                        
                        // link
                        ptr2 = ptr->sig;
                        ptr->sig = nuevo;
                        nuevo->prev = ptr;
                        nuevo->sig = ptr2;
                        ptr2->prev = nuevo;
                        (*num_nodes)++;
                        
                        if (*num_nodes>2 && *num_nodes%2==1)
                            *pivot = (*pivot)->sig;
                    }
                    else  // search from the beginning of the list
                    {
                        for (ptr=*list; ptr->sig!=NULL && ptr->q>=q; ptr=ptr->sig);
                        
                        if (ptr->q < q) ptr = ptr->prev;
                        
                        // link
                        ptr2 = ptr->sig;
                        ptr->sig = nuevo;
                        nuevo->prev = ptr;
                        nuevo->sig = ptr2;
                        ptr2->prev = nuevo;
                        (*num_nodes)++;
                        
                        if (*num_nodes>2 && *num_nodes%2==1)
                            *pivot = (*pivot)->prev;
                    }
                }
            }
        }
    }
    
    return 0;
}



/*
 *    Remove the first element of the list
 */
int removeFirstPivot(int *index, nodeList **list, nodeList **last,
                     nodeList **pivot, int *num_nodes)
{
    nodeList *ptr;
    
    if (*list==NULL)
        return -1;
    else {
        if ((*list)->sig==NULL) {
            ptr = *list;
            *index = ptr->index;
            // link
            *list = NULL;
            *last = *list;
            *pivot = NULL;
            (*num_nodes)--;
        }
        else {
            ptr = *list;
            *index = ptr->index;
            // link
            *list = (*list)->sig;
            (*list)->prev = NULL;
            (*num_nodes)--;
            
            if (*num_nodes>2 && *num_nodes%2==1) {
                if (*pivot!=NULL)
                    *pivot = (*pivot)->sig;
            }
            else if (*num_nodes<3)
                *pivot = NULL;
        }
    }
    
    return 0;
}