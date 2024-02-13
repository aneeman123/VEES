/****************************************************/
/*  list.h                                          */
/*  Header file for Double-Link unordered list      */
/*                                                  */  
/****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifndef _HASH_INT_LIST_
#define _HASH_INT_LIST_

typedef struct node 
{
    unsigned long     key;
    int             value;

    struct node	 *   prev;
    struct node	    *next;
} node;

class HashLink {

public:
  /****************** CONSTRUCTORS-DISTRUCTORS ****************/
  /* Returns a h_listptr which point to a new empty h_list */
  HashLink();
  /* Frees all the memory associated with h_listptr* argument
	 and sets * lptr to NULL */
  ~HashLink();
  
  
  int h_remove_node( unsigned long key, int value);

  /******************* ACCESS FUNCTIONS *********************/
  /* Returns true if h_list contains no element */
  int h_is_empty(); 
  
  /* get value at the given key, if key exists */
  int h_get_value_via_key( unsigned long key);
  
  /* declare whether h_list contains key */
  int h_contains_key( unsigned long key);
  
  /* Returns the first element in the h_list (Pre: !is_empty) */
  int h_get_first();
  
  /* Returns the last element in the h_list (Pre: !is_empty) */
  int h_get_last();
  
  /* Returns the current key in the h_list (Pre: !is_empty Post: !off_end) */
  unsigned long h_get_current_key();
  
  /* Returns the current value in the h_list (Pre: !is_empty Post: !off_end) */
  int h_get_current_value();
  
  /* Returns the length of H_LIST */
  int h_get_length();
  
  /************ MANIPULATION ********************/
/* Adds a new element at the beginning of the h_list 
 * post: !is_empty */
  void h_insert_before_first(int data, unsigned long key_data);
  
  /* Deletes the first element
   * Pre: !is_empty */
  void h_delete_first();
  
  /* Deletes the last element
   * Pre: !is_empty */
  void h_delete_last();
  
  /* Deletes the current element
   * Pre: !is_empty, !off_end
   * Post: off_end */
void h_delete_current();
  
  /************* TRAVERSAL *******************************/
  /* Returns true if current marker referrs to the first element */
  int h_at_first();
  
  /* Returns true if current marker referrs to the last element */
  int h_at_last();
  
  /* Returns true if no element is current */
  int h_off_end();
  
/* Makes first element current 
 * Pre: !is_empty
 * post: !off_end */
  void h_move_first();
  
  /* Makes last element current 
   * Pre: !is_empty
   * post: !off_end */
  void h_move_last();
  
  /* Moves the current pointer towards the beginning of the h_list 
   * Pre: !is_empty &&  !off_end */
  void h_move_prev();
  
  /* Moves the current pointer towards the end of the h_list 
   * Pre: !is_empty && !off_end */
  void h_move_next();
  
  /**************** OTHERS OPERATIONS ******************/
  
  
  /* Print out the h_list keys the input h_list */
  void print_h_list();
  
private:

    int          length;
    struct node	 *current;
    struct node  *head;
    struct node	 *tail;
};

#endif
