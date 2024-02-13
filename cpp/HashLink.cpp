/****************************************************/
/*  list.c                                          */
/*  Header file for Double-Link list                */
/*                                                  */
/****************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "HashLink.h"

/****************** CONSTRUCTORS-DISTRUCTORS ****************/
/* Returns a listptr which point to a new empty list */
HashLink::HashLink()
{

    length   = 0;
    current  = NULL;
    head     = NULL;
    tail     = NULL;

}

/* Frees all the memory associated with listptr* argument
   and sets  lptr to NULL */
HashLink::~HashLink()
{   
    /* walk through list deleting nodes */
    while(!h_is_empty() )
      h_delete_first();
}


/* Return the value going with this key.
 * Returns -1 if key doesn't exist. The assumption is
 * that several keys may hash to this bucket's list,
 * but the keys themselves are unique */
int HashLink::h_get_value_via_key(unsigned long key) {
  struct node * walker;

  if( h_is_empty() )
    return -1;
  else
    walker = head;

    while (walker) {
      if( walker->key == key )
        return walker->value;
      walker = walker->next;
    }
    return -1;
}

/* does this list contain an equivalent key?
 * Returns: boolean true or false */
int HashLink::h_contains_key( unsigned long key) {
  struct node * walker;
  if( h_is_empty() )
    return 0;
  else
    walker = head;
    while (walker) {
      if( walker->key == key )
        return 1;
      walker = walker->next;
    }
    return 0;
}


/* if you have an exact match, remove the node.
 * Both key and value must match */
int HashLink::h_remove_node(unsigned long key, int value) {
  struct node * walker;

  if( h_is_empty() )
    return 0;
  else {

    walker = head;
    while (walker) {

      if( walker->key == key && walker->value == value ) {
	
		if( walker == head ) {
		  h_delete_first(); //calls length--
		}
		else if ( walker == tail ) {
		  h_delete_last(); //calls length--
		}

		else {
		  walker->next->prev = walker->prev;
		  walker->prev->next = walker->next;		  
		  free(walker);		  
		  length--;		 
		}

        return 1;
      }

	  if (walker == tail) {
		return 0;
	  }
      walker = walker->next;
    }
  }
    return 0;
}



/******************* ACCESS FUNCTIONS *********************/
/* Returns 1 if list contains no element */
int HashLink::h_is_empty()
{
    return (length == 0);
}

/* Returns 1 if current marker referrs to the first element */
int HashLink::h_at_first()
{
  if (current->prev == NULL)
	return 1;
  return 0;
}

/* Returns 1 if current marker referrs to the last element */
int HashLink::h_at_last()
{
  if (current->next == NULL) 
	return 1;
  return 0;
}

/* Returns 1 if no element is current */
int HashLink::h_off_end() 
{
    return (current == NULL);
}

/* Returns the first element in the list (Pre: !is_empty) */
int HashLink::h_get_first() 
{
    if (h_is_empty()) { 
	  printf("h_get_first called on an empty linked list\n");
	  return -1;
	}
    return head->value;
}

/* Returns the last element in the list (Pre: !is_empty) */
int HashLink::h_get_last()
{
    if (h_is_empty()) { 
	  printf("h_get_last called on an empty linked list\n");
	  return -1;
	}
    return tail->value;
}

/* Returns the current key in the list (Pre: !is_empty Post: !off_end) */
unsigned long HashLink::h_get_current_key()
{
    if (h_is_empty() || current == NULL) { 
      printf("h_get_current_key called on an empty linked list\n");
	  return 0;
    }
    return current->key;
}

/* Returns the current value in the list (Pre: !is_empty Post: !off_end) */
int HashLink::h_get_current_value()
{
   
    if (h_is_empty() || current == NULL) { 
	  printf("h_get_current_key called on an empty linked list\n");
	  return -1;
    }
    return current->value;
}

/* Returns the length of LIST */
int HashLink::h_get_length() 
{
    return length;
}

/**************** MANUPULATION FUNCTIONS ******************/
/* Makes first element current 
 * Pre: !is_empty
 * post: !off_end */
void HashLink::h_move_first() 
{
    if (!h_is_empty()) 
        current = head;
}

/* Makes last element current 
 * Pre: !is_empty
 * post: !off_end */
void HashLink::h_move_last()
{
    if (!h_is_empty()) 
        current = tail;
}

/* Moves the current pointer towards the beginning of the list 
 * Pre: !is_empty &&  !off_end */
void HashLink::h_move_prev()
{
    if (!h_is_empty() && !h_off_end())
        current = current->prev;
}

/* Moves the current pointer towards the end of the list 
 * Pre: !is_empty && !off_end */
void HashLink::h_move_next() 
{
    if (!h_is_empty() && !h_off_end())
        current = current->next;
}

/* Adds a new element at the beginning of the list 
 * post: !is_empty */
void HashLink::h_insert_before_first(int data, unsigned long key_data) 
{
    struct node *new_node;
    
    if ((new_node = (node *)malloc(sizeof(node))) == NULL) {
        printf("A malloc error in insert_before_h_list()!\n");
    } 
    new_node->key =  key_data;
    new_node->value=     data;
    new_node->prev =     NULL;
    new_node->next =     NULL;

    if (h_is_empty()) {
        head = tail = new_node; 
    }
    else { 
        new_node->next   = head;
        head->prev = new_node; 
        head       = new_node; 
    }
    length++;
}



/* Deletes the first element
 * Pre: !is_empty */
void HashLink::h_delete_first()
{
    node *temp;

    if (!h_is_empty()) {
        if (h_get_length() == 1) {
            free(head);
            length  = 0; 
            head    = NULL;
            tail    = NULL;
            current = NULL; 
        }
        else {
            temp = head;
            head->next->prev = NULL; 
            head             = head->next;

			free (temp);
            temp = NULL;
            length--;
        }
    }
}

/* Deletes the last element
 * Pre: !is_empty */
void HashLink::h_delete_last()
{
    node *temp;

    if (!h_is_empty()) {
        if (h_get_length() == 1) {
            free(tail);
            length  = 0;
            head    = NULL;
            tail    = NULL;
            current = NULL;
        }
        else {
            temp = tail;
            tail->prev->next = NULL;
            tail = tail->prev;

			free (temp);
            temp = NULL;
            length--;
        }
    }
}

/* Deletes the current element
 * Pre: !is_empty, !off_end
 * Post: off_end */
void HashLink::h_delete_current()
{
    node *temp;

    if (!h_is_empty() && !h_off_end()) {
        if (h_at_first())
            h_delete_first();
       
        else if (h_at_last()) 
            h_delete_last();
        
        else if (!h_at_first() && !h_at_last()) {
            temp = current;
            current->next->prev = current->prev;
            current->prev->next = current->next;
	 
            free (temp);
            temp = NULL;  
        }
        length--;
        current = NULL;
    }//non-empty list
}



/* Print out the list elments the input list */
void HashLink::print_h_list() 
{
    struct node *walker;

    walker = head;
    while (walker) {
        printf("(%ld,%d) \n",walker->key, walker->value);
        walker = walker->next;
    }
   
}

