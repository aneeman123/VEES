
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "HashLink.h"
#ifndef _HASH_
#define _HASH_

typedef struct hash *hashptr;
/****************** CONSTRUCTORS-DISTRUCTORS ****************/
/* Returns a hashptr which point to a new hashtable with 
'size' buckets. Each bucket is a linked h_list */

class HashTable {

 public:
  HashTable( long size );
  
  
  /* Frees all the memory associated with h_lists in hashtable
	 and ahstable structure itself */
  ~HashTable();
  
  /***********************************************/
  
  /* get all key words from hash table */
  unsigned long ** key_list(); 
  
  /* get all value words from hash table */
  int* value_list(); 
  
  /* get number of items that have been inserted
   * into the hash table */
  int num_items();
  
  /* hashes key, inserts key AND VALUE into hash entry */
  void insert(unsigned long key, int value);
  
  /* removes h_list at hash */
  void remove_key( unsigned long key);
  
  /* get the value at the key supplied */
  int get_value( unsigned long key );
  
  /* prints each h_list in hash table */
  void print_hashtable();
  
  /* print only items not in filter h_list */
  //void print_filtered_hash(int argc, char ** argv, hashptr hptr);
  
  /* contains an entry that hashes to this key */
  int contains_hash_key( unsigned long key);
  
  
  
  /* Remove the exact key-value pair from the hash table.
   * Returns: boolean indicated whether item was in table to 
   * begin with.
   */
  int remove_item(unsigned long key, int value);
  
  
 private:
  HashLink ** buckets; /* an array of buckets */
  unsigned int length;
  unsigned int numItems; /* number of items currently stored */
};

#endif
