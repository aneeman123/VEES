#include "HashTable.h"


/* hashtable basic structure: a 2D linked h_list */


/* constructor allocates structure and h_list
 * of pointers to bucket h_lists */
HashTable::HashTable(long size) {
    int i; 

    length = size;
    numItems = 0;

    /* allocate array */ 
	buckets = new HashLink *[length]; 
	  //(h_listptr *) malloc(sizeof(h_listptr)*size);

   /* allocate h_lists at each bucket */
    for(i = 0; i < length; i++)
      buckets[i] = new HashLink();
}

/* Frees all the memory associated with h_listptr* argument
   and sets * lptr to NULL */
HashTable::~HashTable() {
  int i;
  for(i = 0; i < length; i++) {
    delete buckets[i];
  }
  delete [] buckets;
}


/* insert a key-value pair into the hash table.
 * may cause chaining at the bucket. Will not insert
 * existing key-value pair.
 * Precondition: this item does not already exist in table,
 * value must be allocated by caller.
 */
void HashTable::insert( unsigned long key, int value) {
  unsigned long hash;
  hash = key % length;
  int val = buckets[hash]->h_get_value_via_key(key);

  if (val != value) {
	buckets[hash]->h_insert_before_first( value, key );
	numItems ++; 
  }
  else printf("HashTable::insert; pair (%ld,%d) already in table\n",
			  key, value);
}

/* dump contents of all items keys to stdout */
void  HashTable::print_hashtable() {
  int i;
 
  for (i = 0; i < length; i++) {
    buckets[i]->print_h_list();
  }
}



/* hash key to bucket, have linked h_list check for key 
 * convention: 1 ist true, 0 is false
*/
int HashTable::contains_hash_key( unsigned long key) {
  unsigned long hash; 
  hash = key % length;
  return (buckets[hash]->h_contains_key(key));
}

/* empty entire bucket */
void HashTable::remove_key( unsigned long key) {
  long hash; 
  hash = key % length;
  delete buckets[hash];
  buckets[hash] = new HashLink();
}

/* returns 1 if found, 0 if not */
int HashTable::remove_item(unsigned long key, int value) {
  long hash;
  int found;

  found = 0;
  hash = key% length;
  found = buckets[hash]->h_remove_node(key, value);
  printf("finished list remove\n");
  if( !found ) {
    fprintf(stderr,"Error-key-value pair %ld,%d is not in hash table",
	    key,value);
  }
  if( found)
    numItems--;

  return found;
}

int HashTable::num_items() {
  return numItems;
}
/***?????????????????????????????????**/


/* returns a h_list of all values from the hash table,
 * in the order read */
int *  HashTable::value_list() {
  int i,j;
  int* word_list;
  
  word_list = (int*) malloc((numItems)*sizeof(int));
  j = 0;
  for(i = 0; i < length; i++) {
    /* if list is not empty, go through list and
       place items in superlist */
    if( buckets[i]->h_get_length() > 0 ) {
       buckets[i]->h_move_first( );
      while( !buckets[i]->h_off_end()) {
		word_list[j] =  buckets[i]->h_get_current_value();
		j++;
		buckets[i]->h_move_next( );
      }
    }
  } 
  return word_list;
}


/* hash to bucket, get value matching key in bucket list */
int HashTable::get_value( unsigned long key ) {
  unsigned long hash;
  int val;

  hash = key % length;
  
  /* hash leads to a linked list...*/
  val = buckets[hash]->h_get_value_via_key( key);
  /*
  if ( val == -1 ) {
    fprintf(stderr, "%ld not contained in hash table\n", key);
    }*/
  return val;
}




