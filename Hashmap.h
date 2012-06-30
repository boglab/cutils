// Implementation extends hashtable description at
// http://www.sparknotes.com/cs/searching/hashmaps/section3.rhtml

#ifndef BCUTILS_HASHMAP
#define BCUTILS_HASHMAP

#include <stdlib.h>
#include <string.h>

// Item in a hash table, with a key/value pair mapping
// the value to a key in the table
typedef struct _hashmap_item_
{
  char *key;
  void *value;
  struct _hashmap_item_ *next;
} HashmapItem;

// Simple hash map structure
typedef struct
{
  int size;
  int usage;
  HashmapItem **table;
} Hashmap;

// Add an item to the hashmap
int hashmap_add(Hashmap *hash, char *key, void *value);

// Destructor: free the memory used by this hashmap
void hashmap_delete(Hashmap *hash, void (*valuefreefunc)(void *));

// Get an item from the hashmap
void *hashmap_get(Hashmap *hash, char *key);

// Hash function
unsigned int hashmap_hash(Hashmap *hash, char *key);

// Get all of the keys associated with this hashmap (user is responsible for freeing memory)
char **hashmap_keys(Hashmap *hash);

// Constructor: allocate some memory for a hashmap
Hashmap *hashmap_new(int size);

// Get the number of key/value pairs in the hashmap
int hashmap_size(Hashmap *hash);

#endif
