#ifndef BCUTILS_ARRAY
#define BCUTILS_ARRAY

#include <assert.h>
#include <stdlib.h>

#define INITIAL_ARRAY_SIZE 25

// A simple dynamically growing array data structure
typedef struct
{
  void **data;
  size_t capacity;
  size_t usage;
} Array;

// Add an element to the array
void array_add(Array *r, void *e);

// Destructor: free the memory previously used by the array
void array_delete(Array *r, void (*valuefreefunc)(void *));

// Get the element at the requested index
void *array_get(Array *r, int index);

// Constructor: allocate memory for a new array
Array *array_new();

// Determine the number of elements in the array
size_t array_size(Array *r);

Array *array_concat(Array *first, Array *second);

#endif
