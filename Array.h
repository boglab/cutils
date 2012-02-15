#ifndef OTFS_ARRAY
#define OTFS_ARRAY

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
void array_add(Array *r, void *e)
{
  if(r->usage == r->capacity)
  {
    r->capacity *= 2;
    r->data = (void **)realloc(r->data, r->capacity * sizeof(void *));
  }
  r->data[r->usage++] = e;
}

// Destructor: free the memory previously used by the array
void array_delete(Array *r)
{
  free(r->data);
  r->data = NULL;
  r->capacity = r->usage = 0;
  free(r);
  r = NULL;
}

// Get the element at the requested index
void *array_get(Array *r, int index)
{
  assert(index < r->usage);
  return r->data[index];
}

// Constructor: allocate memory for a new array
Array *array_new()
{
  Array *r = (Array *)malloc( sizeof(Array) );
  r->capacity = INITIAL_ARRAY_SIZE;
  r->usage = 0;
  r->data = (void **)malloc(r->capacity * sizeof(void *));
  return r;
}

// Determine the number of elements in the array
size_t array_size(Array *r)
{
  return r->usage;
}

#endif
