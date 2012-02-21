#include "Array.h"

void array_add(Array *r, void *e)
{
  if(r->usage == r->capacity)
  {
    r->capacity *= 2;
    r->data = (void **)realloc(r->data, r->capacity * sizeof(void *));
  }
  r->data[r->usage++] = e;
}

void array_delete(Array *r)
{
  free(r->data);
  r->data = NULL;
  r->capacity = r->usage = 0;
  free(r);
  r = NULL;
}

void *array_get(Array *r, int index)
{
  assert(index < r->usage);
  return r->data[index];
}

Array *array_new()
{
  Array *r = (Array *)malloc( sizeof(Array) );
  r->capacity = INITIAL_ARRAY_SIZE;
  r->usage = 0;
  r->data = (void **)malloc(r->capacity * sizeof(void *));
  return r;
}

size_t array_size(Array *r)
{
  return r->usage;
}
