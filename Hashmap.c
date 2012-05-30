#include "Hashmap.h"

int hashmap_add(Hashmap *hash, char *key, void *value) {
  HashmapItem *item;
  unsigned int hashval;

  if (hashmap_get(hash, key) != NULL)
    return 0;

  hashval = hashmap_hash(hash, key);
  item = malloc(sizeof(HashmapItem));
  if (item == NULL)
    return 0;
  item->key = strdup(key);
  item->value = value;
  item->next = hash->table[hashval];
  hash->table[hashval] = item;
  hash->usage += 1;

  return 1;
}

void hashmap_delete(Hashmap *hash, void (*valuefreefunc)(void *)) {
  HashmapItem *item, *temp;
  int i;

  if (hash == NULL)
    return;

  for (i = 0; i < hash->size; i++) {
    item = hash->table[i];
    while (item != NULL) {
      temp = item;
      item = item->next;
      free(temp->key);
      if (valuefreefunc != NULL)
        valuefreefunc(temp->value);
      free(temp);
    }
  }

  free(hash->table);
  free(hash);
}

void *hashmap_get(Hashmap *hash, char *key) {
  HashmapItem *item;
  unsigned int hashval = hashmap_hash(hash, key);
  for (item = hash->table[hashval]; item != NULL; item = item->next) {
    if (strcmp(key, item->key) == 0)
      return item->value;
  }
  return NULL;
}

unsigned int hashmap_hash(Hashmap *hash, char *key) {
  unsigned int hashval = 0;
  for (; *key != '\0'; key++)
    hashval = *key + (hashval << 5) - hashval;
  return hashval % hash->size;
}

char **hashmap_keys(Hashmap *hash) {
  char **keys;
  HashmapItem *item;
  int i, j;

  keys = malloc(sizeof(char *) * hash->usage);
  j = 0;
  for (i = 0; i < hash->size; i++) {
    item = hash->table[i];
    while (item != NULL) {
      keys[j++] = item->key;
      item = item->next;
    }
  }

  return keys;
}

Hashmap *hashmap_new(int size) {
  Hashmap *hash;
  int i;

  if (size < 1)
    return NULL;

  hash = malloc(sizeof(HashmapItem));
  if (hash == NULL)
    return NULL;

  hash->table = malloc(sizeof(HashmapItem *) * size);
  if (hash->table == NULL)
    return NULL;

  for (i = 0; i < size; i++)
    hash->table[i] = NULL;

  hash->size = size;
  hash->usage = 0;
  return hash;
}

int hashmap_size(Hashmap *hash) {
  return hash->usage;
}
