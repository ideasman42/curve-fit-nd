
struct Heap;
struct HeapNode;
typedef struct Heap Heap;
typedef struct HeapNode HeapNode;

typedef void (*HeapFreeFP)(void *ptr);

Heap        *HEAP_new(unsigned int tot_reserve);
bool         HEAP_is_empty(Heap *heap);
void         HEAP_free(Heap *heap, HeapFreeFP ptrfreefp);
void        *HEAP_node_ptr(HeapNode *node);
void         HEAP_remove(Heap *heap, HeapNode *node);
HeapNode    *HEAP_insert(Heap *heap, double value, void *ptr);
void        *HEAP_popmin(Heap *heap);
void         HEAP_clear(Heap *heap, HeapFreeFP ptrfreefp);
unsigned int HEAP_size(Heap *heap);
HeapNode    *HEAP_top(Heap *heap);
double       HEAP_top_value(Heap *heap);
double       HEAP_node_value(HeapNode *node);
