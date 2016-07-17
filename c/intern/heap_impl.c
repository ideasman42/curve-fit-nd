
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>

#include "heap_impl.h"

/* swap with a temp value */
#define SWAP_TVAL(tval, a, b)  {  \
	(tval) = (a);                 \
	(a) = (b);                    \
	(b) = (tval);                 \
} (void)0

#ifdef __GNUC__
//#  define LIKELY(x)       __builtin_expect(!!(x), 1)
#  define UNLIKELY(x)     __builtin_expect(!!(x), 0)
#else
//#  define LIKELY(x)       (x)
#  define UNLIKELY(x)     (x)
#endif

/***/

struct HeapNode {
	void        *ptr;
	double       value;
	unsigned int index;
};

struct HeapNode_Chunk {
	struct HeapNode_Chunk *prev;
	unsigned int    size;
	unsigned int    bufsize;
	struct HeapNode buf[0];
};

/**
 * Number of nodes to include per #HeapNode_Chunk when no reserved size is passed,
 * or we allocate past the reserved number.
 *
 * \note Optimize number for 64kb allocs.
 */
#define HEAP_CHUNK_DEFAULT_NUM \
	(((1 << 16) - sizeof(struct HeapNode_Chunk)) / sizeof(HeapNode))

struct Heap {
	unsigned int size;
	unsigned int bufsize;
	HeapNode **tree;

	struct {
		/* Always keep at least one chunk (never NULL) */
		struct HeapNode_Chunk *chunk;
		/* when NULL, allocate a new chunk */
		HeapNode *free;
	} nodes;
};

/** \name Internal Functions
 * \{ */

#define HEAP_PARENT(i) (((i) - 1) >> 1)
#define HEAP_LEFT(i)   (((i) << 1) + 1)
#define HEAP_RIGHT(i)  (((i) << 1) + 2)
#define HEAP_COMPARE(a, b) ((a)->value < (b)->value)

#if 0  /* UNUSED */
#define HEAP_EQUALS(a, b) ((a)->value == (b)->value)
#endif

static void heap_swap(Heap *heap, const unsigned int i, const unsigned int j)
{

#if 0
	SWAP(unsigned int,  heap->tree[i]->index, heap->tree[j]->index);
	SWAP(HeapNode *,    heap->tree[i],        heap->tree[j]);
#else
	HeapNode **tree = heap->tree;
	union {
		unsigned int  index;
		HeapNode     *node;
	} tmp;
	SWAP_TVAL(tmp.index, tree[i]->index, tree[j]->index);
	SWAP_TVAL(tmp.node,  tree[i],        tree[j]);
#endif
}

static void heap_down(Heap *heap, unsigned int i)
{
	/* size won't change in the loop */
	const unsigned int size = heap->size;

	while (1) {
		const unsigned int l = HEAP_LEFT(i);
		const unsigned int r = HEAP_RIGHT(i);
		unsigned int smallest;

		smallest = ((l < size) && HEAP_COMPARE(heap->tree[l], heap->tree[i])) ? l : i;

		if ((r < size) && HEAP_COMPARE(heap->tree[r], heap->tree[smallest])) {
			smallest = r;
		}

		if (smallest == i) {
			break;
		}

		heap_swap(heap, i, smallest);
		i = smallest;
	}
}

static void heap_up(Heap *heap, unsigned int i)
{
	while (i > 0) {
		const unsigned int p = HEAP_PARENT(i);

		if (HEAP_COMPARE(heap->tree[p], heap->tree[i])) {
			break;
		}
		heap_swap(heap, p, i);
		i = p;
	}
}

/** \} */


/** \name Internal Memory Management
 * \{ */

static struct HeapNode_Chunk *heap_node_alloc_chunk(
        unsigned int tot_nodes, struct HeapNode_Chunk *chunk_prev)
{
	struct HeapNode_Chunk *chunk = malloc(
	        sizeof(struct HeapNode_Chunk) + (sizeof(HeapNode) * tot_nodes));
	chunk->prev = chunk_prev;
	chunk->bufsize = tot_nodes;
	chunk->size = 0;
	return chunk;
}

static struct HeapNode *heap_node_alloc(Heap *heap)
{
	HeapNode *node;

	if (heap->nodes.free) {
		node = heap->nodes.free;
		heap->nodes.free = heap->nodes.free->ptr;
	}
	else {
		struct HeapNode_Chunk *chunk = heap->nodes.chunk;
		if (UNLIKELY(chunk->size == chunk->bufsize)) {
			struct HeapNode_Chunk *chunk_next = heap_node_alloc_chunk(HEAP_CHUNK_DEFAULT_NUM, chunk);
			chunk = chunk_next;
		}
		node = &chunk->buf[chunk->size++];
	}

	return node;
}

static void heap_node_free(Heap *heap, HeapNode *node)
{
	node->ptr = heap->nodes.free;
	heap->nodes.free = node;
}

/** \} */


/** \name Public Heap API
 * \{ */

/* use when the size of the heap is known in advance */
Heap *HEAP_new(unsigned int tot_reserve)
{
	Heap *heap = malloc(sizeof(Heap));
	/* ensure we have at least one so we can keep doubling it */
	heap->size = 0;
	heap->bufsize = tot_reserve ? tot_reserve : 1;
	heap->tree = malloc(heap->bufsize * sizeof(HeapNode *));

	heap->nodes.chunk = heap_node_alloc_chunk((tot_reserve > 1) ? tot_reserve : HEAP_CHUNK_DEFAULT_NUM, NULL);
	heap->nodes.free = NULL;

	return heap;
}

void HEAP_free(Heap *heap, HeapFreeFP ptrfreefp)
{
	if (ptrfreefp) {
		unsigned int i;

		for (i = 0; i < heap->size; i++) {
			ptrfreefp(heap->tree[i]->ptr);
		}
	}

	struct HeapNode_Chunk *chunk = heap->nodes.chunk;
	do {
		struct HeapNode_Chunk *chunk_prev;
		chunk_prev = chunk->prev;
		free(chunk);
		chunk = chunk_prev;
	} while (chunk);

	free(heap->tree);
	free(heap);
}

void HEAP_clear(Heap *heap, HeapFreeFP ptrfreefp)
{
	if (ptrfreefp) {
		unsigned int i;

		for (i = 0; i < heap->size; i++) {
			ptrfreefp(heap->tree[i]->ptr);
		}
	}
	heap->size = 0;

	/* Remove all except the last chunk */
	while (heap->nodes.chunk->prev) {
		struct HeapNode_Chunk *chunk_prev = heap->nodes.chunk->prev;
		free(heap->nodes.chunk);
		heap->nodes.chunk = chunk_prev;
	}
	heap->nodes.chunk->size = 0;
	heap->nodes.free = NULL;
}

HeapNode *HEAP_insert(Heap *heap, double value, void *ptr)
{
	HeapNode *node;

	if (UNLIKELY(heap->size >= heap->bufsize)) {
		heap->bufsize *= 2;
		heap->tree = realloc(heap->tree, heap->bufsize * sizeof(*heap->tree));
	}

	node = heap_node_alloc(heap);

	node->ptr = ptr;
	node->value = value;
	node->index = heap->size;

	heap->tree[node->index] = node;

	heap->size++;

	heap_up(heap, node->index);

	return node;
}

bool HEAP_is_empty(Heap *heap)
{
	return (heap->size == 0);
}

unsigned int HEAP_size(Heap *heap)
{
	return heap->size;
}

HeapNode *HEAP_top(Heap *heap)
{
	return heap->tree[0];
}

void *HEAP_popmin(Heap *heap)
{
	void *ptr = heap->tree[0]->ptr;

	assert(heap->size != 0);

	heap_node_free(heap, heap->tree[0]);

	if (--heap->size) {
		heap_swap(heap, 0, heap->size);
		heap_down(heap, 0);
	}

	return ptr;
}

void HEAP_remove(Heap *heap, HeapNode *node)
{
	unsigned int i = node->index;

	assert(heap->size != 0);

	while (i > 0) {
		unsigned int p = HEAP_PARENT(i);

		heap_swap(heap, p, i);
		i = p;
	}

	HEAP_popmin(heap);
}

double HEAP_node_value(HeapNode *node)
{
	return node->value;
}

void *HEAP_node_ptr(HeapNode *node)
{
	return node->ptr;
}
