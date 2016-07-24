/**
 * Simple memory chunking allocator.
 *
 * Defines need to be set:
 * - TPOOL_IMPL_PREFIX: Prefix to use for the API.
 * - TPOOL_ALLOC_TYPE: Struct type this pool handles.
 * - TPOOL_STRUCT: Name for pool struct name.

 */

/* check we're not building directly */
#if !defined(TPOOL_ALLOC_TYPE) || \
	!defined(TPOOL_STRUCT) || \
	!defined(TPOOL_IMPL_PREFIX)
#  error "This file can't be compiled directly, include in another source file"
#endif

#define _CONCAT_AUX(MACRO_ARG1, MACRO_ARG2) MACRO_ARG1 ## MACRO_ARG2
#define _CONCAT(MACRO_ARG1, MACRO_ARG2) _CONCAT_AUX(MACRO_ARG1, MACRO_ARG2)
#define _TPOOL_PREFIX(id) _CONCAT(TPOOL_IMPL_PREFIX, _##id)

/* local identifiers */
#define pool_create		_TPOOL_PREFIX(pool_create)
#define pool_destroy	_TPOOL_PREFIX(pool_destroy)
#define pool_clear		_TPOOL_PREFIX(pool_clear)

#define pool_elem_free		_TPOOL_PREFIX(pool_elem_free)
#define pool_elem_alloc		_TPOOL_PREFIX(pool_elem_alloc)

#define pool_alloc_chunk	_TPOOL_PREFIX(pool_alloc_chunk)

/* private identifiers (only for this file, undefine after) */
#define TPoolChunk			_TPOOL_PREFIX(TPoolChunk)
#define TPoolChunkElemFree	_TPOOL_PREFIX(TPoolChunkElemFree)


struct TPoolChunk {
	struct TPoolChunk *prev;
	unsigned int    size;
	unsigned int    bufsize;
	struct TPOOL_ALLOC_TYPE buf[0];
};

struct TPoolChunkElemFree {
	struct TPoolChunkElemFree *next;
};

struct TPOOL_STRUCT {
	/* Always keep at least one chunk (never NULL) */
	struct TPoolChunk *chunk;
	/* when NULL, allocate a new chunk */
	struct TPoolChunkElemFree *free;
};

/**
 * Number of elems to include per #TPoolChunk when no reserved size is passed,
 * or we allocate past the reserved number.
 *
 * \note Optimize number for 64kb allocs.
 */
#define _TPOOL_CHUNK_DEFAULT_NUM \
	(((1 << 16) - sizeof(struct TPoolChunk)) / sizeof(TPOOL_ALLOC_TYPE))


/** \name Internal Memory Management
 * \{ */

static struct TPoolChunk *pool_alloc_chunk(
        unsigned int tot_elems, struct TPoolChunk *chunk_prev)
{
	struct TPoolChunk *chunk = malloc(
	        sizeof(struct TPoolChunk) + (sizeof(TPOOL_ALLOC_TYPE) * tot_elems));
	chunk->prev = chunk_prev;
	chunk->bufsize = tot_elems;
	chunk->size = 0;
	return chunk;
}

static struct TPOOL_ALLOC_TYPE *pool_elem_alloc(struct TPOOL_STRUCT *pool)
{
	TPOOL_ALLOC_TYPE *elem;

	if (pool->free) {
		elem = (TPOOL_ALLOC_TYPE *)pool->free;
		pool->free = pool->free->next;
	}
	else {
		struct TPoolChunk *chunk = pool->chunk;
		if (UNLIKELY(chunk->size == chunk->bufsize)) {
			chunk = pool->chunk = pool_alloc_chunk(_TPOOL_CHUNK_DEFAULT_NUM, chunk);
		}
		elem = &chunk->buf[chunk->size++];
	}

	return elem;
}

static void pool_elem_free(struct TPOOL_STRUCT *pool, TPOOL_ALLOC_TYPE *elem)
{
	struct TPoolChunkElemFree *elem_free = (struct TPoolChunkElemFree *)elem;
	elem_free->next = pool->free;
	pool->free = elem_free;
}

static void pool_create(struct TPOOL_STRUCT *pool, unsigned int tot_reserve)
{
	pool->chunk = pool_alloc_chunk((tot_reserve > 1) ? tot_reserve : _TPOOL_CHUNK_DEFAULT_NUM, NULL);
	pool->free = NULL;
}

static void pool_clear(struct TPOOL_STRUCT *pool)
{
	/* Remove all except the last chunk */
	while (pool->chunk->prev) {
		struct TPoolChunk *chunk_prev = pool->chunk->prev;
		free(pool->chunk);
		pool->chunk = chunk_prev;
	}
	pool->chunk->size = 0;
	pool->free = NULL;
}

static void pool_destroy(struct TPOOL_STRUCT *pool)
{
	struct TPoolChunk *chunk = pool->chunk;
	do {
		struct TPoolChunk *chunk_prev;
		chunk_prev = chunk->prev;
		free(chunk);
		chunk = chunk_prev;
	} while (chunk);

	pool->chunk = NULL;
	pool->free = NULL;
}

/** \} */

#undef _TPOOL_CHUNK_DEFAULT_NUM
#undef _CONCAT_AUX
#undef _CONCAT
#undef _TPOOL_PREFIX

#undef TPoolChunk
#undef TPoolChunkElemFree