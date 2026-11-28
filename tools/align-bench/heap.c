/* Injected into the heap benchmark child, never linked into the timed bench. */
#if defined(__APPLE__) && !defined(_DARWIN_C_SOURCE)
#define _DARWIN_C_SOURCE
#elif defined(__linux__) && !defined(_GNU_SOURCE)
#define _GNU_SOURCE
#endif

#include <errno.h>
#include <pthread.h>
#include <stdatomic.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define EXPORT __attribute__((visibility("default")))

#if defined(__APPLE__) && defined(__MACH__)
#include <malloc/malloc.h>
#define HOOK static
#define REPLACE(name) bf_##name
/* dyld does not substitute references in the image that provides the tuple. */
#define REAL(name) name
#define USABLE(ptr) malloc_size(ptr)
#define REALLOC_ZERO_FREES 0
#elif defined(__linux__) && defined(__GLIBC__)
#include <malloc.h>
#define HOOK EXPORT
#define REPLACE(name) name
/* Avoid symbol lookup in a malloc hook, including during loader startup. */
#define REAL(name) __libc_##name
#define USABLE(ptr) malloc_usable_size(ptr)
#define REALLOC_ZERO_FREES 1
extern void *__libc_malloc(size_t size);
extern void *__libc_calloc(size_t count, size_t size);
extern void *__libc_realloc(void *ptr, size_t size);
extern void __libc_free(void *ptr);
extern void *__libc_memalign(size_t alignment, size_t size);
#else
#error "Bitfront heap probing supports macOS and Linux with glibc"
#endif

/* A scope owns all blocks it frees or reallocates, and drops them before end. */
struct scope {
    size_t live;
    size_t peak;
    bool in_allocator;
};

static pthread_once_t key_once = PTHREAD_ONCE_INIT;
static pthread_key_t key;
static atomic_bool key_ready;
static _Thread_local struct scope local_scope;

static _Noreturn void fail(const char *message) {
    const char *parts[] = { "bitfront heap probe: ", message, "\n" };
    for (size_t i = 0; i < sizeof(parts) / sizeof(*parts); ++i) {
        if (write(STDERR_FILENO, parts[i], strlen(parts[i])) < 0) break;
    }
    abort();
}

static void init_key(void) {
    if (pthread_key_create(&key, NULL) != 0) {
        fail("cannot create allocation scope key");
    }
    atomic_store_explicit(&key_ready, true, memory_order_release);
}

static struct scope *active_scope(void) {
    // On Darwin, first access to local_scope may itself allocate. Hooks only
    // consult the pthread key, after it has been created outside any scope.
    if (!atomic_load_explicit(&key_ready, memory_order_acquire)) {
        return NULL;
    }
    return pthread_getspecific(key);
}

EXPORT void bf_heap_begin(void) {
    if (pthread_once(&key_once, init_key) != 0) {
        fail("cannot initialize allocation scope key");
    }
    if (active_scope() != NULL) {
        fail("nested allocation scope");
    }
    struct scope *scope = &local_scope;
    *scope = (struct scope){0};
    if (pthread_setspecific(key, scope) != 0) {
        fail("cannot start allocation scope");
    }
}

EXPORT size_t bf_heap_live(void) {
    struct scope *scope = active_scope();
    return scope == NULL ? 0 : scope->live;
}

EXPORT size_t bf_heap_end(void) {
    struct scope *scope = active_scope();
    if (scope == NULL || scope->in_allocator) {
        fail("no allocation scope to end");
    }
    scope->in_allocator = true;
    if (pthread_setspecific(key, NULL) != 0) {
        fail("cannot end allocation scope");
    }
    if (scope->live != 0) {
        fail("allocation scope ended with live heap");
    }
    return scope->peak;
}

static struct scope *enter_allocator(void) {
    struct scope *scope = active_scope();
    if (scope == NULL || scope->in_allocator) {
        return NULL;
    }
    scope->in_allocator = true;
    return scope;
}

static size_t usable_bytes(void *ptr) {
    if (ptr == NULL) {
        return 0;
    }
    int saved_errno = errno;
    size_t size = USABLE(ptr);
    errno = saved_errno;
    if (size == 0) {
        fail("allocator returned no usable size for a non-null block");
    }
    return size;
}

static size_t release_bytes(const struct scope *scope, void *ptr) {
    size_t bytes = usable_bytes(ptr);
    if (bytes > scope->live) {
        fail("free or realloc exceeds the allocation scope's live heap");
    }
    return bytes;
}

static void *finish_allocation(struct scope *scope, void *ptr) {
    if (scope != NULL) {
        size_t bytes = usable_bytes(ptr);
        if (bytes > SIZE_MAX - scope->live) {
            fail("allocation scope's live heap overflowed");
        }
        scope->live += bytes;
        if (scope->live > scope->peak) {
            scope->peak = scope->live;
        }
        scope->in_allocator = false;
    }
    return ptr;
}

HOOK void *REPLACE(malloc)(size_t size) {
    struct scope *scope = enter_allocator();
    return finish_allocation(scope, REAL(malloc)(size));
}

HOOK void *REPLACE(calloc)(size_t count, size_t size) {
    struct scope *scope = enter_allocator();
    return finish_allocation(scope, REAL(calloc)(count, size));
}

HOOK void *REPLACE(realloc)(void *ptr, size_t size) {
    struct scope *scope = enter_allocator();
    size_t old = scope == NULL ? 0 : release_bytes(scope, ptr);
    void *result = REAL(realloc)(ptr, size);
    // A failed realloc keeps the old block, except glibc's realloc(p, 0).
    // Only the replacement is live once realloc returns, even if it moved.
    if (scope != NULL && (result != NULL || (REALLOC_ZERO_FREES && size == 0))) {
        scope->live -= old;
    }
    return finish_allocation(scope, result);
}

HOOK void REPLACE(free)(void *ptr) {
    struct scope *scope = enter_allocator();
    if (scope != NULL) {
        scope->live -= release_bytes(scope, ptr);
    }
    REAL(free)(ptr);
    if (scope != NULL) {
        scope->in_allocator = false;
    }
}

static int real_posix_memalign(void **out, size_t alignment, size_t size) {
#if defined(__APPLE__)
    return posix_memalign(out, alignment, size);
#else
    if (alignment == 0 || alignment % sizeof(void *) != 0 ||
        (alignment & (alignment - 1)) != 0) {
        return EINVAL;
    }
    int saved_errno = errno;
    void *ptr = __libc_memalign(alignment, size);
    errno = saved_errno;
    if (ptr == NULL) {
        return ENOMEM;
    }
    *out = ptr;
    return 0;
#endif
}

HOOK int REPLACE(posix_memalign)(void **out, size_t alignment, size_t size) {
    struct scope *scope = enter_allocator();
    int status = real_posix_memalign(out, alignment, size);
    finish_allocation(scope, status == 0 ? *out : NULL);
    return status;
}

#if defined(__APPLE__)
#define INTERPOSE(name) \
    __attribute__((used, section("__DATA,__interpose"))) \
    static const struct { const void *replacement; const void *original; } \
    interpose_##name = { \
        (const void *)(uintptr_t)&bf_##name, (const void *)(uintptr_t)&name \
    }

INTERPOSE(malloc);
INTERPOSE(calloc);
INTERPOSE(realloc);
INTERPOSE(free);
INTERPOSE(posix_memalign);
#endif
