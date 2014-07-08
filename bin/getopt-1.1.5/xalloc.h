/*
 * Copyright (C) 2010 Davidlohr Bueso <dave@gnu.org>
 *
 * This file may be redistributed under the terms of the
 * GNU Lesser General Public License.
 *
 * General memory allocation wrappers for malloc, realloc, calloc and strdup
 */

#ifndef UTIL_LINUX_XALLOC_H
#define UTIL_LINUX_XALLOC_H

#include <stdlib.h>
#include <string.h>

#ifndef XALLOC_EXIT_CODE
# define XALLOC_EXIT_CODE 3
#endif

static void *xmalloc(const size_t size)
{
        void *ret = malloc(size);

        if (!ret && size) {
                fprintf(stderr, "%s: error: cannot allocate %zu bytes", program_invocation_short_name, size);
		exit(XALLOC_EXIT_CODE);
	}
        return ret;
}

static void *xrealloc(void *ptr, const size_t size)
{
        void *ret = realloc(ptr, size);

        if (!ret && size)
                fprintf(stderr, "%s: Error: cannot allocate %zu bytes", program_invocation_short_name, size);
        return ret;
}

#endif
