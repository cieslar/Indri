#pragma once
/* strlcpy based on OpenBSDs strlcpy */
#include <stdio.h>
#include <sys/types.h>

size_t strlcpy(char *, const char *, size_t);

/*
 * Copy src to string dst of size siz.  At most siz-1 characters
 * will be copied.  Always NUL terminates (unless siz == 0).
 * Returns strlen(src); if retval >= siz, truncation occurred.
 */

