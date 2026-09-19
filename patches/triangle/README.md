# triangle 1.6 patch for GCC >= 15 (Ubuntu 26.04)

`triangle 1.6` (Shewchuk, 2005) builds fine on Ubuntu 24.04 (GCC 13) but fails on
Ubuntu 26.04 (GCC 15) for a different reason than igraph: **GCC 15 defaults to
C23 (`-std=gnu23`), where `void f();` declares a function taking *no*
arguments, while in C17 and earlier an empty parameter list meant "arguments
unspecified".**

Triangle is written in the legacy K&R style. Its header declares

```c
#ifdef ANSI_DECLARATORS
void triangulate(char const *, struct triangulateio *, struct triangulateio *,
                 struct triangulateio *);
void trifree(VOID *memptr);
#else /* not ANSI_DECLARATORS */
void triangulate();
void trifree();
#endif
```

and the makefile compiles the library with `TRILIBDEFS = -DTRILIBRARY`, i.e.
**without** `ANSI_DECLARATORS`, so the empty-parameter-list branch is the one
used. Under C23 that turns into `void trifree(void)`, which then clashes with
its own definition and with all 9 call sites:

```
./triangle.c:1446:1: error: number of arguments doesn't match prototype
./triangle.h:289:6: error: prototype declaration
./triangle.c:3990:5: error: too many arguments to function 'trifree'; expected 0, have 1
... (9 call sites)
./triangle.c:15729:1: error: number of arguments doesn't match prototype   (triangulate)
```

## Fix

One line, added to the existing `makefile.diff`:

```
CSWITCHES = ... -DLINUX -std=gnu17 ...
```

`-std=gnu17` is exactly the default GCC 13 used on Ubuntu 24.04, so this
restores the previous behaviour without touching a single line of Triangle's
2005 source. The 114 `-Wold-style-definition` warnings are unchanged (warnings
only, as before).

## Why not `-DANSI_DECLARATORS`?

Defining it in `TRILIBDEFS` removes the 13 errors too and activates real
prototypes, but it also activates every other prototype in the file and then
trips on a genuine latent bug introduced by the `char const *` change in
`triangle.c.diff`:

```
triangle.c:15750:23: error: passing argument 2 of 'parsecommandline'
                    from incompatible pointer type
```

(`&triswitches` is `char const **` but `parsecommandline` takes `char **`.)
That would need an extra cast or signature change; `-std=gnu17` avoids the
churn. Note that MechSys' own code is unaffected either way: it includes
`triangle.h` itself with `#define ANSI_DECLARATORS` (see
`mechsys/lib/mesh/unstructured.h`).

## Bonus: `tricall` segfaults (pre-existing, not a GCC 15 issue)

`make` also builds the bundled `tricall` demo, which crashed with SIGSEGV both
before and after the GCC 15 fix (it is not used by MechSys, so nobody noticed).
The `triangle.c.diff` patch added a `triedgemarks` field to
`struct triangulateio` that `triangle.c` allocates lazily:

```c
if (*triedgemarks == (int *) NULL) { *triedgemarks = trimalloc(...); }
```

`tricall.c` initialises every other output pointer of `mid`/`out` to `NULL`,
but was never updated for the new field, so `writeelements()` saw stack garbage
and dereferenced it. `tricall.c.diff` adds the two missing initialisations:

```c
mid.triedgemarks = (int *) NULL;
out.triedgemarks = (int *) NULL;
```

MechSys' own code was always fine here — `mechsys/lib/mesh/unstructured.h`
already sets `Tio.triedgemarks = NULL;`. If you prefer a minimal patch set,
`tricall.c.diff` and its line in `do_patch.bash` can be dropped without
affecting compilation.

