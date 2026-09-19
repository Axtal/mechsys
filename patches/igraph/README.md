# igraph-0.8.2 patch for GCC >= 14 (Ubuntu 26.04)

igraph 0.8.2 (released 2020) builds fine on Ubuntu 24.04 (GCC 13), but fails on
Ubuntu 26.04 (GCC 15) because GCC 14 turned implicit function declarations into
hard errors. Two spots in igraph 0.8.2 rely on that old leniency:

| File | Problem |
|------|---------|
| `src/community_leiden.c` | calls `igraph_i_vector_binsearch_slice()`, which is defined via `FUNCTION(igraph_i_vector, binsearch_slice)` in `src/vector.pmt` (external linkage) but is **declared in no header** |
| `src/f2c/uninit.c` | `#define _GNU_SOURCE 1` sits at line 257, *after* `<stdio.h>`/`<stdlib.h>` were included, so `features.h` has already run and glibc never declares `feenableexcept()`/`fedisableexcept()` |

Errors seen without the patch:

```
community_leiden.c:428:17: error: implicit declaration of function 'igraph_i_vector_binsearch_slice'
f2c/uninit.c:264:9: error: implicit declaration of function 'fedisableexcept'
f2c/uninit.c:266:9: error: implicit declaration of function 'feenableexcept'
```

`do_patch.bash` is called by `mechsys/scripts/install_compile_deps.bash` right
after the tarball is unpacked (the igraph case there sets `DO_PATCH=1`). The
patch only adds a prototype and moves a `#define`; it is a no-op semantically
and is harmless on GCC 13.

Note: `include/igraph_threading.h` needs no patch — it is generated from
`include/igraph_threading.h.in` with `@HAVE_TLS@`, which `./configure
--enable-tls` (already used by the installer) sets to 1.
