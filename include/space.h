#if SPACE < 2
#undef SPACE
#endif

#if SPACE > 3
#undef SPACE
#endif

#ifndef SPACE

#ifdef _TWO_DIM
#define SPACE 2
#else
#define SPACE 3
#endif

#endif

