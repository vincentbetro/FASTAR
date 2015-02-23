/*
** Support functions for writing binary FIELDVIEW unstructured data files.
*/

/* Include system stuff for I/O and string and exit functions. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Include the defines for the FV_* codes and wall_info values. */
#include "fv_reader_tags.h"


/* Don't change these - used by fv_encode_elem_header ! */
#define MAX_NUM_ELEM_FACES     6
#define BITS_PER_WALL  3
#define ELEM_TYPE_BIT_SHIFT    (MAX_NUM_ELEM_FACES*BITS_PER_WALL)


/*
** fv_encode_elem_header:  return an encoded binary element header
**
** Input:
**    elem_type:  integer element type as shown in fv_reader_tags.h
**    wall_info:  array of integer "wall" flags, one for each face of
**                the element.  The wall flags are used during streamline
**                calculation.  Currently, the only meaningful values are
**                A_WALL and NOT_A_WALL as shown in fv_reader_tags.h.
**                Streamlines are forced away from faces marked as
**                "A_WALL", by limiting velocity and position very near
**                the wall.
** Output:
**    Function return value is the encoded binary element header.
*/

#ifndef FV_h
#define FV_h
unsigned int fv_encode_elem_header (int elem_type, int wall_info[]);

/*
** fwrite_str80:  write out a string padded to 80 characters.
**
** Like fwrite, this returns the number of items written, which
** should be 80 if successful, and less than 80 if it failed.
*/
#ifdef __STDC__
size_t fwrite_str80 (char *str, FILE *fp);
#else
int fwrite_str80 (char *str, FILE *fp);
#endif
#endif
