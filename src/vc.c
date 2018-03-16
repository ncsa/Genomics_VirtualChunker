#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>
#include <ctype.h>


/* #defines */

#define FALSE 0
#define TRUE  1
#define VC_VARS_MAGIC		0x1234
#define MAX_PATH_LEN		1024
#define MAX_SEQUENCE_LEN	10240
#define IO_BUF_SIZE		((size_t)(1024 * 1024))

/* Sequence line possible ID flags */
#define AT_LINE		0x01	/* line 1 */
#define SEQ_LINE	0x02	/* line 2 */
#define PLUS_LINE	0x04	/* line 3 */
#define QUAL_LINE	0x08	/* line 4 */

#define IS_SEQ_CHAR(c) ( ( (c) == 'A' ) || ( (c) == 'a' ) || \
                         ( (c) == 'C' ) || ( (c) == 'c' ) || \
                         ( (c) == 'G' ) || ( (c) == 'g' ) || \
                         ( (c) == 'N' ) || ( (c) == 'n' ) || \
                         ( (c) == 'T' ) || ( (c) == 't' ) )

#define IS_QUAL_CHAR(c) ( ( (int)(c) >= 33 ) && ( (int)(c) <= 126 ) )


/* type definitions */

struct vc_vars
{
   unsigned magic;	/* must be set to VC_VARS_MAGIC */

   /* parameters */
   int dry_run;
   unsigned long long chunk_len;
   unsigned long long chunk_index;
   unsigned long long max_chunk_index;
   char input_file_name[MAX_PATH_LEN +1];
   char output_file_name[MAX_PATH_LEN +1];
   char output_fifo_name[MAX_PATH_LEN +1];
   int verbosity;	/* verbosity == 0 ==> output on error only */
                        /* verbosity == 1 ==> terse functional output */
                        /* verbosity == 2 ==> verbose functional output */
                        /* verbosity >= 3 ==> debugging output */

   /* variables */
   off_t input_file_len;
   off_t nominal_chunk_start;
   off_t nominal_chunk_end;
   off_t actual_chunk_start;
   off_t actual_chunk_end;
   int input_fd;
   int output_fd;

   /* variables supporting dry runs / debugging */
   char * first_lines[4];
   char * last_lines[4];
   off_t first_lines_offsets[4];
   off_t last_lines_offsets[4];
};

struct line_header
{
   char * line_ptr;
   off_t offset;
   int line_len;
   unsigned id_flags;
};

struct seq_header
{
   off_t offsets[4];
   char * lines[4];
   int line_len[4];
   unsigned flags[4];
};

struct dry_run_chunk_record
{
   unsigned long long chunk_index;

   off_t nominal_chunk_start;
   off_t nominal_chunk_end;
   off_t actual_chunk_start;
   off_t actual_chunk_end;

   /* the following fields are only used if vars->verbosity > 1 */
   char * first_lines[4];
   char * last_lines[4];
   off_t first_lines_offsets[4];
   off_t last_lines_offsets[4];
};


/* function declarations */

static void close_input_file(struct vc_vars * vars_ptr);
static int compute_actual_chunk_end_point(struct vc_vars * vars_ptr);
static int compute_actual_chunk_start_point(struct vc_vars * vars_ptr);
static int compute_nominal_chunk_endpoints(struct vc_vars * vars_ptr);
static void count_assignment_errors(int offset, struct line_header line_map[], 
                                    int line_count, int * err_count_ptr);
static int dry_run(struct vc_vars * vars_ptr);
static void dry_run__dump_chunk_boundaries(
			              struct dry_run_chunk_record * drc_array,
                                      struct vc_vars * vars_ptr);
static void dump_params(struct vc_vars * vars_ptr);
static void dump_line_map(struct line_header line_map[], int line_count);
static void dump_seq_map(struct seq_header map[], int seq_count);
static int get_input_file_len(struct vc_vars * vars_ptr);
static int get_params(int argc, char *argv[], struct vc_vars * vars_ptr);
static int load_buf_from_file(off_t offset, off_t buf_len, char *buf_ptr,
                              struct vc_vars * vars_ptr);
static int map_test_buf(off_t test_buf_offset, char test_buf[], 
                        int test_buf_len, struct seq_header map[], int map_len,
                        int * seq_count_ptr, struct vc_vars * vars_ptr);
static int map_test_buf_lines(off_t test_buf_offset, char test_buf[], 
                              int test_buf_len, struct line_header line_map[],
                              int map_len, int *line_count_ptr,
                              struct vc_vars * vars_ptr);
static void post_err_mssg(char * err_msg);
static void post_syntax_err(char * err_msg);
static int setup_output_fd(struct vc_vars * vars_ptr);
static void takedown_output_fd(struct vc_vars * vars_ptr);
static void usage(void);
static int write_buf_to_output_fd(off_t buf_len, char *buf_ptr, 
                                  struct vc_vars * vars_ptr);
static int write_chunk(struct vc_vars * vars_ptr);


/*-------------------------------------------------------------------------
 * Function:    close_input_file()
 *
 * Purpose:     Attempt to close the input file.  Generate an error 
 *		error message on failure.
 *
 *		Note that we use file descriptor based I/O here because
 *		we are working with large files, and thus may have to 
 *		seek to locations that cannot be represented with a long.
 * 
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static void
close_input_file(struct vc_vars * vars_ptr)
{
   int result;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert(strlen(vars_ptr->input_file_name) > 0);
   assert(vars_ptr->input_fd != -1);

   errno = 0;
   result = close(vars_ptr->input_fd);
   saved_errno = errno;

   if ( result == -1 )
   {
      fprintf(stderr, "\nCant't close input file \"%s\".\n",
              vars_ptr->input_file_name);
      fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
              strerror(saved_errno));
   }
   else
   {
      vars_ptr->input_fd = -1;
   }

   return;

} /* close_input_file() */


/*-------------------------------------------------------------------------
 * Function:    compute_actual_chunk_end_point()
 *
 * Purpose:     Compute the actual offset of the chunk end point.  This
 *	        is simply the offset of the end of the last complete 
 *		FASTQ sequence that begins within the nominal chunk.
 *
 *		If this is the last chunk, (i.e. chunk_index == max index), 
 *              life is easy, as the noninal end point and the actual end
 *		point are the same.
 *
 *		In all other cases, we load a buffer from the file 
 *		centered on the nominal end point, map all complete
 *		sequences in it, and set our actual start point equal 
 *		to one less than the offset in the input file of the 
 *		earliest sequence whose first character is at offset 
 *		greater than the noninal chunk end.  
 *
 * Return:      TRUE is successful, and FALSE if an error is detected.
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */
static int
compute_actual_chunk_end_point(struct vc_vars * vars_ptr)
{
   int proceed = TRUE;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert(vars_ptr->input_fd >= 0);

   if ( vars_ptr->chunk_index == vars_ptr->max_chunk_index )
   {
      assert(vars_ptr->nominal_chunk_end == vars_ptr->input_file_len - 1);
      vars_ptr->actual_chunk_end = (off_t)(vars_ptr->input_file_len - 1);
   }
   else /* we have work to do */
   {
      char test_buf[4 * MAX_SEQUENCE_LEN];
      int i;
      int seq_count;
      off_t test_buf_start;
      off_t test_buf_len;
      struct seq_header map[MAX_SEQUENCE_LEN / 4];

      test_buf_start = vars_ptr->nominal_chunk_end;

      assert(test_buf_start > (off_t)(MAX_SEQUENCE_LEN *2));

      test_buf_start -= (off_t)(2 * MAX_SEQUENCE_LEN);

      assert(test_buf_start < vars_ptr->input_file_len);

      test_buf_len = (off_t)(4 * MAX_SEQUENCE_LEN);

      assert(test_buf_start + test_buf_len < vars_ptr->input_file_len);

      proceed = load_buf_from_file(test_buf_start, test_buf_len, 
                                   test_buf, vars_ptr);

      if ( proceed )
         proceed = map_test_buf(test_buf_start, test_buf, test_buf_len, map, 
                                MAX_SEQUENCE_LEN / 4, &seq_count, vars_ptr);

      if ( proceed )
      {
         int i = 0;

         assert(seq_count < (MAX_SEQUENCE_LEN / 4));

         while ( ( i < seq_count ) && 
                 ( map[i].offsets[0] <= vars_ptr->nominal_chunk_end ) )
         {
            i++;
         }

         i--; /* to get back to the last sequence that starts before 
               * the nominal end of chunk.
               */

         /* if this assert fails, MAX_SEQUENCE_LEN is probably too small */
         assert(i > 0);

         /* make sure we got the dividing line correctly */
         assert( map[i].offsets[0] <= vars_ptr->nominal_chunk_end );
         assert( map[i+1].offsets[0] > vars_ptr->nominal_chunk_end );

         if ( i >= seq_count )
         {
            proceed = FALSE;

            fprintf(stderr, "\n%s %s\n\n",
                    "ERROR: Can't find last sequence that starts at or ",
                    "before the nominal chunk end.");

            if ( vars_ptr->verbosity >= 3 )
            {
               dump_seq_map(map, seq_count);
            }
         } 
         else /* we found our actual start point */
         {
            assert(i + 1 < seq_count);

            vars_ptr->actual_chunk_end = map[i+1].offsets[0] - 1;

            assert(vars_ptr->actual_chunk_end >= 
                   vars_ptr->nominal_chunk_end);
            assert(vars_ptr->actual_chunk_end < 
                   (vars_ptr->nominal_chunk_end + 
                    (off_t)(2 * MAX_SEQUENCE_LEN)));

            if ( ( vars_ptr->dry_run ) && ( vars_ptr->verbosity >= 2 ) )
            {
               int j;

               assert(vars_ptr->chunk_index < vars_ptr->max_chunk_index);

               for ( j = 0; j < 4; j++ )
               {
                  vars_ptr->last_lines[j] = strdup(map[i].lines[j]);
                  vars_ptr->last_lines_offsets[j] = map[i].offsets[j];
               }
            }
            else if ( ( ! vars_ptr->dry_run ) && ( vars_ptr->verbosity >= 1 ) )
            {
               if ( vars_ptr->verbosity >= 2 )
               {
                  fprintf(stdout, 
                          "==========================================\n");
                  fprintf(stdout, 
                          "Sequence before actual chunk end offset = %llu\n",
                          (unsigned long long)(map[i].offsets[0]));
                  fprintf(stdout, "%s\n%s\n%s\n%s\n",
                          map[i].lines[0], map[i].lines[1],
                          map[i].lines[2], map[i].lines[3]);
                  fprintf(stdout, 
                          "==========================================\n");
               }
           
               fprintf(stdout, "Actual chunk end offset = %llu.\n",
                       (unsigned long long)(vars_ptr->actual_chunk_end));

               if ( vars_ptr->verbosity >= 2 )
               {
                  fprintf(stdout, 
                          "==========================================\n");
                  fprintf(stdout, 
                          "Sequence after actual chunk end offset = %llu\n",
                          (unsigned long long)(map[i + 1].offsets[0]));
                  fprintf(stdout, "%s\n%s\n%s\n%s\n",
                          map[i + 1].lines[0], map[i + 1].lines[1],
                          map[i + 1].lines[2], map[i + 1].lines[3]);
                  fprintf(stdout, 
                          "==========================================\n");
               }
            }
         }
      }
   }

   return(proceed);

} /* compute_actual_chunk_end_point() */


/*-------------------------------------------------------------------------
 * Function:     compute_actual_chunk_start_point()
 *
 * Purpose:      Compute the actual offset of the chunk start point.  This
 *		 is simply the offset of the first complete FASTQ sequence
 *		 that begins on or after the nominal chunk start point.
 *
 *		 If this is the first chunk, (i.e. chunk_index == 0), life
 *		 is easy, as the noninal start point and the actual start
 *		 point are the same.
 *
 *		 In all other cases, we load a buffer from the file 
 *		 centered on the nominal start point, map all complete
 *		 sequences in it, and set our actual start point equal 
 *		 to the offset in the input file of the earliest sequence
 *		 whose first character is at offset greater than or equal 
 *		 to the noninal chunk start.  
 *
 * Return:       TRUE is successful, and FALSE if an error is detected.
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static int
compute_actual_chunk_start_point(struct vc_vars * vars_ptr)
{
   int proceed = TRUE;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert(vars_ptr->input_fd >= 0);

   if ( vars_ptr->chunk_index == 0 )
   {
      assert(vars_ptr->nominal_chunk_start == 0);
      vars_ptr->actual_chunk_start = (off_t)0;
   }
   else /* we have work to do */
   {
      char test_buf[4 * MAX_SEQUENCE_LEN];
      int i;
      int seq_count;
      off_t test_buf_start;
      off_t test_buf_len;
      struct seq_header map[MAX_SEQUENCE_LEN / 4];

      test_buf_start = vars_ptr->nominal_chunk_start;

      assert(test_buf_start > (off_t)(MAX_SEQUENCE_LEN *2));

      test_buf_start -= (off_t)(2 * MAX_SEQUENCE_LEN);

      assert(test_buf_start < vars_ptr->input_file_len);

      test_buf_len = (off_t)(4 * MAX_SEQUENCE_LEN);

      assert(test_buf_start + test_buf_len < vars_ptr->input_file_len);

      proceed = load_buf_from_file(test_buf_start, test_buf_len, 
                                   test_buf, vars_ptr);

      if ( proceed )
         proceed = map_test_buf(test_buf_start, test_buf, test_buf_len, map, 
                                MAX_SEQUENCE_LEN / 4, &seq_count, vars_ptr);

      if ( proceed )
      {
         int i = 0;

         assert(seq_count < (MAX_SEQUENCE_LEN / 4));

         while ( ( i < seq_count ) && 
                 ( map[i].offsets[0] < vars_ptr->nominal_chunk_start ) )
         {
            i++;
         }

         /* if this assert fails, MAX_SEQUENCE_LEN is probably too small */
         assert(i > 0);

         /* make sure we got the dividing line correctly */
         assert( map[i].offsets[0] >= vars_ptr->nominal_chunk_start );
         assert( map[i-1].offsets[0] < vars_ptr->nominal_chunk_start );

         if ( i >= seq_count )
         {
            proceed = FALSE;

            fprintf(stderr, "\n%s %s\n\n",
                    "ERROR: Can't find sequence that starts on or ",
                    "after the nominal chunk start.");

            if ( vars_ptr->verbosity >= 3 )
            {
               dump_seq_map(map, seq_count);
            }
         } 
         else /* we found our actual start point */
         {
            vars_ptr->actual_chunk_start = map[i].offsets[0];

            assert(vars_ptr->actual_chunk_start >= 
                   vars_ptr->nominal_chunk_start);
            assert(vars_ptr->actual_chunk_start < 
                   (vars_ptr->nominal_chunk_start + 
                    (off_t)(2 * MAX_SEQUENCE_LEN)));

            if ( ( vars_ptr->dry_run ) && ( vars_ptr->verbosity >= 2 ) )
            {
               int j;

               assert(vars_ptr->chunk_index > 0 );

               for ( j = 0; j < 4; j++ )
               {
                  vars_ptr->first_lines[j] = strdup(map[i].lines[j]);
                  vars_ptr->first_lines_offsets[j] = map[i].offsets[j];
               }
            }
            else if ( ( ! vars_ptr->dry_run ) && ( vars_ptr->verbosity >= 1 ) )
            {
               if ( vars_ptr->verbosity >= 2 )
               {
                  fprintf(stdout, 
                          "==========================================\n");
                  fprintf(stdout, 
                          "Sequence before actual chunk start offset = %llu\n",
                           (unsigned long long)(map[i - 1].offsets[0]));
                  fprintf(stdout, "%s\n%s\n%s\n%s\n",
                          map[i - 1].lines[0], map[i - 1].lines[1],
                          map[i - 1].lines[2], map[i - 1].lines[3]);
                  fprintf(stdout, 
                          "==========================================\n");
               }
           
               fprintf(stdout, "Actual chunk start offset = %llu.\n",
                       (unsigned long long)(vars_ptr->actual_chunk_start));

               if ( vars_ptr->verbosity >= 2 )
               {
                  fprintf(stdout, 
                          "==========================================\n");
                  fprintf(stdout, 
                          "Sequence after actual chunk start offset = %llu\n",
                          (unsigned long long)(map[i].offsets[0]));
                  fprintf(stdout, "%s\n%s\n%s\n%s\n",
                          map[i].lines[0], map[i].lines[1],
                          map[i].lines[2], map[i].lines[3]);
                  fprintf(stdout, 
                          "==========================================\n");
               }
            }
         }
      }
   }

   return(proceed);

} /* compute_actual_chunk_start_point() */


/*-------------------------------------------------------------------------
 * Function:     compute_nominal_chunk_endpoints()
 *
 * Purpose:      Using the input file length, the chunk length, the 
 *               chunk index, and max index, compute the inclusive 
 *               start and end points for the chunk to be produced.
 *
 *               In general, these are:
 *
 *		    chunk_len * chunk_index
 *
 *		 for the start point, and 
 *
 *		    (chunk_len * (chunk_index +1)) - 1
 *
 *		 for the end point.  However, if chunk_index equals 
 *		 max_chunk_index, the end point is input_file_len - 1.
 *
 *		 Further, if (input_file_len - 1) is less than the noninal
 *               end point of the chunk, set the end point to 
 *               (input_file_len - 1).
 *
 *		 Finally, if (input_file_len - 1) is less than the nominal 
 *		 start point of the chunk, post an error message, and 
 *               return FALSE.
 *
 * Return:       TRUE is successful, and FALSE if an error is detected.
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static int
compute_nominal_chunk_endpoints(struct vc_vars * vars_ptr)
{
   int proceed = TRUE;
   off_t chunk_start;
   off_t chunk_end;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);

   chunk_start = (off_t)(vars_ptr->chunk_index * vars_ptr->chunk_len);
   chunk_end = 
      (off_t)(((vars_ptr->chunk_index + 1) * vars_ptr->chunk_len) - 1);

   assert(chunk_start < chunk_end);

   if ( chunk_start > vars_ptr->input_file_len - 1 ) 
   {
      post_err_mssg("chunk is empty.");
      proceed = FALSE;
   }
   else 
   {
      if ( ( vars_ptr->chunk_index == vars_ptr->max_chunk_index ) ||
           ( chunk_end >= vars_ptr->input_file_len ) )
      {
         chunk_end = vars_ptr->input_file_len - 1;
      }

      assert(chunk_start <= chunk_end);

      vars_ptr->nominal_chunk_start = chunk_start;
      vars_ptr->nominal_chunk_end = chunk_end;

      if ( vars_ptr->verbosity >= 3 )
      {
         off_t nominal_chunk_len;

         nominal_chunk_len = vars_ptr->nominal_chunk_end - 
                             vars_ptr->nominal_chunk_start + 1;

         fprintf(stdout, "nominal chunk start point = %lld\n",
                 (unsigned long long)(vars_ptr->nominal_chunk_start));
         fprintf(stdout, "nominal chunk end point = %lld\n",
                 (unsigned long long)(vars_ptr->nominal_chunk_end));
         fprintf(stdout, "nominal chunk length = %lld\n",
                 (unsigned long long)(nominal_chunk_len));
      }

      if ( (chunk_end - chunk_start + 1) <= (4 * MAX_SEQUENCE_LEN) )
      {
         proceed = FALSE;

         fprintf(stderr, 
                 "\nERROR: Nominal chunk length may not be less than %llu.\n",
                 (unsigned long long)(4 * MAX_SEQUENCE_LEN));
         fprintf(stderr, "      chunk index = %lld\n", vars_ptr->chunk_index);
         fprintf(stderr, "	nominal chunk start point = %lld\n",
                 (unsigned long long)(vars_ptr->nominal_chunk_start));
         fprintf(stderr, "	nominal chunk end point = %lld\n",
                 (unsigned long long)(vars_ptr->nominal_chunk_end));
         fprintf(stderr, "	nominal chunk length = %lld\n",
                 (unsigned long long)(chunk_end - chunk_start + 1));
         fprintf(stderr, 
            "Please revise chunk size and/or max chunk index and rerun.\n\n");
      }
   }

   return(proceed);

} /* compute_nominal_chunk_endpoints() */


/*-------------------------------------------------------------------------
 * Function:    count_assignment_errors()
 *
 * Purpose:     Each sequence must consist of four lines -- the "@" line,
 *		the raw sequence line, the "+" plus line, and the quality
 *		line.  Assuming the given offset for the first "@" line 
 *		in the line map, count the number of contradictions between
 *		the assumed offset and the id flags set for each line.
 *		
 *		Return this value in *err_count_ptr.
 *
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static void
count_assignment_errors(int offset, 
                        struct line_header line_map[], 
                        int line_count, 
                        int * err_count_ptr)
{
   int i;
   int line_type;
   int err_count = 0;
   /* map offset to expected type of first line in line map */
   int init_line_type[4] = {0, 3, 2, 1};

   assert((offset >= 0) && ( offset <= 3));
   assert(line_map != NULL);
   assert(line_count >= 16);
   assert(err_count_ptr != NULL);

   line_type = init_line_type[offset];

   for ( i = 0; i < line_count; i++ )
   {
      switch ( line_type )
      {
         case 0:
            if ( (line_map[i].id_flags & AT_LINE) == 0 )
            {
               err_count++;
            }
            line_type = 1;
            break;

         case 1:
            if ( (line_map[i].id_flags & SEQ_LINE) == 0 )
            {
               err_count++;
            }
            line_type = 2;
            break;

         case 2:
            if ( (line_map[i].id_flags & PLUS_LINE) == 0 )
            {
               err_count++;
            }
            line_type = 3;
            break;

         case 3:
            if ( (line_map[i].id_flags & QUAL_LINE) == 0 )
            {
               err_count++;
            }
            line_type = 0;
            break;

         default:
            /* this should be unreachable */
            assert(FALSE);
            break;
      }
   }

   assert( ( err_count >= 0 ) && ( err_count <= line_count ) );

   *err_count_ptr = err_count;

   return;

} /* count_assignment_errors() */


/*-------------------------------------------------------------------------
 * Function:    dry_run()
 *
 * Purpose:     Find the nominal and actual start and end points of 
 *		all chunks, and display these values if verbosity > 0.  
 *              If verbosity is greater than 1, also display the reads 
 *              just before and just after the end points.  
 *
 *		Note that this function is somewhat memory intensive, and 
 *		is included primarily to support testing -- both of the 
 *		code itself, and of chunks size and number selections.
 *
 *		As a result of these priorities, this function is a bit 
 *		inefficient, as the point is to calculate each of the 
 *		end points just the way they would be calculated for 
 *		each chunk in actual application.  Thus the endpoints 
 *		are calculated twice for each boundary between chunks.
 *
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */
static int
dry_run(struct vc_vars * vars_ptr)
{
   int i;
   int j;
   int proceed = TRUE;
   struct dry_run_chunk_record * drc_array = NULL;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);

   /* allocate and initialize the array of chunk records */
   drc_array = malloc((size_t)(vars_ptr->max_chunk_index + 1) *
                      sizeof(struct dry_run_chunk_record));

   if ( drc_array == NULL ) 
   {
      proceed = FALSE;
      post_err_mssg("malloc of drc_array failed.");
   }
   else
   {
      for ( i = 0; i <= vars_ptr->max_chunk_index; i++ )
      {
         drc_array[i].chunk_index = (unsigned long long)i;

         drc_array[i].nominal_chunk_start = 0;
         drc_array[i].nominal_chunk_end = 0;
         drc_array[i].actual_chunk_start = 0;
         drc_array[i].actual_chunk_end = 0;

         for ( j = 0; j < 4; j++ )
         {
            drc_array[i].first_lines[j] = NULL;
            drc_array[i].last_lines[j] = NULL;
            drc_array[i].first_lines_offsets[j] = 0;
            drc_array[i].last_lines_offsets[j] = 0;
         }
      }
   }

   /* compute chunk start and end for each chunk */
   i = 0;
   while ( ( proceed ) && ( i <= vars_ptr->max_chunk_index ) )
   {
      /* Initialize *vars_ptr for this chunk endpoint calculation. */
      vars_ptr->chunk_index = (unsigned long long)i;

      vars_ptr->nominal_chunk_start = (off_t)0;
      vars_ptr->nominal_chunk_end   = (off_t)0;
      vars_ptr->actual_chunk_start  = (off_t)0;
      vars_ptr->actual_chunk_end    = (off_t)0;

      for ( j = 0; j < 4; j++ )
      {
         vars_ptr->first_lines[j]         = NULL;
         vars_ptr->last_lines[j]          = NULL;
         vars_ptr->first_lines_offsets[j] = (off_t)0;
         vars_ptr->last_lines_offsets[j]  = (off_t)0;
      }

      proceed = compute_nominal_chunk_endpoints(vars_ptr);

      if ( proceed ) proceed = compute_actual_chunk_start_point(vars_ptr);

      if ( proceed ) proceed = compute_actual_chunk_end_point(vars_ptr);

      if ( proceed ) /* save the results for this chunk */
      {
         assert(drc_array[i].chunk_index == (unsigned long long)i);

         drc_array[i].nominal_chunk_start = vars_ptr->nominal_chunk_start;
         drc_array[i].nominal_chunk_end = vars_ptr->nominal_chunk_end;
         drc_array[i].actual_chunk_start = vars_ptr->actual_chunk_start;
         drc_array[i].actual_chunk_end = vars_ptr->actual_chunk_end;

         for ( j = 0; j < 4; j++ )
         {
            drc_array[i].first_lines[j] = vars_ptr->first_lines[j];
            drc_array[i].last_lines[j] = vars_ptr->last_lines[j];
            drc_array[i].first_lines_offsets[j] = 
                  vars_ptr->first_lines_offsets[j];
            drc_array[i].last_lines_offsets[j] = 
                  vars_ptr->last_lines_offsets[j];
         }
      }

      i++;
   }

   /* display start and end of each chunk */
   if ( ( proceed ) && ( vars_ptr->verbosity >= 1 ) )
      dry_run__dump_chunk_boundaries(drc_array, vars_ptr);

   /* free dynamically allocated memory */
   if ( drc_array != NULL )
   {
      for ( i = 0; i <= vars_ptr->max_chunk_index; i++ )
      {
         for ( j = 0; j < 4; j++ )
         {
            if ( drc_array[i].first_lines[j] != NULL )
            {
               free(drc_array[i].first_lines[j]);
               drc_array[i].first_lines[j] = NULL;
            }

            if ( drc_array[i].last_lines[j] != NULL )
            {
               free(drc_array[i].last_lines[j]);
               drc_array[i].last_lines[j] = NULL;
            }
         }
      }

      free(drc_array);
      drc_array = NULL;
   }

   return(proceed);

} /* dry_run() */


/*-------------------------------------------------------------------------
 * Function:    dry_run__dump_chunk_boundaries()
 *
 * Purpose:     Write the nominal and actual chunk boundaries to stdout
 *		in tabular form.
 *
 *		If the verbosity is high enough, also display the sequences
 *		just before and after each chunk boundary.
 *
 * Return:      void
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static void
dry_run__dump_chunk_boundaries(struct dry_run_chunk_record * drc_array,
                               struct vc_vars * vars_ptr)
{
   char *hdr_0   = "   Chunk             Chunk          Chunk         Chunk";
   char *hdr_1   = "   Index             Start          End           Length";
   char *divider = "   =====================================================";

   int i;
   int j;

   assert(drc_array != NULL);
   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);

   fprintf(stdout, "\nChunk Boundary Table for the file \"%s\":\n\n",
           vars_ptr->input_file_name);

   fprintf(stdout, "Input file length = %lld\n\n", 
           (unsigned long long)(vars_ptr->input_file_len));

   fprintf(stdout, "%s\n%s\n%s\n", hdr_0, hdr_1, divider);

   for ( i = 0; i <= vars_ptr->max_chunk_index; i++ )
   {
      assert( i == drc_array[i].chunk_index );
      fprintf(stdout, "    %3d       %12lld   %12lld   %12lld\n",
              i, (unsigned long long)(drc_array[i].actual_chunk_start),
              (unsigned long long)(drc_array[i].actual_chunk_end),
              (unsigned long long)(drc_array[i].actual_chunk_end -
                                   drc_array[i].actual_chunk_start + 1));
   }

   fprintf(stdout, "%s\n\n", divider);

   if ( vars_ptr->verbosity >= 2 ) 
   {
      char * hdr_2 = 
        "Displaying sequences just before and after each chunk boundary.";
      char * hdr_3 =
        "The integer at the beginning of each line is its offset in the file.";

      fprintf(stdout, "%s\n%s\n\n", hdr_2, hdr_3);

      for ( i = 0; i < vars_ptr->max_chunk_index; i++ )
      {
         assert( i     == drc_array[i].chunk_index );
         assert( i + 1 == drc_array[i+1].chunk_index );

         fprintf(stdout, "\n");

         for ( j = 0; j < 4; j++ )
         {
            fprintf(stdout, "%12lld  %s\n", 
                    (unsigned long long)drc_array[i].last_lines_offsets[j],
                    drc_array[i].last_lines[j]);
         }

         fprintf(stdout, 
                 "========= Chunk %d / %d Boundary (%lld/%lld) =========\n", 
                 i, i + 1, (unsigned long long)(drc_array[i].actual_chunk_end),                 (unsigned long long)(drc_array[i+1].actual_chunk_start));

         for ( j = 0; j < 4; j++ )
         {
            fprintf(stdout, "%12lld  %s\n", 
                    (unsigned long long)drc_array[i+1].first_lines_offsets[j],
                    drc_array[i+1].first_lines[j]);
         }

         fprintf(stdout, "\n");
      }
   } 

   return;

} /* dry_run__dump_chunk_boundaries() */


/*-------------------------------------------------------------------------
 * Function:     dump_params()
 *
 * Purpose:      Write the parameters to stdout.
 *
 * Return:       void
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static void
dump_params(struct vc_vars * vars_ptr)
{
   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);

   fprintf(stdout, "\n\nParameters:\n\n");
   fprintf(stdout, "   dry_run          = %d\n", vars_ptr->dry_run);
   fprintf(stdout, "   chunk_len        = %llu\n", vars_ptr->chunk_len);
   fprintf(stdout, "   chunk_index      = %llu\n", vars_ptr->chunk_index);
   fprintf(stdout, "   max_chunk_index  = %llu\n", vars_ptr->max_chunk_index);
   fprintf(stdout, "   input_file_name  = \"%s\"\n", 
           vars_ptr->input_file_name);
   fprintf(stdout, "   output_file_name = \"%s\"\n", 
           vars_ptr->output_file_name);
   fprintf(stdout, "   output_fifo_name = \"%s\"\n", 
           vars_ptr->output_fifo_name);
   fprintf(stdout, "   verbosity        = %d\n\n", vars_ptr->verbosity);

   return;

} /* dump_params() */


/*-------------------------------------------------------------------------
 * Function:    dump_line_map()
 *
 * Purpose:     Write the contents of the supplied line map to stdout.
 *
 * Return:      void
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */
static void
dump_line_map(struct line_header line_map[], int line_count)
{
   int i;

   assert(line_map != NULL);
   assert(line_count > 0);

   fprintf(stdout, "\n\nLine Map:\n\n");

   for ( i = 0; i < line_count; i++ )
   {
      fprintf(stdout, "Line %d at offset %llu:\n", i, 
              (unsigned long long)(line_map[i].offset));
      fprintf(stdout, "	len: %d, flags: 0x%x, body: \"%s\"\n",
              line_map[i].line_len, line_map[i].id_flags, 
              line_map[i].line_ptr);
   }

   return;

} /* dump_line_map() */



/*-------------------------------------------------------------------------
 * Function:    dump_seq_map()
 *
 * Purpose:     Write the contents of the supplied sequence map to 
 *		stdout.
 *
 * Return:      void
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */
static void
dump_seq_map(struct seq_header map[], int seq_count)
{
   int i;

   assert(map != NULL);
   assert(seq_count > 0);

   fprintf(stdout, "\n\nSequence Map:\n\n");

   for ( i = 0; i < seq_count; i++ )
   {
      fprintf(stdout, "Sequence %d at offset %llu:\n", i, 
              (unsigned long long)(map[i].offsets[0]));
      fprintf(stdout, "	line 1: len: %d, flags: 0x%x, body: \"%s\"\n",
              map[i].line_len[0], map[i].flags[0], map[i].lines[0]);
      fprintf(stdout, "	line 2: len: %d, flags: 0x%x, body: \"%s\"\n",
              map[i].line_len[1], map[i].flags[1], map[i].lines[1]);
      fprintf(stdout, "	line 3: len: %d, flags: 0x%x, body: \"%s\"\n",
              map[i].line_len[2], map[i].flags[2], map[i].lines[2]);
      fprintf(stdout, "	line 4: len: %d, flags: 0x%x, body: \"%s\"\n\n",
              map[i].line_len[3], map[i].flags[3], map[i].lines[3]);
   }

   return;

} /* dump_seq_map() */


/*-------------------------------------------------------------------------
 *
 * Function:    get_input_file_len()
 *
 * Purpose:     Stat the input file, load its size (in bytes) into the 
 *              input_file_len field of *vars_ptr, and return TRUE.
 *
 *              If the file does not exist, or if the function fails
 *              for any other reason, *file_len_ptr is undefined, and
 *              the function returns FALSE.
 *
 * Return:      Success: TRUE
 *              Failure: FALSE
 *
 * Changes:     None.
 *
 *-------------------------------------------------------------------------
 */

static int
get_input_file_len(struct vc_vars * vars_ptr)
{
   int proceed = TRUE;
   int saved_errno;
   struct stat buf;

   assert(vars_ptr != NULL);
   assert(strlen(vars_ptr->input_file_name) > 0);

   if ( proceed ) {

      if ( stat(vars_ptr->input_file_name, &buf) != 0 ) 
      {
         saved_errno = errno; /* save a copy of errno */

         proceed = FALSE;

         fprintf(stderr, "\nCant't stat input file \"%s\".\n",
                 vars_ptr->input_file_name);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      } else {

         vars_ptr->input_file_len = buf.st_size;

      }
   }

   if ( ( proceed ) && ( vars_ptr->verbosity >= 3 ) ) 
   {
      fprintf(stdout, "input file length = %llu\n", 
              (unsigned long long)(vars_ptr->input_file_len));
   }

   return(proceed);

} /* get_input_file_len() */


/*-------------------------------------------------------------------------
 * Function:     get_params()
 *
 * Purpose:      Load the user supplied parameters into the supplied
 *		 instance of struct _params.  
 *
 *		 If successful, return TRUE.
 *
 *		 On failure, generate an error message, and return FALSE.
 *
 * Return:       TRUE on success, FALSE otherwise.
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static int
get_params(int argc,
           char *argv[],
           struct vc_vars * vars_ptr)
{
   char *char_ptr;
   int have_chunk_size = FALSE;
   int have_chunk_index = FALSE;
   int have_max_chunk_index = FALSE;
   int opt;
   int proceed = TRUE; 
   unsigned long long multiplier;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);

   if ( argc < 5 ) 
   {
      post_syntax_err("too few arguments.");
      proceed = FALSE;
   }

   while ( ( proceed ) &&
           ( (opt = getopt(argc, argv, "ds:i:m:f:p:v?")) != -1 ) )
   {
      switch ( opt ) 
      {
         case 'd':
            vars_ptr->dry_run = TRUE;
            break;

         case 's':
	    char_ptr = optarg;
            assert(optarg != NULL);
            while ( ( *char_ptr != '\0' ) && ( isdigit(*char_ptr) ) )
            {
               char_ptr++;
            }

            if ( char_ptr == optarg ) 
            {
               post_syntax_err(
                  "Chunk size must be specified with an integer.");
               proceed = FALSE;
            }

            if ( ( proceed ) && ( *char_ptr != '\0' ) )
            {
               switch ( *char_ptr )
               {
                  case 'k':
                     multiplier = 1024;
                     break;

                  case 'm':
                     multiplier = 1024 * 1024;
                     break;

                  case 'g':
                     multiplier = 1024 * 1024 * 1024;
                     break;

                  default:
                     post_syntax_err(
                        "chunk size suffix must be either 'k', 'm', or 'g'.");
                     proceed = FALSE;
                     break;
               }

               if ( ( proceed ) && ( *(char_ptr+1) != '\0' ) )
               {
                  post_syntax_err("Extra characters after chunk size suffix.");
                  proceed = FALSE;
               }
            }
            else
            {
               multiplier = 1024 * 1024;
            }

            if ( proceed ) /* convert string to integer */
            {
               errno = 0;
               vars_ptr->chunk_len = strtoull(optarg, NULL, 10);
               assert(errno == 0);
               if ( vars_ptr->chunk_len <= 0 )
               {
                  post_syntax_err("Chunk size must be greater than 0.");
                  proceed = FALSE;
               }
            }

            if ( proceed ) 
            {
               /* apply multiplier -- should have overflow check here */
               vars_ptr->chunk_len *= multiplier;

               if ( vars_ptr->chunk_len < (5 * MAX_SEQUENCE_LEN) )
               {
                  post_syntax_err("Chunk size must be greater 50 KB.");
                  proceed = FALSE;
               }
               else
               {
                  have_chunk_size = TRUE;
               }
            }
            break;

         case 'i':
	    char_ptr = optarg;
            assert(optarg != NULL);
            while ( *char_ptr != '\0' )
            {
               if ( ! isdigit(*char_ptr) )
               {
                  post_syntax_err(
                        "chunk index must be an unsigned decimal integer.");
                  proceed = FALSE;
               }
               char_ptr++;
            }
            errno = 0;
            vars_ptr->chunk_index = strtoull(optarg, NULL, 10);
            assert(errno == 0);
            have_chunk_index = TRUE;
            break;

         case 'm':
	    char_ptr = optarg;
            assert(optarg != NULL);
            while ( *char_ptr != '\0' )
            {
               if ( ! isdigit(*char_ptr) )
               {
                  post_syntax_err(
                      "max chunk index must be an unsigned decimal integer.");
                  proceed = FALSE;
               }
               char_ptr++;
            }
            errno = 0;
            vars_ptr->max_chunk_index = strtoull(optarg, NULL, 10);
            assert(errno == 0);
            have_max_chunk_index = TRUE;
            break;

         case 'f':
            if ( strlen(optarg) >= MAX_PATH_LEN ) 
            {
               post_syntax_err("output file path too long.");
               proceed = FALSE;
            }
            else
            {
               strcpy(vars_ptr->output_file_name, optarg);
            }
            break;

         case 'p':
            if ( strlen(optarg) >= MAX_PATH_LEN ) 
            {
               post_syntax_err("output FIFO (AKA named pipe) path too long.");
               proceed = FALSE;
            }
            else
            {
               strcpy(vars_ptr->output_fifo_name, optarg);
            }
            break;

         case 'v':
            vars_ptr->verbosity += 1;
            break;

         case '?':
            proceed = FALSE;
            usage();
            break;

         default: /* this should be unreachable */
            assert(0);
            break;
      }
   }

   if ( proceed ) 
   {
      if ( ! have_chunk_size ) 
      {
         post_syntax_err("Chunk size not specified.");
         proceed = FALSE;
      } 
      else if ( ( ! have_chunk_index ) && ( ! vars_ptr->dry_run ) ) 
      {
         post_syntax_err("Chunk index not specified.");
         proceed = FALSE;
      }
      else if ( ! have_max_chunk_index )
      {
         post_syntax_err("Maximum chunk index not specified.");
         proceed = FALSE;
      }
   }

   if ( ( proceed ) && ( vars_ptr->chunk_index > vars_ptr->max_chunk_index ) )
   {
      post_syntax_err("chunk index may not exceed max chunk index.");
      proceed = FALSE;
   }

   if ( ( proceed ) && ( strlen(vars_ptr->output_file_name) > 0 ) &&
        ( strlen(vars_ptr->output_fifo_name) > 0 ) )
   {
      post_syntax_err("both output file and output fifo defined?");
      proceed = FALSE;
   }

   if ( ( proceed ) && ( optind >= argc) )
   {
      post_syntax_err("missing argument after option?");
      proceed = FALSE;
   }

   if ( ( proceed ) && ( optind < (argc - 1)) )
   {
      post_syntax_err("excess parameters?");
      proceed = FALSE;
   }

   if ( proceed )
   {
      assert(argv[optind] != NULL);

      if ( strlen(argv[optind]) >= MAX_PATH_LEN ) 
      {
         post_syntax_err("input file path too long.");
         proceed = FALSE;
      }
      else if ( strlen(argv[optind]) == 0 )
      {
         post_syntax_err("zero length input file path.");
         proceed = FALSE;
      }
      else
      {
         strcpy(vars_ptr->input_file_name, argv[optind]);
      }
   }

   return(proceed);

} /* get_params() */


/*-------------------------------------------------------------------------
 *
 * Function:    load_buf_from_input_file()
 *
 * Purpose:     Seek to the specified location in the input file and
 *              then read the indicated number of bytes into the supplied
 *              buffer.
 *
 *              Return TRUE if completely successful and FALSE otherwise.
 *
 *              If successful, return TRUE.  On failure for any reason,
 *              return FALSE.
 *
 * Return:      Success: TRUE
 *              Failure: FALSE
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static int
load_buf_from_file(off_t offset,
                   off_t buf_len,
                   char *buf_ptr,
                   struct vc_vars * vars_ptr)
{
   char *local_buf_ptr;
   int read_itteration = 0;
   int saved_errno;
   int proceed = TRUE;
   size_t local_buf_len;
   ssize_t bytes_read = 0;
   ssize_t result;
   off_t seek_result;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert(vars_ptr->input_fd >= 0);
   assert(vars_ptr->input_file_len > 0);
   assert(offset < vars_ptr->input_file_len);
   assert(buf_len > 0);
   assert(vars_ptr->input_file_len >= offset + buf_len);

   if ( proceed )
   {
      errno = 0;
      if ( (off_t)-1 == lseek(vars_ptr->input_fd, offset, SEEK_SET) )
      {
         saved_errno = errno;

         proceed = FALSE;

         fprintf(stderr, "\nCant't seek to %llu in input file.\n", offset);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      }
   }

   local_buf_ptr = buf_ptr;
   local_buf_len = (size_t)buf_len;

   while ( ( proceed ) && ( bytes_read < (ssize_t)buf_len ) )
   {
      errno = 0;
      result = read(vars_ptr->input_fd, (void *)local_buf_ptr, local_buf_len);
      saved_errno = errno;

      if ( result == -1 )
      {
         proceed = FALSE;

         fprintf(stderr, 
                "\nCant't read from in input file on itteration %d.\n", 
                read_itteration);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      }
      else 
      {
         assert(result <= local_buf_len);

         read_itteration++;

         bytes_read += result;
         local_buf_ptr += result;
         local_buf_len -= result;
      }
   }

   if ( proceed ) 
      assert(bytes_read == (ssize_t)buf_len);

   return(proceed);

} /* load_buf_from_file() */


/*-------------------------------------------------------------------------
 * Function:    map_test_buf()
 *
 * Purpose:     Scan the supplied test buffer to find the offset into
 *		the buffer of the first character of every FASTQ sequence
 *		that has its first character in the buffer.
 *
 *		Unfortunately, the structure of a FASTQ file was not 
 *		designed to make this easy, and the documentation on 
 *		on the FASTQ file format is not very good.  
 *
 *		I have extracted the following BNF for FASTQ:
 *
 *		<FASTQ file> ::= <sequence>*
 *
 *		<sequence> ::= <line 1> <line 2> <line 3> <line 4>
 *
 *		<line 1> ::= '@' <seq_id> [<seq_desc>]  <EOL>
 *
 *		<line 2> ::= <raw_sequence_char>+ <EOL>
 *
 *		<line 3> ::= '+' [<seq_desc>] <EOL>
 *
 *		<line 4> ::= <quality values> <EOL>
 *
 *		<raw_sequence_char> ::= 
 *                ('A' | 'a' | 'C' | 'c' | 'G' | 'g' | 'N' | 'n' | 'T' | 't')+
 *
 *		<quality_values> ::= <quality char>+
 *
 *		<quality char> ::= and ASCII character in the set:
 *                  !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLM
 *                  NOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
 *
 *		Key:	<>  indicates non-terminal
 *			[]  indicates option
 *		        ()  indicates grouping
 *			*   closure
 *			+   positive closure 
 *
 *		Unfortunately, I haven't been able to find any restrictions
 *		on the <seq_id> or the <seq_desc> -- other than to observe 
 *		that white space seems to be allowed, as do all printing
 *		characters 
 *
 *		Note the following semantic requirements:
 *
 *		In any given sequence, the strings replacing 
 *		the <raw_sequence_data> and <quality_values> non
 *		terminals must have the same length, and may contain 
 *		no white space.
 *
 *		If <seq_id> appears in <line 3>, then the strings 
 *		replacing the <seq_id> in <line 1> and <line 3> must
 *		be identical.
 *
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes      none.
 *
 *-------------------------------------------------------------------------
 */
static int
map_test_buf(off_t test_buf_offset,
             char test_buf[], 
             int test_buf_len,
             struct seq_header map[],
             int map_len,
             int * seq_count_ptr,
             struct vc_vars * vars_ptr)
{
   int best_offset;
   int min_errs;
   int i;
   int j;
   int bad_line_count = 0;
   int line_count = 0;
   int assignment_errors[4];
   int proceed = TRUE;
   struct line_header line_map[2 * MAX_SEQUENCE_LEN];

   proceed = map_test_buf_lines(test_buf_offset, 
                                test_buf, 
                                test_buf_len,
                                line_map,
                                2 * MAX_SEQUENCE_LEN,
                                &line_count,
                                vars_ptr);

   /* At this point, we have broken the test buffer into lines, scanned
    * each line, and marked it with flags indicating the possible type 
    * of FASTQ line it could be.  
    *
    * Scan the line map, and look for lines that can't appear in a FASTQ
    * file.  If we find any, post an error message and quit.
    */
   if ( proceed )
   {
      for ( i = 0; i < line_count; i++ )
      {
         if ( line_map[i].id_flags == 0 )
         {
            bad_line_count++;

            if ( vars_ptr->verbosity >= 3)
            {
                fprintf(stderr, "\nERROR: invalid line at offset %llu.\n",
                        (unsigned long long)(line_map[i].offset));
                fprintf(stderr, "	line body: \"%s\"\n",
                        line_map[i].line_ptr);
            }
         }
      }

      if ( bad_line_count > 0 )
      {
         proceed = FALSE;

         fprintf(stderr, 
                "\nERROR: Encountered invalid line(s) in input file.\n");
      }
      else
      {
         /* If this assertion fails, MAX_SEQUENCE_LEN is probably too small */
         assert(line_count >= 16);
      }
   }

   /* If we get this far without errors, we know that all the lines in 
    * the line map are lines that could conceiveably appear in a FASTQ
    * file.  Must now assemble a sequence map from the line map.
    *
    * Since the lines must appear in a well defined order, there are 
    * only four possible assignments.  We will examine them all, and 
    * compute how many errors each one produces.  
    */

   if ( proceed )
   {
      best_offset = -1;
      min_errs = line_count + 1;

      for ( i = 0; i < 4; i++ )
      {
         count_assignment_errors(i, line_map, line_count, 
                                 &(assignment_errors[i]));

         if ( assignment_errors[i] < min_errs )
         {
            min_errs = assignment_errors[i];
            best_offset = i;
         }
      }

      assert((best_offset >= 0) && ( best_offset <= 3));

      /* At some point, we may wish to allow error tollerance here. But 
       * for now, scream and die if there is no assignememtn without 
       * errors.
       */

      if ( min_errs > 0 ) 
      {
         proceed = FALSE;

         fprintf(stderr, 
                 "\n\nERROR: Unable to construct error free sequence map.\n");
         fprintf(stderr, "	Best offset = %d, Best err count = %d.\n\n",
                 best_offset, min_errs);

         if ( vars_ptr->verbosity >= 3 ) 
         {
            dump_line_map(line_map, line_count);
         }
      }
   }

   if ( proceed ) /* verify that the best offset is unique */
   {
      /* if it turns out that the best offset isn't unique, we have 
       * ways of disambiguating the situation.  However, for now
       * just test to see if the situation occurs.  
       */
      for ( i = 0; i < 4 ; i++ )
      {
         assert((i == best_offset) ||
                (assignment_errors[best_offset] < assignment_errors[i]));
      }
   }

   if ( proceed ) /* construct the sequence map */
   {
      j = 0;
      for (i = best_offset; i < line_count; i += 4)
      {
         if ( i + 3 < line_count )
         {
            assert(j < map_len);

            map[j].offsets[0]  = line_map[i].offset;
            map[j].lines[0]    = line_map[i].line_ptr;
            map[j].line_len[0] = line_map[i].line_len;
            map[j].flags[0]    = line_map[i].id_flags;

            map[j].offsets[1]  = line_map[i + 1].offset;
            map[j].lines[1]    = line_map[i + 1].line_ptr;
            map[j].line_len[1] = line_map[i + 1].line_len;
            map[j].flags[1]    = line_map[i + 1].id_flags;

            map[j].offsets[2]  = line_map[i + 2].offset;
            map[j].lines[2]    = line_map[i + 2].line_ptr;
            map[j].line_len[2] = line_map[i + 2].line_len;
            map[j].flags[2]    = line_map[i + 2].id_flags;

            map[j].offsets[3]  = line_map[i + 3].offset;
            map[j].lines[3]    = line_map[i + 3].line_ptr;
            map[j].line_len[3] = line_map[i + 3].line_len;
            map[j].flags[2]    = line_map[i + 3].id_flags;

            j++;
         }
      }

      *seq_count_ptr = j + 1;
   }

   return(proceed);

} /* map_test_buf() */


/*-------------------------------------------------------------------------
 * Function:    map_test_buf_lines()
 *
 * Purpose:     Scan the test buffer, and replace all new lines with 
 *		null characters.  Construct an index of all lines that 
 *		start in the test buffer in the supplied array of 
 *		struct line_header.  In passing, examine characters that
 *		appear in the lines, and set flags indicating possible
 *		identifications.
 *
 * Return:       TRUE if successful, and FALSE otherwise.
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static int
map_test_buf_lines(off_t test_buf_offset,
                   char test_buf[], 
                   int test_buf_len,
                   struct line_header line_map[],
                   int map_len,
                   int *line_count_ptr,
                   struct vc_vars * vars_ptr)
{
   enum scan_state 
   { 
      finding_first_line, 
      starting_line, 
      scanning_line, 
      nulling_new_line
   } state = finding_first_line;
   int could_be_line_1;
   int could_be_line_2;
   int could_be_line_3;
   int could_be_line_4;
   int i = 0;
   int j = 0;
   int line_count = 0;
   int proceed = TRUE;

   assert(test_buf != NULL);
   assert(test_buf_len > 0);
   assert(line_map != NULL);
   assert(map_len > 0);
   assert(line_count_ptr != NULL);
   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);

   if ( ( test_buf[i] == '\r' ) || ( test_buf[i] == '\n' ) )
   {
      state = nulling_new_line;
   }

   while ( ( proceed ) && ( i < test_buf_len ) )
   {
      switch ( state )
      {
         case finding_first_line:
            if ( ( i + 1 < test_buf_len ) &&
                 ( ( test_buf[i + 1] == '\r' ) ||
                   ( test_buf[i + 1] == '\n' ) ) )
            {
               state = nulling_new_line;
            }
            i++;
            break;

         case starting_line:
            assert(j < map_len);
            
            /* initialze the map entry for this line.  Setting line_len
             * to zero below is not a typo, as we are going to change state
             * to scanning_line and fall through to that state.
             */
            line_map[j].line_ptr = &(test_buf[i]);
            line_map[j].offset = test_buf_offset + (off_t)i;
            line_map[j].line_len = 0;
            line_map[j].id_flags = 0;

            if ( test_buf[i] == '@' )
            {
               could_be_line_1 = TRUE;
            }
            else
            {
               could_be_line_1 = FALSE;
            }

            could_be_line_2 = TRUE;

            if ( test_buf[i] == '+' )
            {
               could_be_line_3 = TRUE;
            }
            else
            {
               could_be_line_3 = FALSE;
            }

            could_be_line_4 = TRUE;

            /* set state to case scanning_line and fall through 
             * to that state 
             */
            state = scanning_line;

            /**** NOTE FALL THROUGH *****/

         case scanning_line:
            if ( ( test_buf[i] == '\r' ) || ( test_buf[i] == '\n' ) )
            {
               /* We have found the end of the line -- finish setting 
                * up the line map entry, increment j, and line_count,
                * and set state to nulling_new_line.
                *
                * Note that we must not increment i, as i must be the
                * index of the initial new line character on entry to 
                * the nulling_new_line state.
                */
               assert(j == line_count);
               assert(line_map[j].line_len > 0);

               if ( could_be_line_1 ) line_map[j].id_flags |= AT_LINE;
               if ( could_be_line_2 ) line_map[j].id_flags |= SEQ_LINE;
               if ( could_be_line_3 ) line_map[j].id_flags |= PLUS_LINE;
               if ( could_be_line_4 ) line_map[j].id_flags |= QUAL_LINE;

               j++;
               line_count++;

               assert(j < map_len);

               state = nulling_new_line;
            }
            else if ( (i + 1) >= test_buf_len )
            {
               /* unfinished line -- NULL out the current line map entry
                * and increment i.
                */
               line_map[j].line_ptr = NULL;
               line_map[j].offset = (off_t)0;
               line_map[j].line_len = 0;
               line_map[j].id_flags = 0;

               /* increment i to get us out of the while loop */
               i++;
            }
            else 
            {
               /* update possible line id flags */
               could_be_line_2 = could_be_line_2 && IS_SEQ_CHAR(test_buf[i]);
               could_be_line_4 = could_be_line_4 && IS_QUAL_CHAR(test_buf[i]);

               /* increment the line length */
               line_map[j].line_len += 1;

               /* increment i */
               i++;
            }
            break;

         case nulling_new_line:
            /* No telling what machine this file was prepared on. Thus
             * accept either unix / MacOS X style new lines "\n", DOS
             * style new lines ("\r\n"), or old Macintosh style new
             * lines ("\r").
             *
             * Note that in a FASTQ file, we can have no empty lines,
             * so flag an error if we don't see one of the above 
             * sequences.
             */
            assert((test_buf[i] == '\r') || (test_buf[i] == '\n'));

            if ( test_buf[i] == '\n' )
            {
               /* Unix new line -- null this character and increment i.
                * Will verify that the next character is printable 
                * of off the buffer later.
                */
               test_buf[i] = '\0';
               i++;
            } 
            else if ( test_buf[i] == '\r' )
            {
               /* should be either an old Macintosh new line or an 
                * MSDOS newline.
                */
               if ( i + 1 >= test_buf_len )
               {
                  /* it doesn't matter as we are at the end of the 
                   * buffer -- just null the character and increment i
                   */
                  test_buf[i] = '\0';
                  i++;
               } 
               else if ( test_buf[i] == '\n' )
               {
                  /* It is an MSDOS new line -- null both characters
                   * and increment i by 2.  Will verify that the next 
                   * character is printable or off the buffer later.
                   */
                 test_buf[i] = '\0';
                 test_buf[i+1] = '\0';
                 i += 2;
               }
               else 
               {
                  /* it is either an old Macintosh new line or an error.
                   * null the character and check for error later. 
                   */
                  test_buf[i] = '\0';
                  i++;
               }
            }
               
            /* At this point, we have nulled the new line character(s),
             * and must verify that either we have reached the end of 
             * the buffer, or that the next character is printable.
             * In either case, we can set state to starting_line.
             *
             * Otherwise, flag a file syntax error indicating ill formed
             * end of line marker.
             */
            if ( ( i >= test_buf_len ) || ( isgraph(test_buf[i]) ) )
            {
               state = starting_line;
            }
            else
            {
               proceed = FALSE;

               fprintf(stderr, 
                  "\nIll formed EOL at offset %llu in input file \"%s\".\n",
                  (test_buf_offset + i - 1), vars_ptr->input_file_name);
            }
            break;

         default: /* should be unreachable */
            assert(0);
            break;

      } /* end switch */
   } /* end while */

   /* If this assertion failes, the MAX_SEQUENCE_LEN is too small */
   assert(line_count >= 8);

   /* null out the rest of the line map */
   for ( j = line_count; j < map_len; j++ )
   {
      line_map[j].line_ptr = NULL;
      line_map[j].offset = (off_t)0;
      line_map[j].line_len = 0;
      line_map[j].id_flags = 0;
   }

   *line_count_ptr = line_count;

   return(proceed);

} /* map_test_buf_lines() */


/*-------------------------------------------------------------------------
 * Function:     open_input_file()
 *
 * Purpose:      Attempt to open the input file.  If successful, return
 *               TRUE.  Otherwise, issue an error message, and return
 *               FALSE.
 *
 *		 Note that we use file descriptor based I/O here because
 *		 we are working with large files, and thus may have to 
 *		 seek to locations that cannot be represented with a long.
 *
 * Return:       TRUE if successful, and FALSE otherwise.
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static int
open_input_file(struct vc_vars * vars_ptr)
{
   int proceed = TRUE;
   int fd;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert(strlen(vars_ptr->input_file_name) > 0);
   assert(vars_ptr->input_fd == -1);

   errno = 0;
   fd = open(vars_ptr->input_file_name, O_RDONLY, S_IRUSR);
   saved_errno = errno;

   if ( fd == -1 )
   {
      proceed = FALSE;

      fprintf(stderr, "\nCant't open input file \"%s\".\n",
              vars_ptr->input_file_name);
      fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
              strerror(saved_errno));
   }
   else
   {
      vars_ptr->input_fd = fd;
   }

   return(proceed);

} /* open_input_file() */


/*-------------------------------------------------------------------------
 * Function:     post_err_mssg()
 *
 * Purpose:      Post the supplied error message, and return.
 *
 * Return:       void
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static void
post_err_mssg(char * err_msg)
{
   assert(err_msg != NULL);

   fprintf(stderr, "\nERROR: %s\n\n", err_msg);

   return;

} /* post_err_mssg() */


/*-------------------------------------------------------------------------
 * Function:     post_syntax_err()
 *
 * Purpose:      Post the supplied syntax error message, display the
 *		 usage, and return.
 *
 * Return:       void
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static void
post_syntax_err(char * err_msg)
{
   assert(err_msg != NULL);
   fprintf(stderr, "\nSYNTAX ERROR: %s\n\n", err_msg);
   usage();

   return;

} /* post_syntax_err() */


/*-------------------------------------------------------------------------
 * Function:    setup_output_fd
 *
 * Purpose:	Setup the file descriptor to be used for output.
 *
 *		If neither vars_ptr->output_file_name nor 
 *		vars_ptr->output_fifo_name are defined, simply
 *		set vars_ptr->output_fd to STDOUT_FILENO.
 *		
 *		If vars_ptr->output_file_name is defined, assert
 *		that vars_ptr->output_fifo_name is undefined.  Then
 *		attempt to open the file O_WRONLY / O_APPEND O_CREAT.
 *		If successful, set vars_ptr->output_fd to the file
 *		descriptor of the newly opened file, and return TRUE.
 *		On failure, generate an error message and return FALSE.
 *
 *		If vars_ptr->output_fifo_name is defined, assert that 
 *		vars_ptr->output_file_name is undefined.  Then attempt
 *		to create the fifo -- as it may already exist, do not
 *		flag an error mkfifo() returns EEXIST.  On all other 
 *		errors, generate an error messages and return FALSE.
 *		Assuming no fatal errors in mkfifo(), open the fifo
 *		O_WRONLY.  On failure, generate an error message and
 *		return FALSE.  Otherwise, set vars_ptr->output_fd to 
 *		the file descriptor associated with the newly opened
 *		FIFO, and return TRUE.
 *
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */
static int
setup_output_fd(struct vc_vars * vars_ptr)
{
   int proceed = TRUE;
   int fd;
   int result;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert( ! ( ( strlen(vars_ptr->output_file_name) > 0 ) && \
               ( strlen(vars_ptr->output_fifo_name) > 0 ) ) );
   assert(vars_ptr->output_fd == -1);

   if ( ( strlen (vars_ptr->output_file_name) == 0 ) &&
        ( strlen (vars_ptr->output_fifo_name) == 0 ) )
   {
      /* send output to stdout */
      vars_ptr->output_fd = STDOUT_FILENO;
   }
   else if ( strlen (vars_ptr->output_file_name) > 0 )
   {
      /* send output to a normal file */

      errno = 0;
      fd = open(vars_ptr->output_file_name, 
                O_WRONLY|O_APPEND|O_CREAT, S_IRUSR|S_IWUSR);
      saved_errno = errno;

      if ( fd == -1 )
      {
         proceed = FALSE;

         fprintf(stderr, "\nCant't open output file \"%s\".\n",
                 vars_ptr->output_file_name);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      }
      else
      {
         vars_ptr->output_fd = fd;
      }
   }
   else
   {
      /* send output to a FIFO (AKA named pipe) */

      assert(strlen(vars_ptr->output_fifo_name) > 0);

      /* first, try to create the FIFO.  As it may already exist,
       * ignore a failure with errno == EEXIST.
       */
      errno = 0;
      result = mkfifo(vars_ptr->output_fifo_name, S_IRUSR|S_IWUSR);
      saved_errno = errno;

      if ( ( result == -1 ) && ( saved_errno != EEXIST ) )
      {
         proceed = FALSE;

         fprintf(stderr, "\nCant't create output fifo \"%s\".\n",
                 vars_ptr->output_fifo_name);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      }

      if ( proceed )
      {
         errno = 0;
         fd = open(vars_ptr->output_fifo_name, O_WRONLY, O_APPEND);
         saved_errno = errno;

         if ( fd == -1 )
         {
            proceed = FALSE;

            fprintf(stderr, "\nCant't open output fifo \"%s\".\n",
                    vars_ptr->output_fifo_name);
            fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                    strerror(saved_errno));
         }
         else
         {
            vars_ptr->output_fd = fd;
         }
      }
   }

   return(proceed);

} /* setup_output_fd() */


/*-------------------------------------------------------------------------
 * Function:    takedown_output_fd()
 *
 * Purpose:     Takedown the output file descriptor.  What has to be 
 *		done depends on its value.
 *
 *		if vars_ptr->output_fd == STDOUT_FILENO just set 
 *		vars_ptr->output_fd = -1, as we don't want to close
 *		stdout.
 *
 *		Any other value of vars_ptr->output_fd indicates that 
 *		it points to either a regular file of a FIFO (AKA 
 *		named pipe).  In either case we must use the close()
 *		system call to close it.  Generate an error 
 *		message on failure.
 *
 *		Note that we use file descriptor based I/O here because
 *		we are working with large files, and thus may have to 
 *		seek to locations that cannot be represented with a long.
 * 
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static void
takedown_output_fd(struct vc_vars * vars_ptr)
{
   int result;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert( ! ( ( strlen(vars_ptr->output_file_name) > 0 ) && \
               ( strlen(vars_ptr->output_fifo_name) > 0 ) ) );
   assert(vars_ptr->output_fd != -1);

   if ( vars_ptr->output_fd == STDOUT_FILENO )
   {
      vars_ptr->output_fd = -1;
   }
   else if ( strlen(vars_ptr->output_file_name) > 0 )
   {
      errno = 0;
      result = close(vars_ptr->output_fd);
      saved_errno = errno;

      if ( result == -1 )
      {
         fprintf(stderr, "\nCant't close output file \"%s\".\n",
                 vars_ptr->output_file_name);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      }
   }
   else
   {
      assert( strlen(vars_ptr->output_fifo_name) > 0 );

      errno = 0;
      result = close(vars_ptr->output_fd);
      saved_errno = errno;

      if ( result == -1 )
      {
         fprintf(stderr, "\nCant't close output fifo \"%s\".\n",
                 vars_ptr->output_fifo_name);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      }
   }

   return;

} /* takedown_output_fd() */


/*-------------------------------------------------------------------------
 * Function:     usage()
 *
 * Purpose:      Print usage message.
 *
 * Return:       void
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static void
usage (void)
{
 fprintf(stdout, "\
Virtual Chunker (vc)\n\
\n\
purpose: Extract the specified chunk from the supplied FASTQ file, and \n\
   write it to stdout, a file, or a FIFO (AKA named pipe) as specified.\n\
\n\
usage: vc [OPTIONS] <input FASTQ file>\n\
\n\
   OPTIONS\n\
      -d                 Dry run.  Compute all chunk boundaries. \n\
                         Display if verbosity > 0. No output to file. \n\
\n\
      -s n[k|m|g]        Chunk size in kilo-bytes (if \"k\" suffix appears),\n\
                         mega-bytes (if no suffix, or \"m\" suffix), or \n\
                         giga-bytes (if \"g\" suffix appears). \n\
\n\
      -i n               Index of the chunk to extract. \n\
                         Ignored if dry run is specified. \n\
\n\
      -m n               Max chunk index. \n\
\n\
      -f <file name>     Name of the file to which to write the chunk.\n\
\n\
      -p <FIFO name>     Name of the FIFO (AKA named pipe) to which to  \n\
                         write the chunk.\n\
\n\
      -v                 Verbosity.  Repeat for more output. \n\
                           verbosity == 0 ==> output on error only \n\
                           verbosity == 1 ==> terse functional output \n\
                           verbosity == 2 ==> verbose functional output \n\
                           verbosity >= 3 ==> debugging output \n\
\n\
   If chunk index is equal to max chunk index, then the delivered chunk \n\
   will extend to the end of file, regardless of size. \n\
\n\
   The -f and -p options are mutually exclusive.  If neither is provided, \n\
   output is to stdout. \n\
\n\
   If specified, both the output file and the output FIFO will be created \n\
   vc if they do not exist on entry.  Unlinking the FIFO is the callers \n\
   responsibility. \n\n");

   return;

} /* usage() */


/*-------------------------------------------------------------------------
 *
 * Function:    write_buf_to_output_fd()
 *
 * Purpose:     Write the contents of the supplied buffer to the file 
 *		descriptor listed in vars_ptr->output_fd.  This may be 
 *		either sdtoud, a regular file, or a FIFO (AKA named pipe).
 *
 *              Return TRUE if completely successful and FALSE otherwise.
 *
 * Return:      Success: TRUE
 *              Failure: FALSE
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static int
write_buf_to_output_fd(off_t buf_len,
                       char *buf_ptr,
                       struct vc_vars * vars_ptr)
{
   char *local_buf_ptr;
   int write_itteration = 0;
   int saved_errno;
   int proceed = TRUE;
   size_t local_buf_len;
   ssize_t bytes_written = 0;
   ssize_t result;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert(vars_ptr->output_fd != -1);
   assert(buf_len > 0);
   assert(buf_ptr != NULL);

   local_buf_ptr = buf_ptr;
   local_buf_len = (size_t)buf_len;

   while ( ( proceed ) && ( bytes_written < (ssize_t)buf_len ) )
   {
      errno = 0;
      result = write(vars_ptr->output_fd, (void *)local_buf_ptr, 
                     local_buf_len);
      saved_errno = errno;

      if ( result == -1 )
      {
         proceed = FALSE;

         fprintf(stderr, 
                "\nCant't write to output fd on itteration %d.\n", 
                write_itteration);
         fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
                 strerror(saved_errno));
      }
      else 
      {
         assert(result <= local_buf_len);

         write_itteration++;

         bytes_written += result;
         local_buf_ptr += result;
         local_buf_len -= result;
      }
   }

   if ( proceed ) 
      assert(bytes_written == (ssize_t)buf_len);

   return(proceed);

} /* write_buf_to_output_fd() */


/*-------------------------------------------------------------------------
 * Function:    write_chunk
 *
 * Purpose:	Load the desired chunk from the input file, and write 
 *		it to the output specified by vars_ptr->output_fd.
 *
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */
static int
write_chunk(struct vc_vars * vars_ptr)
{
   char * buf_ptr = NULL;
   int proceed = TRUE;
   int saved_errno;
   off_t offset;
   off_t offset_this_read;
   off_t bytes_remaining;
   off_t bytes_this_read;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == VC_VARS_MAGIC);
   assert(vars_ptr->input_fd != -1);
   assert(vars_ptr->output_fd != -1);

   buf_ptr = malloc(IO_BUF_SIZE);

   if ( buf_ptr == NULL )
   {
      post_err_mssg("Can't allocate I/O buffer.");
      proceed = FALSE;
   }

   if ( proceed )
   {
      offset = vars_ptr->actual_chunk_start;
      bytes_remaining = vars_ptr->actual_chunk_end - 
                        vars_ptr->actual_chunk_start + 1;
   }

   while ( ( proceed ) && ( bytes_remaining > 0 ) )
   {
      if ( bytes_remaining >= IO_BUF_SIZE )
      {
         offset_this_read = offset;
         offset += IO_BUF_SIZE;
         bytes_this_read = IO_BUF_SIZE;
         bytes_remaining -= IO_BUF_SIZE;
      }
      else
      {
         bytes_this_read = bytes_remaining;
         bytes_remaining -= bytes_this_read;
         offset_this_read = offset;
         offset += bytes_this_read;
      }

      proceed = load_buf_from_file(offset_this_read, bytes_this_read, 
                                   buf_ptr, vars_ptr);

      if ( proceed ) 
         proceed = write_buf_to_output_fd(bytes_this_read, buf_ptr, vars_ptr);
   }

   assert(offset == vars_ptr->actual_chunk_end + 1);
   assert(bytes_remaining == 0 );

   if ( buf_ptr != NULL )
   {
      free(buf_ptr);
      buf_ptr = NULL;
   } 

   return(proceed);

} /* write_chunk() */


/*-------------------------------------------------------------------------
 * Function:     main()
 *
 * Purpose:      main routine for vc (virtual chunker).
 *
 * Changes:      none.
 *
 *-------------------------------------------------------------------------
 */

int
main(int argc,
     char *argv[])
{
   int proceed = TRUE;
   struct vc_vars vars = 
   {
      /* magic               = */ VC_VARS_MAGIC,
      /*                       */ 
      /* dry_run             = */ FALSE,
      /* chunk_len           = */ (unsigned long long)0,
      /* chunk_index         = */ (unsigned long long)0,
      /* max_chunk_index     = */ (unsigned long long)0,
      /* input_file_name     = */ "",
      /* output_file_name    = */ "",
      /* output_fifo_name    = */ "", 
      /* verbosity           = */ 0,
      /*                       */
      /* input_file_len      = */ (off_t)0,
      /* nominal_chunk_start = */ (off_t)0,
      /* nominal_chunk_end   = */ (off_t)0,
      /* actual_chunk_start  = */ (off_t)0,
      /* actual_chunk_end    = */ (off_t)0,
      /* input_fd            = */ -1,
      /* output_fd           = */ -1,
      /*                       */
      /* first_lines         = */ {NULL, NULL, NULL, NULL},
      /* last_lines          = */ {NULL, NULL, NULL, NULL},
      /* first_lines_offsets = */ {(unsigned long long)0, 
      /*                       */  (unsigned long long)0,
      /*                       */  (unsigned long long)0,
      /*                       */  (unsigned long long)0},
      /* last_lines_offsets  = */ {(unsigned long long)0, 
      /*                       */  (unsigned long long)0,
      /*                       */  (unsigned long long)0,
      /*                       */  (unsigned long long)0}
   };

   /* We will be working with large files (10's of GB), and thus 
    * we will be using the OS file system calls for I/O (open(),
    * read(), lseek(), etc) instead of the file I/O calls in the 
    * C library since they use off_t for offsets instead of longs.
    *
    * However, in some C library implementations, off_t is a 4 byte
    * value.  To catch this, assert that sizeof(off_t) >= 8.
    */
   assert(sizeof(off_t) >= 8);

   proceed = get_params(argc, argv, &vars);

   if ( vars.verbosity >= 3 ) dump_params(&vars);

   if ( proceed ) proceed = get_input_file_len(&vars);

   if ( ( proceed ) && 
        ( vars.chunk_len > (unsigned long long)(vars.input_file_len) ) )
   {
      post_err_mssg("chunk size greater than input file size.");
      proceed = FALSE;
   }

   if ( ( proceed ) &&
        ( ((vars.chunk_len * vars.max_chunk_index) + (5 * MAX_SEQUENCE_LEN))
          >= vars.input_file_len ) )
   {
      post_err_mssg("Last chunk too small.  Decrease max chunk index?");
      proceed = FALSE;
   }

   if ( proceed ) proceed = open_input_file(&vars);

   if ( vars.dry_run ) /* compute all chunk end points, then exit */
   {
      if ( proceed ) proceed = dry_run(&vars);
   }
   else /* normal execution */
   {
      if ( proceed ) proceed = compute_nominal_chunk_endpoints(&vars);

      if ( proceed ) proceed = compute_actual_chunk_start_point(&vars);

      if ( proceed ) proceed = compute_actual_chunk_end_point(&vars);

      if ( proceed ) proceed = setup_output_fd(&vars);

      if ( proceed ) proceed = write_chunk(&vars);

      if ( vars.output_fd != -1 ) takedown_output_fd(&vars);
   }

   if ( vars.input_fd != -1 ) close_input_file(&vars);

   return(!proceed);

} /* main */
