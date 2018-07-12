#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <assert.h>


/* #defines */

#define FALSE 0
#define TRUE  1
#define MAX_PATH_LEN		1024
#define MAX_LINE_LEN		10240
#define IO_BUF_SIZE		((size_t)(2 * 1024))
#define FIFO_2_FILE_VARS_MAGIC	0x2345

struct fifo2file_vars
{
   unsigned magic;  /* always set to FIFO_2_FILE_VARS_MAGIC */

   /* parameters */
   char input_fifo_name[MAX_PATH_LEN + 1];
   char output_file_name[MAX_PATH_LEN + 1];
   int verbosity;

   /* variables */
   int input_fd;
   int output_fd;
};

/* function declarations */
static void close_input_fifo(struct fifo2file_vars * vars_ptr);
static int copy_data_from_fifo_to_file(struct fifo2file_vars * vars_ptr);
static void dump_params(struct fifo2file_vars * vars_ptr);
static int get_params(int argc, char *argv[], 
                      struct fifo2file_vars * vars_ptr);
static int load_buf_from_fifo(off_t buf_len, char *buf_ptr, 
                              off_t *bytes_read_ptr, int *eof_ptr,
                              struct fifo2file_vars * vars_ptr);
static int open_input_fifo(struct fifo2file_vars * vars_ptr);
static void post_err_mssg(char * err_msg);
static void post_syntax_err(char * err_msg);
static int setup_output_fd(struct fifo2file_vars * vars_ptr);
static void takedown_output_fd(struct fifo2file_vars * vars_ptr);
static void usage(void);
static int write_buf_to_output_fd(off_t buf_len, char *buf_ptr, 
                                  struct fifo2file_vars * vars_ptr);



/*-------------------------------------------------------------------------
 * Function:    close_input_fifo()
 *
 * Purpose:     Attempt to close the input fifo.  Generate an error 
 *		error message on failure.
 * 
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static void
close_input_fifo(struct fifo2file_vars * vars_ptr)
{
   int result;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);
   assert(strlen(vars_ptr->input_fifo_name) > 0);
   assert(vars_ptr->input_fd != -1);

   errno = 0;
   result = close(vars_ptr->input_fd);
   saved_errno = errno;

   if ( result == -1 )
   {
      fprintf(stderr, "\nCant't close input fifo \"%s\".\n",
              vars_ptr->input_fifo_name);
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
 * Function:    copy_data_from_fifo_to_file()
 *
 * Purpose:     While not EOF or error, read data from the input fifo,
 *		and write it to the output file.  
 *
 * Return:      TRUE if no errors detected, and FALSE otherwise.
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */
static int
copy_data_from_fifo_to_file(struct fifo2file_vars * vars_ptr)
{
   char buf[IO_BUF_SIZE + 1];
   int eof = FALSE;
   int proceed = TRUE;
   off_t bytes_read;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);

   while ( ( proceed ) && ( ! eof ) )
   {
      proceed = load_buf_from_fifo((off_t)IO_BUF_SIZE, &(buf[0]), 
                                   &bytes_read, &eof, vars_ptr);

      if ( ( proceed ) && ( ! eof ) )
         proceed = write_buf_to_output_fd(bytes_read,  &(buf[0]), vars_ptr);
   }

   return(proceed);

} /* copy_data_from_fifo_to_file */


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
dump_params(struct fifo2file_vars * vars_ptr)
{
   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);

   fprintf(stdout, "\n\nParameters:\n\n");
   fprintf(stdout, "   input_fifo_name  = \"%s\"\n", 
           vars_ptr->input_fifo_name);
   fprintf(stdout, "   output_file_name = \"%s\"\n", 
           vars_ptr->output_file_name);
   fprintf(stdout, "   verbosity        = %d\n\n", vars_ptr->verbosity);

   return;

} /* dump_params() */


/*-------------------------------------------------------------------------
 * Function:     open_input_fifo()
 *
 * Purpose:      Attempt to open the input fifo.  If successful, return
 *               TRUE.  Otherwise, issue an error message, and return
 *               FALSE.
 *
 * Return:       TRUE if successful, and FALSE otherwise.
 *
 * Changes:	 none.
 *
 *-------------------------------------------------------------------------
 */
static int
open_input_fifo(struct fifo2file_vars * vars_ptr)
{
   int proceed = TRUE;
   int fd;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);
   assert(strlen(vars_ptr->input_fifo_name) > 0);
   assert(vars_ptr->input_fd == -1);

   errno = 0;
   fd = mkfifo(vars_ptr->input_fifo_name, S_IRUSR|S_IWUSR);
   saved_errno = errno;

   if ( fd == -1 )
   {
      proceed = FALSE;

      fprintf(stderr, "\nCant't create input fifo \"%s\".\n",
              vars_ptr->input_fifo_name);
      fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
              strerror(saved_errno));
   }

   errno = 0;
   fd = open(vars_ptr->input_fifo_name, O_RDONLY);
   saved_errno = errno;

   if ( fd == -1 )
   {
      proceed = FALSE;

      fprintf(stderr, "\nCant't open input fifo \"%s\".\n",
              vars_ptr->input_fifo_name);
      fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
              strerror(saved_errno));
   }
   else
   {
      vars_ptr->input_fd = fd;
   }

   return(proceed);

} /* open_input_fifo() */


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
           struct fifo2file_vars * vars_ptr)
{
   char *char_ptr;
   int opt;
   int proceed = TRUE; 

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);

   if ( argc < 1 ) 
   {
      post_syntax_err("too few arguments.");
      proceed = FALSE;
   }

   while ( ( proceed ) &&
           ( (opt = getopt(argc, argv, "f:v?")) != -1 ) )
   {
      switch ( opt ) 
      {
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
         post_syntax_err("input fifo path too long.");
         proceed = FALSE;
      }
      else if ( strlen(argv[optind]) == 0 )
      {
         post_syntax_err("zero length input fifo path.");
         proceed = FALSE;
      }
      else
      {
         strcpy(vars_ptr->input_fifo_name, argv[optind]);
      }
   }

   return(proceed);

} /* get_params() */


/*-------------------------------------------------------------------------
 *
 * Function:    load_buf_from_fifo()
 *
 * Purpose:     Attempt to read a buffer from the input fifo.  
 *
 *		If either a full or partial buffer is read, simply return
 *		with *bytes_read_ptr set appropriately.
 *
 *		Since we are using blocking I/O, read should only
 *		return 0 when the fifo is closed at the other end.
 *		Thus set *eof_ptr to TRUE, *bytes_read_ptr to 0, 
 *		and return if this happens.
 *
 *		If an error is reported, generate an error message,
 *		and return FALSE.
 *
 * Return:      TRUE if no errors are detected, and FALSE otherwise.
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static int
load_buf_from_fifo(off_t buf_len,
                   char *buf_ptr,
                   off_t *bytes_read_ptr,
                   int *eof_ptr,
                   struct fifo2file_vars * vars_ptr)
{
   char *local_buf_ptr;
   int saved_errno;
   int proceed = TRUE;
   size_t local_buf_len;
   ssize_t bytes_read = 0;
   off_t seek_result;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);
   assert(vars_ptr->input_fd >= 0);
   assert(buf_len > 0);
   assert(buf_ptr != NULL);
   assert(bytes_read_ptr != NULL);
   assert(eof_ptr != NULL);
   assert(!(*eof_ptr));

   errno = 0;
   bytes_read = read(vars_ptr->input_fd, (void *)buf_ptr, buf_len);
   saved_errno = errno;

   if ( bytes_read == -1 )
   {
      proceed = FALSE;

      fprintf(stderr, "\nCant't read from in input fifo.\n");
      fprintf(stderr, "errno = %d -- \"%s\"\n", saved_errno,
              strerror(saved_errno));
   }
   else if ( bytes_read == 0 )
   {
      *eof_ptr = TRUE;
   }
   else
   {
      assert(bytes_read > 0);
      assert((off_t)bytes_read <= buf_len);

      *bytes_read_ptr = (off_t)bytes_read;
   }

   return(proceed);

} /* load_buf_from_file() */


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
 *		If vars_ptr->output_file_name is not defined, simply
 *		set vars_ptr->output_fd to STDOUT_FILENO.
 *		
 *		If vars_ptr->output_file_name is defined, assert
 *		that vars_ptr->output_fifo_name is undefined.  Then
 *		attempt to open the file O_WRONLY / O_APPEND O_CREAT.
 *		If successful, set vars_ptr->output_fd to the file
 *		descriptor of the newly opened file, and return TRUE.
 *		On failure, generate an error message and return FALSE.
 *
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:	none.
 *
 *-------------------------------------------------------------------------
 */

static int
setup_output_fd(struct fifo2file_vars * vars_ptr)
{
   int proceed = TRUE;
   int fd;
   int result;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);
   assert(vars_ptr->output_fd == -1);

   if ( strlen (vars_ptr->output_file_name) == 0 ) 
   {
      /* send output to stdout */
      vars_ptr->output_fd = STDOUT_FILENO;
   }
   else /* strlen (vars_ptr->output_file_name) > 0 */
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
 *		it referes to a regular file.  In this case we must use 
 *		the close() system call to close it.  Generate an error 
 *		message on failure.
 *
 * Return:      TRUE if successful, and FALSE otherwise.
 *
 * Changes:     none.
 *
 *-------------------------------------------------------------------------
 */
static void
takedown_output_fd(struct fifo2file_vars * vars_ptr)
{
   int result;
   int saved_errno;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);
   assert(vars_ptr->output_fd != -1);

   if ( vars_ptr->output_fd == STDOUT_FILENO )
   {
      vars_ptr->output_fd = -1;
   }
   else /* strlen(vars_ptr->output_file_name) > 0 */
   {
      assert( strlen(vars_ptr->output_file_name) > 0);

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
FIFO 2 File (fifo2file)\n\
\n\
purpose: Read data from the supplied FIFO (AKA named pipe) until it is \n\
   closed and the other end, and write it to the supplied file, or to \n\
   stdout. \n\
\n\
usage: fifo2file [OPTIONS] <input FIFO>\n\
\n\
   OPTIONS\n\
\n\
      -f <file name>     Name of the file to which to write the data \n\
			 read from the FIFO. \n\
\n\
      -v                 Verbosity.  Repeat for more output. \n\
                           verbosity == 0 ==> output on error only \n\
                           verbosity == 1 ==> terse functional output \n\
                           verbosity == 2 ==> verbose functional output \n\
                           verbosity >= 3 ==> debugging output \n\
\n\
   If the -f option is not provided, output is to stdout. \n\n");

   return;

} /* usage() */


/*-------------------------------------------------------------------------
 *
 * Function:    write_buf_to_output_fd()
 *
 * Purpose:     Write the contents of the supplied buffer to the file 
 *		descriptor listed in vars_ptr->output_fd.  This may be 
 *		either sdtoud, or a regular file.
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
                       struct fifo2file_vars * vars_ptr)
{
   char *local_buf_ptr;
   int write_itteration = 0;
   int saved_errno;
   int proceed = TRUE;
   size_t local_buf_len;
   ssize_t bytes_written = 0;
   ssize_t result;

   assert(vars_ptr != NULL);
   assert(vars_ptr->magic == FIFO_2_FILE_VARS_MAGIC);
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
 * Function:     main()
 *
 * Purpose:      main routine for fifo2file.
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
   struct fifo2file_vars vars = 
   {
      /* magic               = */ FIFO_2_FILE_VARS_MAGIC,
      /*                       */ 
      /* input_fifo_name     = */ "",
      /* output_file_name    = */ "",
      /* verbosity           = */ 0,
      /*                       */
      /* input_fd            = */ -1,
      /* output_fd           = */ -1
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

   if ( proceed ) proceed = open_input_fifo(&vars);

   if ( proceed ) proceed = setup_output_fd(&vars);

   if ( proceed ) proceed = copy_data_from_fifo_to_file(&vars);

   if ( vars.input_fd != -1 ) close_input_fifo(&vars);

   if ( vars.output_fd != -1 ) takedown_output_fd(&vars);

   return(!proceed);

} /* main */

