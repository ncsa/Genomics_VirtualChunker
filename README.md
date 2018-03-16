# Genomics_VirtualChunker

Virtual Chunker

purpose: Extract the specified chunk from the supplied FASTQ file, and
   write it to stdout, a file, or a FIFO (AKA named pipe) as specified.

usage: vc [OPTIONS] <input FASTQ file>

   OPTIONS
      -d                 Dry run.  Compute all chunk boundaries.
                         Display if verbosity > 0. No output to file.

      -s n[k|m|g]        Chunk size in kilo-bytes (if "k" suffix appears),
                         mega-bytes (if no suffix, or "m" suffix), or
                         giga-bytes (if "g" suffix appears).

      -i n               Index of the chunk to extract.
                         Ignored if dry run is specified.

      -m n               Max chunk index.

      -f <file name>     Name of the file to which to write the chunk.

      -p <FIFO name>     Name of the FIFO (AKA named pipe) to which to
                         write the chunk.

      -v                 Verbosity.  Repeat for more output.
                           verbosity == 0 ==> output on error only
                           verbosity == 1 ==> terse functional output
                           verbosity == 2 ==> verbose functional output
                           verbosity >= 3 ==> debugging output

   If chunk index is equal to max chunk index, then the delivered chunk
   will extend to the end of file, regardless of size.

   The -f and -p options are mutually exclusive.  If neither is provided,
   output is to stdout.
