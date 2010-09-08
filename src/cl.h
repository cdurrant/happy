/*  Last edited: Feb 13 13:49 1996 (rmott) */
/* CL.H is the header file for the command-line parsing functions in CL.C */
/* (C) Richard Mott, ICRF */


FILE *openfile( char *filename, char *mode );
     
     /* attempts to open filename, using mode to determine function
	(read/write/append), and stops the program on failure */

FILE *openfile_in_searchpath(char *base, char *mode, char *searchpath, char *fullname );
FILE *openfile_in_envpath(char *base, char *mode, char *env, char *fullname );

int skip_comments( FILE *file, char *string );
int getfloat( char *format, float *variable, int argc, char **argv ); /* gets a float from command line */
int getint( char *format, int *variable, int argc, char **argv ); /* gets an int from the command line */
int clcheck( char *format, int argc, char **argv ); /* checks if a string is on the command line*/
int getbool( char *format, int *bool, int argc, char **argv );
     
     /* checks for boolean switch, eg
	
	-gap=yes  -gap=y -gap=1 -gap=true -gap=t -gap	TRUE
	-gap=no -gap=n -gap=0 -gap=false -gap=f		FALSE
	
	format should be the string "-gap" 
	
	*/
int getboolean( char *format, int *bool, int argc, char **argv ); /* similar to getbool(), except that function returns 1 and sets *bool = 1 if clcheck is true also . Therefore the switch

-switch

is equivalent to

-switch=y

*/
int getarg( char *format, char *arg, int argc, char **argv ); /* gets a string argument*/

FILE *nextfile( char *filename, int argc, char **argv );
     
     /* searches the command line for arguments not beginning with "-", and tries
	to open the first one it finds for reading. The corresponding argument in the command-line list is then blanked
	out so that on the next call a different file is opened. The name of the
	file opened is returned in filename. On failure the function returns NULL */

char *extension ( char *filename, char *ext);
     
     /* changes the extension of filename to ext.
	
	e.g.
	extension( "file.dat", "out" );
	
	produces file.out
	
	extension( "file", "out");
	
	produces file.out
	
	extension( "file", ".out");
	
	produces file.out
	
	extension( "file.out", "");
	
	produces file
	
	Note: filename MUST be long enough to cope ! 

        Returns the modified string */


char *my_basename (char *filename);
     
     /* strips any directory from filename
	
	e.g.
	basename("/gea/users/rmott/.cshrc");
	
	produces .cshrc
	
	basename("/usr");
	
	produces usr
	
	*/
     
char *dirname( char *pathname );

/* strips off the filename from the path, leaving the directory */

char *directory( char *filename, char *dir );
     
     /* changes the directory part of filename to dir.
	
	e.g.
	directory( "/home/gea/fred/file.dat", "/usr/local/" );
	
	produces /usr/local/file.dat
	
	directory("/usr/lib", "var");
	
	produces var/lib, i.e. a directory slash is inserted if necessary.
	
	*/
     
char *rootname( char *filename );
     
     /* trims off the directory and the extension from filename */
     
char *make_legal( char *filename );

FILE *argfile( char *format, char *mode, int argc, char **argv, char *filename );
     
     /* parses the command-line using format, attempts to get a filename and open
	it. If format does not match any argument then NULL is returned. If format
	matchs an argument but the file will not open then the program stops*/
     
FILE *argfile_force( char *format, char *mode, int argc, char **argv, char *filename );
     
     /* parses the command-line using format, attempts to get a filename and open
	it. If format does not match any argument or if format matchs an argument
	but the file will not open then the program stops*/

int file_time( char *filename );
     
     /* returns the time that filename was last modified, or 0 if filename
	cannot be opened for reading or if the call to stat fails */
     
char *file_date( char *filename );
     
     /* similar to file_time() execpt that a pointer to the date (in English) */

int legal_string( char *string, char **strings, int size, int *value );

/* checks if string is a member os strings, and sets value to the index in the array strings
returns 1 on success and 0 on failure */

char *next_arg( char *format, int argc, char **argv );
char *cl_stub( char *format );

int add_commands_from_file( int argc, char **argv, int *nargc, char ***nargv );
void print_usage( int argc, char **argv, int stop );
void append_usage( char *format, char *text, char *def, int overide );
void gethelp( int argc, char **argv );
int check_usage( FILE *fp, int argc, char **argv );


char 
**split_on_separator( char *string, char separator, int *items);

int get_restricted_arg( char *format, char *arg, int argc, char ** argv, int stringc, char **strings, int *value );

char *Now(void);
int EndsWith( char *text, char *suffix );

