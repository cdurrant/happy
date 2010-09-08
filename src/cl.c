/*  Last edited: Sep 10 18:26 1996 (rmott) */
/* COMMAND_LINE PARSING AND FILE UTILITY ROUTINES */
/*
   (C) Richard Mott, Genome Analysis Lab, ICRF
   */
#define _GNU_SOURCE

#include<time.h>
#include<sys/time.h>
#include<sys/types.h> 
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"cl.h"

#ifdef ALPHA
extern void exit(int);
#endif

/*command line option to suppress all comments*/
int comment = 1;

FILE *
openfile( char *filename, char *mode )
     
     /* attempts to open filename, using mode to determine function
	(read/write/append), and stops the program on failure */
{
  FILE *fp;
  char *make_legal();

  if ( ( fp = fopen( make_legal(filename), mode ) ) == NULL )
    {
      fprintf( stderr, "\n ERROR: file %s will not open for function %s\n\n", filename, mode );
      exit(1);
    }
  else
    return fp;
}

int 
getfloat( char *format, float *variable, int argc, char **argv ) /* gets a float from command line */
     
/* format can be a string like "-a=%f", or "-a". In the former case is
looks for entries line -a=10.3 AND -a 10.3, in the latter just for -a
10.3 */

{ float t; char *s;
  char Format[256];

  sprintf(Format,"%g",*variable);
  append_usage( format, "float", Format, 0 );

  if ( (s = next_arg( format, argc, argv ) ) && sscanf( s, "%f", &t ) == 1 )
    {
      *variable = t;
      return 1;
    }


  s = format;
  while(*s && *s != '=')
    s++;
  if ( *s == '=' )
    sprintf(Format,"%s", format);
  else
    sprintf(Format,"%s%s", format, "=%f");

  while ( --argc > 0 )
    {
      if ( sscanf( argv[argc], Format, &t ) == 1 ) 
	{
	  *variable = t;
	  return 1;
	}
    }

  return 0;
}

char *
next_arg( char *format, int argc, char **argv )

/* scan the command lines for argument matching stub, and return the
NEXT arg 

eg if format == "-max=%d", then argv == "-max 10" returns "10" */

{
  char *s = cl_stub( format );

  argc--;
  while( --argc > 0 )
    {
      if ( ! strcmp( s, argv[argc] ) )
	return argv[argc+1];
    }
  return NULL;
}

char *
cl_stub( char *format ) /* extract stub from format, eg "-files=%s" -> "-files" */
{
  static char stub[256];
  char *s=stub;

  while ( *format && *format != '=' )
    *s++ = *format++;
  *s = 0;

  return stub;
}


int  
getint( char *format, int *variable, int argc, char **argv ) /* gets an int from the command line */
     
{
  int t;
  char *s;
  char Format[256];

  sprintf(Format,"%d",*variable);
  append_usage( format, "integer", Format, 0 );

  if ( (s = next_arg( format, argc, argv ) ) && sscanf( s, "%d", &t ) == 1 )
    {
      *variable = t;
      return 1;
    }

  s = format;
  while(*s && *s != '=')
    s++;
  if ( *s == '=' )
    sprintf(Format,"%s", format);
  else
    sprintf(Format,"%s%s", format, "=%d");

  while ( --argc > 0 )
    {
      if ( sscanf( argv[argc], Format, &t ) == 1 ) 
	{
	  *variable = t;
	  return 1;
	}
    }

  return 0;
} 


int 
clcheck( char *format, int argc, char **argv ) /* checks if a string is on the command line*/
{
  append_usage( format, "switch", "", 0 );

  while ( --argc > 0 ) 
    {
      if ( strcmp( format, argv[argc] ) == 0 )
	return 1;
    }
	
  return 0;
}

int 
getbool( char *format, int *bool, int argc, char **argv ) 
     
     /* checks for boolean switch, eg
	
	-gap=yes  -gap=y -gap=1 -gap=true -gap=t -gap	TRUE
	-gap=no -gap=n -gap=0 -gap=false -gap=f		FALSE
	
	format should be the string "-gap" 
	
	*/
     
{
  int arg;

  if ( *bool )
    append_usage( format, "switch", "true", 0 );
  else
    append_usage( format, "switch", "false", 0 );

  if ( getint( format, bool, argc, argv) )
    return 1;
  else
    {
      char reply[256];
      memset(reply,0,256);
      if ( getarg( format, reply, argc, argv ) )
	{
	  if ( ! strcasecmp( reply, "yes" ) ||
	      ! strcasecmp( reply, "y" ) ||
	      ! strcasecmp( reply, "true" ) ||
	      ! strcasecmp( reply, "t" )  )
	    {
	      *bool = 1;
	      return 1;
	    }
	  else if ( ! strcasecmp( reply, "no" ) ||
		   !  strcasecmp( reply, "n" ) ||
		   !  strcasecmp( reply, "false" ) ||
		   !  strcasecmp( reply, "f" )  )
	    {
	      *bool = 0;
	      return 1;
	    }
	}
    }

  if ( clcheck( format, argc, argv ) )
    return 1;
  else
    return 0;
}

int 
getboolean( char *format, int *bool, int argc, char **argv )

/* similar to getbool(), except that function returns 1 and sets *bool = 1 if clcheck is true also . Therefore the switch

-switch

is equivalent to

-switch=y

-noswitch is equivalent to 

-switch=n
*/

{
  char noformat[256];
  int Argc = argc;
  int n;

  if ( *bool )
    append_usage( format, "switch", "true", 1 );
  else
    append_usage( format, "switch", "false", 1 );

  if ( *format == '-' )
    sprintf(noformat,"-no%s", &format[1] );
  else
    sprintf(noformat,"-no%s", format );
  
  while ( --Argc > 0 ) 
    {
      if ( ! strcmp( format, argv[Argc] ) )
	{
	  *bool = 1;
	  return 1;
	}
      else if ( ! strcmp( noformat,argv[Argc]) )
	{
	  *bool = 0;
	  return 1;
	}
    }

/* NB there is a slight bug here in that  a command line of the form

   program -notest -test=y

   will yield FALSE

   because it searches for -notest before -test=y 

   THerefore you should stick to one convention within a command-line

*/

  return getbool( format, bool, argc, argv );
}


int 
getarg( char *format, char *arg, int argc, char ** argv ) /* gets a string argument*/
{
  char *s;
  char Format[256];

  append_usage( format, "text", arg, 0 );

  s = format;
  while(*s && *s != '=')
    s++;
  if ( *s == '=' )
    sprintf(Format,"%s", format);
  else
    sprintf(Format,"%s%s", format, "=%s");

  if ( s = next_arg( format, argc, argv ) )
    {
      strcpy( arg, s );
      return 1;
    }
  while ( --argc > 0 ) 
    {
      if ( sscanf( argv[argc], Format, arg ) > 0 ) return 1;
    }
  return 0;
}

int get_restricted_arg( char *format, char *arg, int argc, char ** argv, int stringc, char **strings, int *value ) {
  getarg( format, arg, argc, argv );
  return legal_string( arg, strings, stringc, value );
}

FILE *
nextfile( char *filename, int argc, char **argv )
     
     /* searches the command line for arguments not beginning with "-", and tries
	to open the first one it finds for reading. The corresponding argument in the command-line list is then blanked
	out so that on the next call a different file is opened. The name of the
	file opened is returned in filename. On failure the function returns NULL */
{
  FILE *fp;
  int i;

  *filename = 0;

  while ( --argc > 0 )
    {
      if ( argv[argc][0] != 0 && argv[argc][0] != '-' )
	{
	  if (fp = fopen( filename, "r" ))
	    {
	      strcpy( filename, argv[argc] );
	      for(i=0;i<=strlen(argv[argc]);i++)
		argv[argc][i] = 0;
	      return fp;
	    }
	}
    }

  return NULL;
}

char *
extension ( char *filename, char *ext)
     
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
{
  int n;

  if ( ! ext )
    return NULL;

  if ( ext[0] == '.' )
    ext++;

  n = strlen(filename);

  while ( filename[n] != '.' && n > 0 )
    n--;

  if ( filename[n] != '.' )
    {
      n = strlen(filename);
      filename[n] = '.';
    }

  n++;

  strcpy( &filename[n], ext );

  if ( filename[n=strlen(filename)-1] == '.' )
    filename[n] = 0;

  return filename;
}

char *
my_basename (char *filename)
     
     /* strips any directory from filename
	
	e.g.
	my_basename("/gea/users/rmott/.cshrc");
	
	produces .cshrc
	
	my_basename("/usr");
	
	produces usr
	
	*/
     
{
  int n, m;

  n = strlen(filename);

  while( filename[n] != '/' && n > 0 )
    n--;

  if ( filename[n] == '/' )
    n++;

  m = 0;
  while ( filename[m] )
    filename[m++] = filename[n++];

  return filename;
}

char *
dirname( char *pathname )

/* strips off the filename from the path, leaving the directory */
{
  char *s = &pathname[strlen(pathname)-1];

  while( *s && s > pathname && *s != '/' )
    s--;

  if ( s == pathname && *s != '/' )
    strcpy(pathname, "./");
  else if ( s == pathname && *s == '/' )
    strcpy(pathname,"/");
  else
    *s = 0;

  return pathname;
}

char *
directory( char *filename, char *dir )
     
     /* changes the directory part of filename to dir.
	
	e.g.
	directory( "/home/gea/fred/file.dat", "/usr/local/" );
	
	produces /usr/local/file.dat
	
	directory("/usr/lib", "var");
	
	produces var/lib, i.e. a directory slash is inserted if necessary.
	
	*/
     
{
  char name[512];
  
  my_basename( filename );
  
  if ( dir && *dir  )
    {
      if ( dir[strlen(dir)-1] != '/' )
	sprintf( name, "%s/%s", dir, filename );
      else
	sprintf( name, "%s%s", dir, filename );

      strcpy( filename, name );
    }
  return filename;
}	

char *
rootname( char *filename )
     
     /* trims off the directory and the extension from filename */
     
{
  return extension(my_basename(filename),"");
}


/* use compiler switch -DVAX to use the form of make_legal which converts . to _ */

#ifdef VAX

char *
make_legal( char *filename )
     
     /* converts files like file.in.out to file_in.out */
     
{
  char *s = filename, *t;
  int n=0;

  while( *s )
    {
      if ( *s == '.' )
	{
	  n++;
	  t = s;
	}
      s++;
    }

  if ( n > 1 )
    {
      t--;
      while ( t != filename )
	{
	  if ( *t == '.' )
	    *t = '_';
	  t--;
	}
    }
  return filename;
}


#else

char *
make_legal( char *filename )
     
{
  return filename;
}

#endif

FILE *
argfile( char *format, char *mode, int argc, char **argv, char *filename )
     
     /* parses the command-line using format, attempts to get a filename and open
	it. If format does not match any argument then NULL is returned. If format
	matchs an argument but the file will not open then an error message is printed to stderr */
     
{
  FILE *fp=NULL;

  *filename = 0;
  if ( getarg( format, filename, argc, argv ) )
    {
      if ( ( fp = fopen( filename, mode ) ) == NULL )
	fprintf( stderr, "\n ERROR: file %s will not open for function %s\n\n", filename, mode );
    }
  if ( !strcmp( mode, "w" ) )
    append_usage( format, "Writeable File", "", 1 );
  else if ( !strcmp( mode, "r" ) )
    append_usage( format, "Readable File", "", 1 );
  else
    append_usage( format, "File", "", 1 );

  return fp;
}


FILE *
argfile_force( char *format, char *mode, int argc, char **argv, char *filename )
     
     /* parses the command-line using format, attempts to get a filename and open
	it. If format does not match any argument or if format matchs an argument
	but the file will not open then the program stops*/
     
{
  FILE *fp;

  if ( getarg( format, filename, argc, argv ) )
    {
      if ( ( fp = fopen( filename, mode ) ) == NULL )
	{
	  fprintf( stderr, "\n ERROR: file %s will not open for function %s\n\n", filename, mode );
	  exit(1);
	}
      else
	return fp;
    }
  else
    {
      fprintf( stderr, "\n ERROR: argument %s not on command line\n\n", format );
      exit(1);
    }
}


#include<sys/types.h>
#include<sys/stat.h>
#include<sys/times.h>


int 
file_time( char *filename )
     
     /* returns the time that filename was last modified, or 0 if filename
	cannot be opened for reading or if the call to stat fails */
     
{

  FILE *fp;
  struct stat buf;

  if ( ! (fp = fopen(filename,"r")) )
    return 0;
  else
    {
      fclose(fp);
      if ( ! stat( filename, &buf ) )
	return (int)buf.st_mtime;
      else
	return 0;
    }
}

char *
file_date( char *filename )
     
     /* similar to file_time() execpt that a pointer to the date (in English) */

{
  FILE *fp;
  struct stat buf;
  char date[256], *t;

  strcpy( date, "?");

  if ( (fp = fopen(filename,"r")) )
    {
      fclose(fp);
      if ( ! stat( filename, &buf ) )
	{
	  sprintf( date, "%s", ctime(&buf.st_mtime) );
	  t = date;
	  while ( *t )
	    {
	      if ( *t == '\n')
		*t = 0;
	      t++;
	    }

	}
    }
  t = date;
  return t;
}

typedef struct use_list
{
  char *format;
  char *text;
  char *def;
  struct use_list *next;
}
USE_LIST;

static USE_LIST *use_begin=NULL, *use_end=NULL;

void 
append_usage( char *format, char *text, char *def, int overide )
     
{
  USE_LIST *use;

  if ( ! use_begin )
    use_begin = use_end = (USE_LIST *)calloc( 1, sizeof(USE_LIST) ); 
  else 
    {
      use = use_begin;
      while( use )
	{
	  if  ( ! strcmp( use->format, format ) )
	    {
	      if ( overide )
		break;
	      else
		return;
	    }
	  use = use->next;
	}

      if ( ! use )
	use_end = use_end->next = (USE_LIST *)calloc( 1, sizeof(USE_LIST) );
      else
	use_end = use;
    }

/*  printf("append %s %s %s\n", format, text, def ); */
  use_end->format = (char*)strdup( format );  
  use_end->text = (char*)strdup( text );
  use_end->def = (char*)strdup( def );
}

void
gethelp( int argc, char **argv )
{

  int want_help = clcheck("-help", argc, argv );
  int errors = check_usage( stderr, argc, argv );
  if ( want_help || errors )
    print_usage( argc, argv, 1 );
}

void 
print_usage( int argc, char **argv, int stop )
     
{
  USE_LIST *u;

  u = use_begin;

  fprintf( stderr, "\nusage: %s\n", argv[0] );
  while( u )
    {
      fprintf( stderr, "	%-15s %-20s [ %s ]\n", u->format, u->text, u->def );
      u = u->next;
    }
  fprintf(stderr, "\n\n");
  if ( stop ) 
    exit(1);
}

int 
check_usage( FILE *fp, int argc, char **argv )
{
  USE_LIST *u;
  int errors=0;
  char *arg, *noarg;

  while( --argc > 0 )
    {
      arg = argv[argc];
      if ( *arg == '-' )
	{
	  if ( strlen(arg) > 3  && (arg[1] == 'n') && (arg[2] == 'o') )
	    noarg = &arg[3];
	  else 
	    noarg = NULL;

	  u = use_begin;
	  while( u && strncmp( argv[argc], u->format, strlen(argv[argc]) ) && ( noarg == NULL || strncmp( noarg, &u->format[1], strlen(noarg) ) ) )
	    u = u->next;

	  if ( ! u && ! isdigit(arg[1]) )
	    {
	      if ( fp )
		fprintf(fp, "WARNING: unknown argument %s\n", argv[argc] );
	      errors++;
	    }
	}
    }
  
  return errors;
}

#define COM "-comfile="
#define MAX_COM_DEPTH 10
static int com_depth=0;

int 
add_commands_from_file( int argc, char **argv, int *nargc, char ***nargv )
     
     /* looks for a command file argument and parses the command file, adding at to the command line */
     
{
  FILE *fp;
  char comfile[256], line[256], *s, **targv = argv;
  int n = argc, targc=argc;

  *nargc = targc = argc;
  *nargv = targv = argv;

  if ( ++com_depth < MAX_COM_DEPTH  && (fp = argfile( "-comfile=%s", "r", argc, argv, comfile ) ))
    {
      /*      printf( "! %2d reading commands from file %s\n",  com_depth, comfile );*/

      while( skip_comments( fp, line ) != EOF  )
	targc++;

      rewind(fp);

      targv = (char**)calloc( targc, sizeof(char*) );

      while( n-- )
	{
	  if ( argv[n] && strncmp( COM, argv[n] , strlen(COM) ) )
	    targv[n] = argv[n];
	  else targv[n] = (char*)strdup(" ");
	}
      targc = argc;
      while( skip_comments( fp, s=line ) != EOF )
	{
	  while( isspace(*s) )
	    *s++;
	  targv[targc++] = (char*)strdup(s);
	}
      add_commands_from_file( targc, targv, nargc, nargv );
      com_depth--;
      return 1;
    }
  else
    {
      com_depth--;
      return 0;
    }
}

int 
legal_string( char *string, char **strings, int size, int *value )

/* checks if string is a member os strings, and sets value to the index in the array strings
returns 1 on success and 0 on failure */
{
  int i;

  if ( string )
    for(i=0;i<size;i++)
      if ( ! strcmp( string, strings[i] ) )
	{
	  *value = i;
	  return 1;
	}

  return 0;
}

	

char **
split_on_separator( char *string, char separator, int *items)

/* split a string into substrings with separator 
Returns array of null-terminated strings, or null if none found
items returns the number of strings;
*/

{
  char **item=NULL;
  char *s, *t;
  int len;

  *items = 0;

  if ( string )
    {
      *items = 1;
      s = string;
      while(*s)
	{
	  (*items) += ( *s == separator );
	  s++;
	}

      item = (char**)calloc(*items,sizeof(char*));
      s = string;
      *items = 0;
      while( *s )
	{
	  t = s;
	  while( *s != separator && *s )
	    s++;
	  
	  if ( s != t )
	    {
	      len = s-t;
	      item[*items] = (char*)calloc(len+1,sizeof(char));
	      strncpy(item[*items],t,len);
	      (*items) ++;
	    }
	  if ( *s == separator )
	    s++;
	}
    }
  return item;
}

char *Now(void) /* returns a character string containing the current date and time */
{

	time_t time(), *tloc;
	char *ctime(), *t;
	int i;
	static char temp[256];

	tloc = (time_t*)calloc(1,sizeof(time_t));
	time(tloc);
	sprintf( temp, "%s", ctime(tloc) );

	t = temp;
	while ( *t )
	{
		if ( *t == '\n')
			*t = 0;
		t++;
	}

	return temp;
}

int EndsWith( char *text, char *suffix ) {

  int n = strlen(text)-1;
  int m = strlen(suffix)-1;

  if ( m > n )
    return 0;
  else {
    char *t, *s;
    while( m >= 0 ) {
      if ( text[n] != suffix[m] )
	return 0;
      n--; m--;
    }
  }
  return 1;
}
