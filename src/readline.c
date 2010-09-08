/*  Last edited: Oct 18 14:40 1995 (rmott) */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include"readline.h"

/* functions for reading in lines (C) Richard Mott,, ICRF */


void uncomment( char *string ) /* truncates string at the first ! */
{
  while ( *string != '!' && *string != '#' && *string != 0 ) string++;
  *string = 0;
}

int read_line( FILE *file, char *string ) {
  int	c;
  int	i=0;
	
  if ( file != NULL )  {
    while(c=getc(file))  {
      if (!i && c==EOF)  return EOF;
      if (i && c==EOF ) return i;
      if (c=='\n') return i;
      string[i] = c;
      string[++i] = '\0';
    }
  }
  return EOF;
}

int next_line( FILE *file ) {
  int c;
  if ( file != NULL ) {
    while(c=getc(file)) {
      if (feof(file)) return 0;
      if (c=='\n') return 1;
    }
  }
  return EOF;
}

int not_blank( char *string ) /* checks whether string is full of white space */
{
  while ( *string != 0 ) {
    if ( ! isspace(*string) ) return 1;
    string++;
  }
    
  return 0;
}

/* reads in successive lines, truncating
   comments and skipping blank lines */

int skip_comments( FILE *file, char *string ) { 
  int n = EOF;
  
  *string = 0;
  
  if  ( file ) {
    while ( ( n = read_line( file, string ) ) != EOF ) {
      uncomment( string );
      if ( not_blank( string ) ) return n;
    }
  }
  return n;
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

