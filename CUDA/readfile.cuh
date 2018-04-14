#define READFILE_H

#define MAX_LINE_LENGTH 1000
#define MAX_COL 50
///////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//------------------------------------------------
#ifndef READ
#define READ
class readfile{

public:
	readfile();
	void   openinput(char*);
	void   closeinput( void );
	int    setinput( char* );
    char*  getinput( char* );
    char*  setget( char*, char* );
	int    read_one_line( void );
  
private:
	  int  already_open;
	  FILE *fd;
	  char *buffer;
	  char *result;

};


#endif
