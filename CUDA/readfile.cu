#include <iostream>
#include "readfile.cuh"
using namespace std;
readfile::readfile()
{
  already_open = 0;
  buffer = new char [MAX_LINE_LENGTH];
  if(!buffer){ cout<< "allocation error in readfile"<<endl; exit(0);}
  result = new char [MAX_LINE_LENGTH];
  if(!result){ cout<<"allocation error in readfile"<<endl; exit(0);}
}

void readfile::openinput(char *file)
{
	
	if(!already_open){
		fd=fopen(file,"r+");
      if(fd==NULL){
         cout<<"readfile::openinput: can't open file"<<endl;
         exit(1);
      }
      already_open=1;
   }
   else{
      rewind(fd);
   }
}

void readfile::closeinput( void )
{
   if(already_open){
     int value=fclose( fd );
     already_open=0;
     if(value!=0)
	 {
		 cout<<"file not correctly closed"<<endl;
		 exit(0);
	 }
   }
}

int readfile::setinput(char *a)
{
   int m,n;

   n = (int)strlen(a);

   rewind(fd);

   while(read_one_line()){
      m=(int)strlen(buffer);
      if(m==n){
          if(strncmp(buffer,a,n)==0)return(1);
      }
   }

   return(0);
}

char* readfile::setget(char *key, char *a)
{
   int m,i,n,j=0;

   n = (int)strlen(a);
                                                   // reset file pointer to the key word
   if (!setinput(key)) {
     printf( "\n readfile::setget: key word '%s' missing\n", key );
     exit(-1);
     }

   while(read_one_line()){                        // read lines following the key
     m=(int)strlen(buffer);
     if (strchr(buffer,38)) break;                // '&' contained -> break
     if(m>n+1) {                                  // length sufficient
       for(int i=0;i<m-n;i++) {                       // scan the line for variable name
	 if(strncmp(buffer+i,a,n)==0){            // if found, write it to result[]
	   if(buffer[n+i]=='='){
	     i++;
	     while(buffer[n+i+j]!=',' && n+i+j<m ) {
	       result[j]=buffer[n+i+j];
	       j++;
	     }
	     result[j]=0;
	     return(result);                        // and return pointer to result
	   }
	 }
       }
     }
   }
   printf(" readfile::setget: can't find name ");   // otherwise: send error message
   for(i=0;i<n;i++)putchar(a[i]);
   printf(" following key word %s \n\n",key);
   exit(1);
   return(result);
}

char* readfile::getinput(char *a)
{
   int m,n,i=0,j=0;

   n = (int)strlen(a);

   rewind(fd);

   while(read_one_line()){                        // read lines
     m=(int)strlen(buffer);
     if(m>n+1) {                                  // length sufficient
       for(i=0;i<m-n;i++) {                       // scan the line for variable name
	 if(strncmp(buffer+i,a,n)==0){            // if found, write it to result[]
	   if(buffer[n+i]=='='){
	     i++;
	     while(buffer[n+i+j]!=',' && n+i+j<m ) {
	       result[j]=buffer[n+i+j];
	       j++;
	     }
	     result[j]=0;
	     return(result);                      // and return pointer to result
	   }
	 }
       }
     }
   }
   printf("readfile::getinput: can't find name ");
   for(i=0;i<n;i++)putchar(a[i]);
   printf(" in input file \n");
   exit(1);
   return(result);
}

int readfile::read_one_line( void )
{
   int i=0,c;
   while(i<MAX_LINE_LENGTH){
      c=getc(fd);
      if(c==EOF)return(0);
      else if(c=='\n'){
         buffer[i++]=0;
         return(1);
      }
      else if(c=='#'){
         buffer[i++]=0;
         while(getc(fd)!='\n');
         return(1);
      }
      else if(c!=' '){
         buffer[i++]=c;
      }
   }
   printf("readfile::read_one_line: line too long\n");
   exit(-1);
   return(-1);
}



