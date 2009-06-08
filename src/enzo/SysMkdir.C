// mkdir system call - George Kremenek / SDSC SRB group
//   default umask is 022 (mode 755)
//   leading slash is required!
 
 
#include <sys/stat.h>
#include <stdio.h>
#include <errno.h>
#include <string.h>
 
#define MAX_TOKEN 512
 
int SysMkdir (char *startDir, char *directory)
{
    int status;
    int startLen;
    int pathLen, tmpLen;
    char tmpPath[MAX_TOKEN];
    struct stat statbuf;
 
    startLen = strlen (startDir);
    pathLen = strlen (directory);
 
    strcpy (tmpPath, directory);
 
    tmpLen = pathLen;
 
    while (tmpLen > startLen) {
        status = stat (tmpPath, &statbuf);
        if (status == 0) {
            if (statbuf.st_mode & S_IFDIR) {
                break;
            } else {
                fprintf(stderr, "SysMkdir: A local non-directory %s already exists \n", tmpPath);
                return (1);
            }
        }
 
        // Go backward
 
        while (tmpLen && tmpPath[tmpLen] != '/')
            tmpLen --;
        tmpPath[tmpLen] = '\0';
    }
 
    // Go forward and mk the required dir
 
    while (tmpLen < pathLen) {
 
        // Put back the '/'
 
        tmpPath[tmpLen] = '/';
 
        status = mkdir(tmpPath, 0755);
 
        if (status < 0) {
            fprintf (stderr, "SysMkdir: mkdir failed for %s, errno = %d\n", 
		     tmpPath, errno);
            return (-1);
        }
 
        while (tmpLen && tmpPath[tmpLen] != '\0')
            tmpLen ++;
 
    }
 
    return (0);
 
}
 
/*
int main ()
{
  int ierr;
 
  ierr=SysMkdir ("/tmp", "/tmp/a/b/c/d");
  printf("ierr=%"ISYM", dir=(/tmp/a/b/c/d)\n",  ierr);
  return(ierr);
}
*/
