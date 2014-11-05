#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <math.h>



int main()
{
  int fid, n=6, m=4;
  double key[n], value[m];



  if((fid=open("cache_file", O_RDONLY))==-1){
    perror("Open cache file error\n");
    return 1;
  }
  
  while(read(fid,key, n*sizeof(double))==n*sizeof(double)){
    read(fid,value, m*sizeof(double));
    printf("x=[%lf,%lf,%lf,%lf,%lf,%lf]-> [%lf,%lf,%lf,%lf]\n",key[0],
	   key[1], key[2], key[3], key[4], key[5],
	   value[0], value[1], value[2], value[3]);
  }

  close(fid);

  return 0;
}
