/* Converts values between arg3 and arg4 to arg5 from source 16-bit PGM file arg1 to tarǵet file arg2 */

# define _BSD_SOURCE
# define _XOPEN_SOURCE
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include <stdint.h>
# include <unistd.h>

int main(int argc, char *argv[])
{
   int Xdim,Ydim,i;
   char hdr[300];
   long N;
   size_t arrsize; 
   uint16_t *pgmarr, val, repval, minval, maxval;
   FILE *INPGM, *OUTPGM;

   INPGM=fopen(argv[1],"r");
   OUTPGM=fopen(argv[2],"w");
   minval=atoi(argv[3]);
   maxval=atoi(argv[4]);
   repval=atoi(argv[5]);

  for(i=0;i<3;i++)
  {
     memset(hdr,0,300);
     fgets(hdr,299,INPGM);
     fprintf(OUTPGM,"%s",hdr);
     if(hdr[0]=='#') { i--; continue; }
     if(i==1) sscanf(hdr,"%d %d",&Xdim,&Ydim);
  }
  arrsize=Xdim*Ydim;
  pgmarr=malloc(arrsize*2);
  fread(pgmarr,1,arrsize*2,INPGM);
  swab(pgmarr,pgmarr,arrsize*2);
  fclose(INPGM);

  for(N=0;N<arrsize;N++)
  {
     val=pgmarr[N];    
     if(val>=minval && val<=maxval) pgmarr[N]=repval;
  }
  
  swab(pgmarr,pgmarr,arrsize*2);
  fwrite(pgmarr,1,2*arrsize,OUTPGM);
  fclose(OUTPGM);
  free(pgmarr);
  return(0);
}

