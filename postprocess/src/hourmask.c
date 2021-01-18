# define _XOPEN_SOURCE
# include <stdio.h>
# include <stdlib.h>
# include <time.h>
# include <stdint.h>
# include <string.h>
# include <limits.h>
# include <unistd.h>
# include <math.h>
# include <sys/stat.h>

time_t sec_from_date(char *date);
void date_from_sec(char *date,time_t secs);

int main(int argc, char *argv[])
{
  uint8_t *Zarr;
  long N,Xdim,Ydim,X,Y;
  uint16_t *maskarr,maskval=65535;
  int START=1,NZ,i,len;
  size_t arrsize;
  time_t dT=300,T,endT,startT;
  char buf[500],maskfile[250],*pgmdir,*pgmfile, pgmname[100]={0};
  char *p,*maskname,startstamp[15],endstamp[15],stamp[15]={0},sf[250];
  FILE *ZPGM,*MASKF;

  /*
  sprintf(startstamp,"%s00",argv[1]);
  endstamp=argv[2];

  */
  pgmdir=argv[1];
  pgmfile=argv[2];
  maskname=argv[3];
  if(argc==5) dT=(time_t)atol(argv[4]);

  p=strrchr(pgmfile,'_');
  len=p-pgmfile+1;
  strncpy(pgmname,pgmfile,len);
  sprintf(startstamp,"%.12s00",p+1);
  p=strchr(p,'-');
  sprintf(endstamp,"%.14s",p+1);

  /*  printf("%s %s %s\n",pgmname,startstamp,endstamp); */

  /* dBZ_M00_202011041445-20201104145030.pgm */

  startT=sec_from_date(startstamp);
  endT=sec_from_date(endstamp);

  for(T=startT ; T<endT ; T+=dT)
  {
    date_from_sec(stamp,T);
    sprintf(sf,"%s/%s%.12s-%s.pgm",pgmdir,pgmname,startstamp,stamp);
    ZPGM=NULL;
    ZPGM=fopen(sf,"r");
    if(ZPGM)
    {
       for(i=0;i<3;i++)
       {
         memset(buf,0,500);
         fgets(buf,499,ZPGM);
         if(buf[0]=='#') i--;
         if(i==1) if(START)
         {
           sscanf(buf,"%ld %ld",&Xdim,&Ydim);
           arrsize=Xdim*Ydim;
           Zarr=malloc(arrsize);
           maskarr=calloc(arrsize,2);
           START=0;
         }
       }
       fread(Zarr,arrsize,1,ZPGM);
       fclose(ZPGM);
    }

    for(N=0;N<arrsize;N++)
    {
       NZ = Zarr[N];
       if(NZ==255) maskarr[N]=maskval;
    }

    if((T%3600 == 0) && (T!=startT))
    {
      sprintf(maskfile,"%s/%s_%.12s-%.12s.pgm",pgmdir,maskname,startstamp,stamp);
      MASKF=fopen(maskfile,"w");
      printf("%s\n",maskfile);
      fprintf(MASKF,"P5\n%ld %ld\n65535\n",Xdim,Ydim);
      fwrite(maskarr,arrsize*2,1,MASKF);
      fclose(MASKF);
      memset(maskarr,0,arrsize*2);
    }
  }
  
  free(Zarr);
  free(maskarr);

  return(0);
}

time_t sec_from_date(char *date)
{
   struct tm Sd,*Sdd;
   int y,m;
   time_t secs;
 
   sscanf(date,"%4d%2d%2d%2d%2d%2d",&y,&m,&Sd.tm_mday,&Sd.tm_hour,&Sd.tm_min,&Sd.tm_sec);
   Sd.tm_year=y-1900;
   Sd.tm_mon=m-1;
   Sd.tm_isdst=0;
   Sd.tm_sec=0;
   Sdd=&Sd;
   secs=mktime(Sdd);
   return(secs);
}

void date_from_sec(char *date,time_t secs)
{
   struct tm *Sdd;

   Sdd=gmtime(&secs);
   sprintf(date,"%d%02d%02d%02d%02d%02d",Sdd->tm_year+1900,Sdd->tm_mon+1,
                                         Sdd->tm_mday,Sdd->tm_hour,
                                         Sdd->tm_min,Sdd->tm_sec);
   return;
}
