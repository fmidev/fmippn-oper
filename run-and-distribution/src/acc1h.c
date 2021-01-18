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

time_t sec_from_date(char *stim);
char *date_from_sec(char *datestr,time_t secs);

int main(int argc, char *argv[])
{
  uint8_t *Zarr;
  uint16_t *Accarr;
  long N,Xdim,Ydim,X,Y;
  int START=1,accmins,NZ,i;
  size_t arrsize;
  time_t dT=300,T,endT,startT;
  char buf[500],*pgmdir,*pgmname,*endstamp,stamp[13]={0},sf[250],*accfile;
  double *fAccarr,kT,Af,Ac,Bf,Bc,limZ,dBZ,NtoR[256],R,fNoData=655.35;
  FILE *ZPGM,*ACCF;
 
  Af=atof(getenv("ZR_AF"));
  Ac=atof(getenv("ZR_AC"));
  Bf=atof(getenv("ZR_BF"));
  Bc=atof(getenv("ZR_BC"));
  accmins=atoi(getenv("ACCMINS"));
  dT=60*accmins;
  kT = (double)dT/3600.0;

  if(Af==Ac) limZ = 10.0*log10(Af);
  else 
  {
    R = pow(Af/Ac,1/(Bc-Bf));
    limZ = 10.0 * log10(Af*pow(R,Bf));
  }

  /*  printf("# lim dBZ = %.1f\n",limZ); */

  {
     double A,B,C;

     NtoR[0] = NtoR[255] = 0.0;
     A = Af;
     B = Bf;
     C = -log10(A)/B;
     for(N=1;N<255;N++)
     {
        dBZ = 0.5*(double)N - 32.0;
        if(dBZ > limZ) { A=Ac; B=Bc; C = -log10(A)/B; }
        R = pow(10.0, 0.1*dBZ/B + C);
        NtoR[N] = R;
	/*        printf("%.1f\t%.2f\n",dBZ,R); */
     }
  }

   
  endstamp=argv[1];
  pgmdir=argv[2];
  pgmname=argv[3];
  accfile=argv[4];

  endT=sec_from_date(endstamp);
  startT=endT-3600;

  for(T=startT ; T<endT ; T+=dT)
  {
    date_from_sec(stamp,T);
    sprintf(sf,"%s/%s%s",pgmdir,stamp,pgmname);
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
           fAccarr=calloc(arrsize,sizeof(double));
           Accarr=calloc(arrsize,sizeof(uint16_t));
           START=0;
	 }
       }
       fread(Zarr,arrsize,1,ZPGM);
       fclose(ZPGM);
    } else printf("MISSING %s\n",sf);

    for(N=0;N<arrsize;N++)
    {
       NZ = Zarr[N];
       if(NZ<255)
       {
          R = NtoR[NZ];
          fAccarr[N] += R*kT;
       }  else fAccarr[N] = fNoData;
    }
  }
  
  for(N=0;N<arrsize;N++) Accarr[N] = (uint16_t)(100.0 * fAccarr[N]);
  swab(&Accarr[0],&Accarr[0],arrsize*sizeof(uint16_t));

  ACCF=fopen(accfile,"w");
  fprintf(ACCF,"P5\n%ld %ld\n65535\n",Xdim,Ydim);
  fwrite(Accarr,arrsize*2,1,ACCF);
  fclose(ACCF);
  free(Zarr);
  free(fAccarr);
  free(Accarr);
  printf("Generated %s\n",accfile);
  

  return(0);
}

/*---------------------------------------------------------------------------------------*/

time_t sec_from_date(char *stim)
{
   struct tm Sdd;
   time_t secs;
   int res;

   res=sscanf(stim,"%4d%2d%2d%2d%2d%2d",&Sdd.tm_year,&Sdd.tm_mon,
                                    &Sdd.tm_mday,&Sdd.tm_hour,
                                    &Sdd.tm_min,&Sdd.tm_sec);
   Sdd.tm_year-=1900;
   Sdd.tm_mon--;
   if(res<6) Sdd.tm_sec=0;
   secs=mktime(&Sdd);
   return(secs);
}

char *date_from_sec(char *datestr,time_t secs)
{
   struct tm Sdd;

   Sdd=*localtime(&secs);
   strftime(datestr,13,"%Y%m%d%H%M",&Sdd);
   return(datestr);
}
