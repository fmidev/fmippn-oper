/* 

Converts radar data projections and geometry. 
Accepted input formats are cartesian IRIS, ODIM HDF5, PGM and PPM.
Accepted output formats are GeoTIFF, ODIM HDF5, PPM and PGM.
20180524 / Harri Hohti, FMI                                  

*/

# define _BSD_SOURCE
# include <stdio.h>
# include <stdlib.h>
# include <endian.h>
# include <hdf5.h>
# include <hdf5_hl.h>
# include <projects.h>
# include <string.h>
# include <unistd.h>
# include <getopt.h>
# include <stdint.h>
# include <ctype.h>
# include <inttypes.h>
# include <sys/types.h>
# include <sys/stat.h>
# include <time.h>
# ifndef NO_GEOTIFF
# include <gdal.h>
# include <cpl_string.h>
# include <ogr_srs_api.h>
# endif

# define IN 0
# define OUT 1
# define GEO 0
# define CART 1
# define SW 0
# define NE 1
# define CE 2
# define NW 3
# define SE 4
# define XD 0
# define YD 1
# define KOR 0
# define VAN 1
# define ANJ 2
# define IKA 3
# define KUO 4
# define VIM 5
# define UTA 6
# define LUO 7
# define KES 8
# define PET 9
# define RADN 10

# define P_aeqd  0
# define P_merc  1
# define P_stere 2
# define P_utm   3
# define P_nsper 4
# define P_eqc   5
# define P_gnom  6
# define P_tmerc 7
# define P_lcc   8

# define IRIS_CAPPI 3
# define IRIS_TOPS  5

# define ODIM_TH        0
# define ODIM_DBZH      1
# define ODIM_VRAD      2
# define ODIM_HCLASS    3
# define ODIM_ETOP      4
# define ODIM_RATE      5
# define ODIM_ACRR      6
# define ODIM_DBZHC     7
# define ODIM_HAILPROB  8
# define ODIM_SNOWPROB  9
# define ODIM_HCLASSBG 10
# define ODIM_ACCPROB  11
# define ODIM_THEIGHT  12
# define ODIM_PHEIGHT  13
# define ODIM_DBZH_QC  14
# define ODIM_AMVU     15
# define ODIM_AMVV     16
# define ODIM_AMVEL    17
# define ODIM_AMDIR    18
# define ODIM_AMVELDEV 19
# define ODIM_AMDIRDEV 20


typedef struct { char ODIM_qstr[10]; char GTIFF_qstr[100]; char unit[20]; } QuantityStrings;
typedef struct { int has_arg; const char *name; const char *expl; } OwnOption;

static QuantityStrings QuantStr[] = 
{
  {"TH",       "Total reflectivity",                           "dBZ"   },
  {"DBZH",     "Corrected reflectivity",                       "dBZ"   },
  {"VRAD",     "Radial velocity",                              "m/s"   },
  {"HCLASS",   "Hydrometeor classification",                   "class" },
  {"ETOP",     "Echo top height",                              "km"    },
  {"RATE",     "Precipitation intensity",                      "mm/h"  },
  {"ACRR",     "Precipitation accumulation",                   "mm"    },
  {"DBZHC",    "Attenuation corrected reflectivity",           "dBZ"   },
  {"HAILPROB", "Probability of hail",                          "%"     },
  {"SNOWPROB", "Probability of snow",                          "%"     },
  {"HCLASSBG", "Hydrometeor classification background",        "class" },
  {"ACCPROB",  "Probability of accumulation",                  "%"     },
  {"THEIGHT",  "Height of isotherm",                           "m"     },
  {"PHEIGHT",  "Height of isobar",                             "m"     },
  {"DBZH_QC",  "Quality controlled reflectivity",              "dBZ"   },
  {"AMVU",     "U-component of atmospheric motion vector",     "m/s"   },
  {"AMVV",     "V-component of atmospheric motion vector",     "m/s"   },
  {"AMVEL",    "Velocity of atmospheric motion vector",        "m/s"   },
  {"AMDIR",    "Direction of atmospheric motion vector",       "deg"   },
  {"AMVELDEV", "Deviation of velocity of AMV",                 "m/s"   },
  {"AMDIRDEV", "Deviation of direction of AMV",                "deg"   },
  {"\0","\0","\0"}
};

static OwnOption optionlist[] =
{
     {1, "cfgfile",            "Configuration file. Can contain any of following options (without preceding --)." },
     {1, "LUTdir",             "Directory to store geometry conversion LUTs. Default is working directory." },
     {1, "inSWlon",            "SW corner longitude of input data" },
     {1, "inSWlat",            "SW corner latitude of input data" },
     {1, "inNElon",            "NE corner longitude of input data" },
     {1, "inNElat",            "NE corner latitude of input data" },
     {1, "inCElon",            "Center longitude of input data" },
     {1, "inCElat",            "Center latitude of input data" },
     {1, "inBBOXgeo",          "Bounding box of input image in geographic coordinates (CSV format: SWlon,SWlat,NElon,NElat)" },
     {1, "lon_0",              "Reference longitude (if not defined otherwise)" },
     {1, "lat_0",              "Reference latitude (if not defined otherwise)" },
     {1, "indim",              "Common X and Y dimension in pixels of input image (square image)" },
     {1, "inXdim",             "X dimension of input image" },
     {1, "inYdim",             "Y dimension of input image" },
     {1, "inres",              "Common resolution of input image in meters/pixel (square image)" },
     {1, "inXres",             "X resolution of input image in meters/pixel or degrees/pixel if lonlat projection" },
     {1, "inYres",             "Y resolution of input image in meters/pixel or degrees/pixel if lonlat projection" },
     {1, "inrange",            "Common range of input image in kilometers (metric dimension of image)" },
     {1, "inXrange",           "X range of input image in kilometers along line going thru reference point" },
     {1, "inYrange",           "Y range of input image in kilometers" },
     {1, "inproj",             "PROJ4 string of input projection. Format +proj=projname +lon_0=34.3 etc. Overrides possible metadata found in input file." },
     {1, "infile",             "Input file. IRIS, PGM, PPM and ODIM HDF5 are accepted. Can be given also as the first non-option argument." },
     {1, "outSWlon",           "SW corner longitude of output data" },
     {1, "outSWlat",           "SW corner latitude of output data" },
     {1, "outNElon",           "NE corner longitude of output data" },
     {1, "outNElat",           "NE corner latitude of output data" },
     {1, "outCElon",           "Center longitude of output data" },
     {1, "outCElat",           "Center latitude of output data" },
     {1, "outBBOXgeo",         "Bounding box of output image in geographic coordinates (CSV format: SWlon,SWlat,NElon,NElat)" },
     {1, "outBBOXcart",        "Bounding box of output image in cartesian coordinates of output projection (CSV format: SWx,SWy,NEx,NEy)" },
     {1, "outdim",             "Common X and Y dimension in pixels of output image (square image)" },
     {1, "outXdim",            "X dimension of output image" },
     {1, "outYdim",            "Y dimension of output image" },
     {1, "outres",             "Common resolution of output image in meters/pixel (square image)" },
     {1, "outXres",            "X resolution of output image in meters/pixel" },
     {1, "outYres",            "Y resolution of output image in meters/pixel" },
     {1, "outrange",           "Common range of output image in kilometers (metric dimension of image)" },
     {1, "outXrange",          "X range of output image in kilometers along line going thru reference point" },
     {1, "outYrange",          "Y range of output image in kilometers" },
     {1, "outproj",            "PROJ4 string of output projection. Format +proj=projname +lon_0=34.3 etc." },
     {1, "outfile",            "Output file. Possible formats are PGM, PPM, ODIM HDF5 and GeoTIFF. Can also be the second non-option argument." },
     {1, "outbytes",           "Forced output data quantity resolution in bytes. Usually taken from input data." },
     {1, "ODIMsource",         "ODIM source string, e.g. originating centre, radar node code etc." },
     {1, "ODIMquantity",       "Radar quantity as ODIM code" },
     {1, "ODIMdatasetnum",     "ODIM dataset number (e.g. sweep number or cartesian field)" },
     {1, "ODIMdatanum",        "ODIM data number (e.g. quantity number)" },
     {1, "ODIMobject",         "ODIM object type (defaut is IMAGE, other are POLAR or COMP)" },
     {1, "ODIMproduct",        "ODIM product type" },
     {1, "ODIMACCnum",         "ODIM number of images used in accumulation product" },
     {1, "sweepnumber",        "Sweep number of IRIS RAW input data" },
     {1, "polarzoom",          "Intermediate zoom factor when converting polar data to cartesian" },
     {1, "polarrange",         "Maximum range in kilometers from radar when converting polar data to cartesian" },
     {1, "EPSG",               "EPSG code for output projection" },
     {1, "TIFFcompression",    "Compression method for GeoTIFF output file. Default is DEFLATE." },
     {1, "Zlevel",             "Compression level for GeoTIFF output file, from 0 (no) to 9 (best, slowest). Default is 6." },
     {1, "tiling",             "Sets GeoTIFF tiling. Default is 1 (YES). If tiling>1, the X and Y block dimensions are set to the value." },
     {1, "Xtile",              "X dimension of GeoTIFF tile (default is 256)." },
     {1, "Ytile",              "Y dimension of GeoTIFF tile (default is 256)." },
     {1, "obstime",            "Observation time of data if not in input files. Format YYYYMMDDhhmm." },
     {1, "overviews",          "Set overviews to GeoTIFF. Format: overviews=CUBIC,2,4,8,16,32 (default). If downsampling method is not known NONE is used.\n                        Available methods are NEAREST, GAUSS, CUBIC, AVERAGE, MODE, AVERAGE_MAGPHASE or NONE. overviews=OFF disables overviews."},
     {1, "starttimes",         "Time stamp of beginning of the polar scan(s) in HDF5 output file. Format YYYYMMDDhhmmss." },
     {1, "endtimes",           "Time stamp of end of the polar scan(s) in HDF5 output file. Format YYYYMMDDhhmmss." },
     {1, "metafile",           "Contents of this file will be inserted to GeoTIFF GDALMetadata tag. Format should be XML" },
     {1, "QDparam",            "Quantity parameter for QD formatted PGM header, e.g. PrecipitationRate, CorrectedReflectivity." },
     {0, "QDheaders",          "Flag for printing QD header information to PGM header output file as comments" },
     {0, "RtodBZ",             "Flag for converting input data from rain rate to dBZ using Z=AR^B relation" },
     {0, "dBZtoR",             "Flag for converting input data from dBZ to rain rate using Z=AR^B relation" },
     {1, "ZRA",                "Parameter A of Z=AR^B relation" },
     {1, "ZRB",                "Parameter B of Z=AR^B relation" },
     {1, "ZSA",                "Parameter A of Z=AS^B relation for snow" },
     {1, "ZSB",                "Parameter B of Z=AS^B relation for snow" },
     {1, "ZRAc",               "Parameter A of Z=AR^B relation for convective rain" },
     {1, "ZRBc",               "Parameter B of Z=AR^B relation for convective rain" },
     {1, "Rfactor",            "Factor to multiply R in R(Z) conversion when adjustment with rain gauge is needed" },
     {1, "isotherm",           "Isotherm temperature for THEIGHT quantity [C]"},
     {1, "isobar",             "Isobar pressure for PHEIGHT quantity [hPa]"},
     {0, "nodefaultmeta",      "Flag for not to insert default GDALMetadata tag" },
     {0, "use_incorners",      "Flag for using corner coordinates of input data also for output data" },
     {0, "use_inSW",           "Flag for using SW corner coordinates of input data also for output data" },
     {0, "use_inNE",           "Flag for using NE corner coordinates of input data also for output data" },
     {0, "ignore_incenter",    "Flag for ignoring center coordinates of polar coordinate input data" },
     {0, "use_incenter",       "Flag for using center coordinates of input data to define geometry of output data" },
     {0, "flip_output",        "Flag for flipping output data over X-axis" },
     {0, "fill_nodata",        "Flag for replacing nodata values with undetect value in polar data or during cartesian dBZ to rate data conversion" },
     {0, "overwrite",          "Flag for overwriting existing GeoTIFF output file" },
     {0, "IRISrain",           "Flag for converting rain rate or accumulation data to IRIS 16-bit float data type." },
     {0, "verbose",            "Flag for printing metadata information and messages during execution. This should be the first option." },
     {0, "quiet",              "Flag for no stdout. This should be the first option." },
     {0, "help",               "" },
     {-1,"",""}
};

static struct option long_options[300],null_option={0,0,0,0}; 

char   projtypes[9][20]={"aeqd","merc","stere","utm","nsper","eqc","gnom","tmerc","lcc"};
char   *projtype,InLUT[300],OutLUT[300];

projUV coords[2][2][5]={{{{0,0}}}},
       DEG[5],
       IRISgeo;
 
double inrange[2]={0},
       outrange[2]={0},
       outres[2]={0},
       inres[2]={0},
       inoffs[2]={0},
       polarzoom=1,
       polarrange=-1,
       ERad,
       ERflatten,
       lon_0=190,
       lat_0=190,
       binlen,
       elangle,
       NI=0.0,
       TOPthr=0.0,
       Goffset=0.0,
       Ggain=1.0,
       inUndetect, /* data transformation parameters of input HDF5 data */
       inNodata,
       inGain,
       inOffset,
       ZRA=200.0,ZRB=1.6,ZRB10, /* Default (Marshall-Palmer) A and B in Z=AR^B */
       ZSA=0.0,ZSB=0.0,ZSB10, /* A and B in Z=AS^B in snow */
       ZRC, /* constant for R(Z) conversion,  -log10(ZRA)/ZRB */
       ZRAc=0,ZRBc=0, ZRB10c, ZRCc, /* convective A and B. Same as frontal If not explicitely set */
       conv_dBZlim=500.1, /* dBZ value where rain changes from frontal to convective. See calc_convlim. */
       conv_Rlim=-1.0, /* Rain rate value where rain changes from frontal to convective. See calc_convlim. */
       Rfactor=1.0, /* factor to multiply R when adjusting to rain gauge */
       acc_threshold, /* threshold of accumulation of precipitation in mm when probabilistic forecasts */
       isotherm=-300.0, /* isotherm temperature of THEIGHT quantity [C] */
       isobar=-1.0;   /* isobar pressure of PHEIGHT quantity [hPa]    */

int    indim[2]={0},
       outdim[2]={0},
       sweepnumber=1,
       datasetnumber=1,
       datanumber=1,
       IRISDATA=0,
       CARTESIAN_HDF5=0,
       POLAR_HDF5=0,
       PGM=0,
       PPM=0,
       PNM=0,
       HDF5=0,
       HOST_LE=1,
       OUTCENTER=0,
       USE_INSW=0,
       USE_INNE=0,
       USE_INCE=0,
       IGN_INCE=0,
       OUT_BBOX_CART=0,
       VERB=0,
       QUIET=0,
       PNMOUT=0,
       QDHEADERS=0,
       GTIFF=0,
       tiling=1,
       HDFOUT=0,
       OVERWRITE=0,
       Zlevel=0,
       RATE_TO_DBZ=0,
       DBZ_TO_RATE=0,
       EPSG=0,
       ODIM_quantcode=-1,
       IRISdata=0,
       Gnodata=255,
       Gundetect=0,
       Acchours=0,
       Forecast=0,
       RetVal=0,
       FLIP_OUTPUT_DATA=0,
       FILL_NODATA=0,
       INSERT_METAFILE=0,
       DEFAULTMETA=1,
       OUTBYTES=0,
       accmins,
       overviews=5,
       def_OverviewList[5]={2,4,8,16,32},
       *eff_OverviewList,
       *OverviewList=NULL,
       ODIM_NODES=0,
       IRISRAIN=0;


long   Chans=1,
       ODIMACCnum=0;

char   projdef[2][1000]={{0}},
       DMS[20]={0},
       cfgfile[1000]={0},
       LUTdir[300]=".",
       sitecode[4]={0},
       infile[1000]={0},
       outformat[100]={0},
       outfile[1000]={0},
       pnmoutfile[1000]={0},
       GTiffoutfile[1000]={0},
       HDFoutfile[1000]={0},
       ODIMquantity[30]={0},
       scaninfofile[300]={0},
       TIFFcompression[200]="DEFLATE",
       sZlevel[3]="6",
       Xtile[6]="256",
       Ytile[6]="256",
       overview_method[50]="CUBIC",
       obstime[20]={0},
       starttimes[50][20]={{0}},
       endtimes[50][20]={{0}},
       QDparam[50]={0},
       ODIMobject[50]="IMAGE",
       ODIMproduct[50]={0},
       ODIMsource[200]="ORG:EFKL",
       fortime_start[13],
       fortime_end[13],
       nodes[5000]={0}, /* list of radar nodes in ODIM HDF5 format */
       metafile[300]; /* contents of this file will be inserted to GeoTIFF GDALMetadata tag */

uint16_t mBZtoR[15000];

int64_t raybins,rays;
uint16_t IRIS_prodtype=0;

char radars[RADN][4]={"KOR","VAN","ANJ","IKA","KUO","VIM","UTA","LUO","KES","PET"};
char coord[3]="XY";

int configreader(char *cfgfile);
int option_handler(int argc, char **argv);
char *toDMS(double degs);
void give_IRIS_geocorners(projUV *IRIScorners, double reflon, double reflat,
     double a, double f, char *projtype, projUV geoC, double *IRISrange);
double geodist(projUV geo1, projUV geo2,double a,double f);
int option_solver(struct option option, char *optarg);
void print_options(void);
double dBZtoR(double dBZ);
double get_csv_item(char *csv, int it);
uint16_t IrisRain(double F);
void date_from_sec(char *date,time_t secs);
time_t sec_from_date(char *date);

int num_options;
size_t readsize;
int reti;
char *retp;

int main(int argc, char *argv[])
{
  projUV lonlat, tocoord, fromcoord;
  int i,j,args,LUT=0;
  int32_t irissize=0;
  long N,insize,outsize,oN,cN,X,Y,relY,oX,oY,C,maxoutN;
  long *LUTarr,lutsize;
  PJ *projref[2]={NULL};
  double halfXres,halfYres,inXoff=0,inYoff=0;
  char buf[300]={0},*p,indatatype[100]={0};
  int shades,optres;
  size_t bytes=1;
  uint8_t *datarr=NULL,*inarr=NULL;
  int T,R,Ci,Co;
  int SW_CORNER = 0,NE_CORNER = 0;
    unsigned char irishdr[640];
  struct stat filestat;
  char inpnmcode[3]={0},outpnmcode[3]="P5";
  char lutname[500]={0},inlutproj[300]={0},outlutproj[300]={0};
  hid_t H5_INF;
  FILE *INF=NULL, *OUTF=NULL, *LUTF=NULL;
 
  setbuf(stdout,NULL);
  setbuf(stderr,NULL);

  /* Endianess of the host
  { 
       uint16_t h=1,c;

       c=htole16(h);
       if(c==h) HOST_LE=1; else HOST_LE=0;
  } 
  */
  eff_OverviewList=def_OverviewList;

  /*================================== INITIALIZATION ====================================================================*/

   /* Set options */
   {
     int oi;

     for(oi=0;;oi++)
     {
       if(optionlist[oi].has_arg<0) { long_options[oi]=null_option; break; }
       long_options[oi].name=optionlist[oi].name;
       long_options[oi].has_arg=optionlist[oi].has_arg;
       long_options[oi].val='A';
       long_options[oi].flag=NULL;
     }
     num_options=oi;
   }

   /* Geographic coords in radians < PI, so initializing them over PI to mark them unused */

   for(T=IN;T<=OUT;T++) for(R=SW;R<=SE;R++) coords[GEO][T][R].u=coords[GEO][T][R].v = 5;

  /*====================== READING  ENVIRONMENT, CONFIGURATION OR OPTIONS  ===========================================*/

   for(i=1;i<argc;i++) 
   { 
     if(!strcmp(argv[i],"--verbose")) { VERB=1;  break; }
     if(!strcmp(argv[i],"--quiet"))   { QUIET=1; break; }
   }

   /* check if any environment variables according to long options list is set */
   if(VERB) printf("OPTIONS");
   for(i=0;;i++)
   {
     char *envstr;
     struct option option;
    
     option=long_options[i];
     if(!option.name) break;
     envstr = getenv(option.name);
     if(envstr) { if(!QUIET) printf("\nEnviron "); option_solver(option,envstr); }
   }
   /*    if(VERB) printf("\n\n"); */


  /* check command line options */
  optres = option_handler(argc,argv);
  if(optres==1) return(0);

  if(VERB) printf("\nFILES\n\n");
  args=argc-optind;

  /* only arguments accepted are input and output files in this order */
  if(args)
  {
    if(args>=2)
    {
      int a;

      sprintf(infile,"%s",argv[optind]); 
      for(a=0;a<args-1;a++)
      {
         sprintf(outfile,"%s",argv[optind+1+a]);
         if(!QUIET) printf("Output file #%d = %s\n",a+1,outfile);
         p=strrchr(outfile,'.');
         if(p) 
	 {
	   if(strstr(p+1,"pgm") || strstr(p+1,"ppm")) { PNMOUT=1; sprintf(pnmoutfile,"%s",outfile); }
	   if(strstr(p+1,"tif")) { GTIFF=1;  sprintf(GTiffoutfile,"%s",outfile); }
	   if(strstr(p+1,"hdf") || strstr(p+1,"h5")) { HDFOUT=1; sprintf(HDFoutfile,"%s",outfile); }
	 }
         if(! (p || GTIFF || PNMOUT)) { if(!QUIET) printf("Unknown output format, use extensions pgm, ppm, tif, hdf or h5!"); return(1); }
      }
    }
  

    while(args==1)
    {
      if(infile[0])  { sprintf(outfile,"%s",argv[optind]); break; }
      if(outfile[0]) { sprintf(infile,"%s",argv[optind]);  break; }
      break;
    }
  }

  if(!infile[0] || !outfile[0]) 
  { 
     if(args && !QUIET) printf("Give input and at least one output file!\n\n"); 
     else print_options();  
     return(1);
  }

  if(VERB) printf("Input file %s\nOutput file %s\n\nINPUT GEOMETRY\n==============\n\n",infile,outfile);

  ZRB10=ZRB10c = 10.0*ZRB;
  ZRC=ZRCc = -log10(ZRA)/ZRB;

  if(ZRAc>0.0 && ZRBc>0.0)
  { 
     ZRB10c = 10.0*ZRBc;
     ZRCc = -log10(ZRAc)/ZRBc;

     /* Calculating the dBZ value where frontal and convective rain rate are the same. 
        Beyond this value the convective ZRAc and ZRBc are used in dBZ -> RR conversion */
 
     if(ZRA==ZRAc) 
     { 
        conv_dBZlim = 10.0*log10(ZRA);
        conv_Rlim=dBZtoR(conv_dBZlim);
     }
     else 
     {
        double R;

        R = pow(ZRA/ZRAc,1.0/(ZRBc-ZRB));
        conv_dBZlim = 10.0 * log10(ZRA*pow(R,ZRB));
        conv_Rlim=R;
     }
  }


  /*================================== INPUT DATA FORMAT TEST ====================================================================*/

  /* Checking the input format and open the input file */

  {
     stat(infile,&filestat);

     INF=NULL;
     H5Eset_auto(H5E_DEFAULT,NULL,NULL);
     H5_INF = H5Fopen(infile, H5F_ACC_RDONLY, H5P_DEFAULT);
     if(H5_INF >= 0) { HDF5=1; } else
     {  
       HDF5=0;
       INF=fopen(infile,"r");
       if(!INF && !QUIET) { printf("Unable to open file %s\n",infile); return(1); }
       readsize=fread(irishdr,640,1,INF);
       while(1)
       {
          memcpy(&irissize,irishdr+4,4);
          if(irishdr[0]==27 && irishdr[12]==26 && irissize==(int32_t)filestat.st_size)
          {
            IRISDATA=1;
            break;
          }
          memcpy(inpnmcode,irishdr,2);
          if(strcmp(inpnmcode,"P5")==0) {PGM=1; rewind(INF); break; }
          if(strcmp(inpnmcode,"P6")==0) {PPM=1; rewind(INF); break; }
    
          if(!QUIET) printf("Unknown input data format!\n");
          return(1);
  
       } 
     }
  }
  

  /*================================== PROCESSING OF HDF5 INPUT  ==========================================*/

  if(HDF5) 
  {
      char objectStr[100];
      char sourceStr[200]={0};
      hid_t datatype=0;
      herr_t  ret;   
      hsize_t dims[2];
      size_t typesize;
      double rgain,roffset,fingain=0.5,finoffset=-32.0;

      ret=H5LTget_attribute_string(H5_INF, "what", "object", objectStr);
      if(strcmp(objectStr,"COMP") == 0 || strcmp(objectStr,"IMAGE") == 0)
      {
	  H5T_class_t dataclass;
          int64_t xsize,ysize;
          unsigned char *floatarr;
          long values;
          double LL_lon,LL_lat,UR_lon,UR_lat;
          char tmp_projdef[300]={0},real_projdef[300]={0},datawhat[50],setdata[50],quantity[50]={0};
          unsigned char undetect,nodata;
          double fnodata,fundetect;

          CARTESIAN_HDF5=1;
          
          sprintf(ODIMobject,"%s",objectStr);
          ret=H5LTget_attribute_string(H5_INF, "what", "source", sourceStr);
#if __SIZEOF_LONG__==4
          ret=H5LTget_attribute_long_long(H5_INF,"where","xsize",&xsize);          
          ret=H5LTget_attribute_long_long(H5_INF,"where","ysize",&ysize);          
#else
          ret=H5LTget_attribute_long(H5_INF,"where","xsize",&xsize);          
          ret=H5LTget_attribute_long(H5_INF,"where","ysize",&ysize);          
#endif
          ret=H5LTget_attribute_string(H5_INF,"where","projdef",tmp_projdef);
          ret=H5LTget_attribute_double(H5_INF,"where","LL_lon",&LL_lon);          
          ret=H5LTget_attribute_double(H5_INF,"where","LL_lat",&LL_lat);          
          ret=H5LTget_attribute_double(H5_INF,"where","UR_lon",&UR_lon);          
          ret=H5LTget_attribute_double(H5_INF,"where","UR_lat",&UR_lat);          
          ret=H5LTget_attribute_double(H5_INF,"where","xscale",&inres[XD]);          
          ret=H5LTget_attribute_double(H5_INF,"where","yscale",&inres[YD]);          
          ret=H5LTget_attribute_string(H5_INF,"how","nodes",nodes);
          if(ret>=0) ODIM_NODES = 1;          

          indim[XD]=(int)xsize;
          indim[YD]=(int)ysize;
          coords[GEO][IN][SW].u=DEG_TO_RAD*LL_lon;
          coords[GEO][IN][SW].v=DEG_TO_RAD*LL_lat;
          coords[GEO][IN][NE].u=DEG_TO_RAD*UR_lon;
          coords[GEO][IN][NE].v=DEG_TO_RAD*UR_lat;
          if(strstr(tmp_projdef,"proj=lonlat") || strstr(tmp_projdef,"proj=latlon") ||
             strstr(tmp_projdef,"proj=longlat") || strstr(tmp_projdef,"proj=latlong"))
	  {
	    /*
	     inres[XD] *= DEG_TO_RAD;
	     inres[YD] *= DEG_TO_RAD;
	    */
	    inres[XD] = (coords[GEO][IN][NE].u - coords[GEO][IN][SW].u)/(double)indim[XD];
	    inres[YD] = (coords[GEO][IN][NE].v - coords[GEO][IN][SW].v)/(double)indim[YD];
	  }
          values=xsize*ysize;
          for(Co=XD;Co<=YD;Co++) inrange[Co]=indim[Co]*inres[Co];

	  {
	    char *tp,*p,*s,param[100];
            tp=tmp_projdef;
            while(p=strchr(tp,'+'))
	    { 
	       memset(param,0,100);
               s=strchr(p,' ');
               if(!s) s=strchr(p,0);
               strncpy(param,p,s-p+1);
               tp=s;
               if(strstr(param,"x_0=")) continue;
               if(strstr(param,"y_0=")) continue;
               strcat(real_projdef,param);
               /* sprintf(real_projdef,"%s%s",real_projdef,param); */
	    }
	  }
          
          sprintf(datawhat,"/dataset%d/data%d/what",datasetnumber,datanumber);
          sprintf(setdata,"/dataset%d/data%d/data",datasetnumber,datanumber);

          ret=H5LTget_attribute_string(H5_INF,datawhat,"quantity",quantity);
          if(ret<0)
	  {
             sprintf(datawhat,"/dataset%d/what",datasetnumber);
             ret=H5LTget_attribute_string(H5_INF,datawhat,"quantity",quantity);
	  }
          ret=H5LTget_attribute_double(H5_INF,datawhat,"gain",&inGain);          
          ret=H5LTget_attribute_double(H5_INF,datawhat,"offset",&inOffset);          
          ret=H5LTget_attribute_double(H5_INF,datawhat,"nodata",&inNodata);          
          ret=H5LTget_attribute_double(H5_INF,datawhat,"undetect",&inUndetect);          
          sprintf(indatatype,"Cartesian HDF5, quantity %s\n",quantity);
          sprintf(ODIMquantity,"%s",quantity);

          ret = H5LTget_dataset_info(H5_INF,setdata,dims,&dataclass,&typesize);

          if(dataclass == H5T_INTEGER)
	  { 
   	      inarr=calloc(values,typesize);
	      Ggain=fingain;
              Goffset=finoffset;
              Gnodata=255;
              Gundetect=0;

	      if(typesize==1)
	      {
		 if(!QUIET) printf("Reading %d byte integer data as unsigned char\n",typesize);
		 ret = H5LTread_dataset(H5_INF,setdata,H5T_NATIVE_UCHAR,inarr); 
    	      }

	      if(typesize==2)
	      { 
		 if(!QUIET) printf("Reading %d byte integer data as unsigned short\n",typesize);
		 ret = H5LTread_dataset(H5_INF,setdata,H5T_NATIVE_USHORT,inarr); 
    	      }
              bytes=typesize;

              while(1)
	      {
	         if(strcmp(quantity,"DBZH")==0 || strcmp(quantity,"TH")==0)
	         {
		    ODIM_quantcode=ODIM_DBZH;
		    if(inGain!=fingain || inOffset!=finoffset) 
		    { 
		       if(bytes == 1)
		       {
		          unsigned char N,M;
		          long i;             
		          double gain,offset,fM;

		          undetect=(unsigned char)inUndetect;          
		          nodata=(unsigned char)inNodata;
		          gain=inGain/fingain;
		          offset=(inOffset-finoffset)/fingain;
		          for(i=0;i<values;i++)
		          {
			     N=inarr[i];
			     while(1)
			     {
			        if(N==undetect) { M=0; break; }
			        if(N==nodata) { M=255; break; }
			        fM=gain*(double)N+offset;
			        if(fM<0) { M=0; break; }
			        if(fM>=255.0) { M=254; break; }
			        M=(unsigned char)fM;
			        break;
			     }
			     inarr[i]=M;
		          }
		       }
		    }
                    break;
		 }

	         if(strcmp(quantity,"RATE")==0)
		 {
		    char date[9],time[7];
                    time_t startsecs,endsecs;

		    Ggain=inGain;
                    Goffset=inOffset;
                    Gundetect=inUndetect;
                    Gnodata=inNodata;
                    ODIM_quantcode=ODIM_RATE;

                    ret=H5LTget_attribute_string(H5_INF,"what","date",date);
                    ret=H5LTget_attribute_string(H5_INF,"what","time",time);
                    sprintf(obstime,"%s%.4s",date,time);

                    sprintf(datawhat,"/dataset%d/what",datasetnumber);

                    ret=H5LTget_attribute_string(H5_INF,datawhat,"startdate",date);
                    ret=H5LTget_attribute_string(H5_INF,datawhat,"starttime",time);
                    sprintf(fortime_start,"%s%.4s",date,time);
                    startsecs=sec_from_date(fortime_start);

                    ret=H5LTget_attribute_string(H5_INF,datawhat,"enddate",date);
                    ret=H5LTget_attribute_string(H5_INF,datawhat,"endtime",time);
                    sprintf(fortime_end,"%s%.4s",date,time);
                    endsecs=sec_from_date(fortime_end);

                    accmins=(endsecs-startsecs)/60;

                    sprintf(datawhat,"/dataset%d/data%d/what",datasetnumber,datanumber);
                    break;
		 }

	         if(strcmp(quantity,"PROB")==0)
		 {
		    char date[9],time[7];
                    time_t startsecs,endsecs;

		    Ggain=inGain;
                    Goffset=inOffset;
                    Gundetect=inUndetect;
                    Gnodata=inNodata;
                    ODIM_quantcode=ODIM_ACCPROB;

                    ret=H5LTget_attribute_string(H5_INF,"what","date",date);
                    ret=H5LTget_attribute_string(H5_INF,"what","time",time);
                    sprintf(obstime,"%s%.4s",date,time);

                    sprintf(datawhat,"/dataset%d/what",datasetnumber);

                    ret=H5LTget_attribute_string(H5_INF,datawhat,"startdate",date);
                    ret=H5LTget_attribute_string(H5_INF,datawhat,"starttime",time);
                    sprintf(fortime_start,"%s%.4s",date,time);
                    startsecs=sec_from_date(fortime_start);

                    ret=H5LTget_attribute_string(H5_INF,datawhat,"enddate",date);
                    ret=H5LTget_attribute_string(H5_INF,datawhat,"endtime",time);
                    sprintf(fortime_end,"%s%.4s",date,time);
                    endsecs=sec_from_date(fortime_end);

                    accmins=(endsecs-startsecs)/60;

                    sprintf(datawhat,"/dataset%d/data%d/what",datasetnumber,datanumber);
		    ret=H5LTget_attribute_double(H5_INF,datawhat,"threshold_value",&acc_threshold);
                    break;
		 }

	         if(strcmp(quantity,"HCLASS")==0)
		 {
		    char date[9],time[7];
                    time_t startsecs,endsecs;

		    Ggain=inGain;
                    Goffset=inOffset;
                    Gundetect=inUndetect;
                    Gnodata=inNodata;
                    ODIM_quantcode=ODIM_HCLASS;

                    ret=H5LTget_attribute_string(H5_INF,"what","date",date);
                    ret=H5LTget_attribute_string(H5_INF,"what","time",time);
                    sprintf(obstime,"%s%.4s",date,time);

                    sprintf(datawhat,"/dataset%d/what",datasetnumber);

                    ret=H5LTget_attribute_string(H5_INF,datawhat,"startdate",date);
                    ret=H5LTget_attribute_string(H5_INF,datawhat,"starttime",time);
                    sprintf(fortime_start,"%s%.4s",date,time);
                    startsecs=sec_from_date(fortime_start);

                    ret=H5LTget_attribute_string(H5_INF,datawhat,"enddate",date);
                    ret=H5LTget_attribute_string(H5_INF,datawhat,"endtime",time);

                    break;
		 }

	         if(!QUIET) printf("Cartesian HDF5 quantity is not DBZH, TH, RATE, PROB or HCLASS, exiting...\n"); 
                 H5Fclose(H5_INF); 
                 return(1); 
             }
	  }
      
          if(dataclass == H5T_FLOAT)
	  { 
            unsigned char N;
            uint16_t N2;
            float F;
            double D,la,dBZ,fN,RR;
            int RATE=0,ACC=0,DBZ=0;

            la=10.0*log10(ZRA);
            bytes=2;
            if(strcmp(quantity,"RATE")==0) RATE=1;
            if(strcmp(quantity,"DBZH")==0) { DBZ=1; RATE=1; RATE_TO_DBZ=1; printf("DBZ DATA READ\n"); }
            if(strcmp(quantity,"ACCR")==0 || strcmp(quantity,"ACRR")==0)
            { 
               ACC=1; 
               ODIM_quantcode = ODIM_ACRR; 
               sprintf(ODIMquantity,"ACRR");
            }

            if(RATE && RATE_TO_DBZ) bytes=1;

    	    floatarr=calloc(values,typesize);
	    if(typesize == 8)
	    {
    		 if(!QUIET) printf("Reading %d byte float data as double\n",typesize);
		 ret = H5LTread_dataset(H5_INF,setdata,H5T_NATIVE_DOUBLE,floatarr);
    	    }
    
	    if(typesize == 4)
	    {
    		 if(!QUIET) printf("Reading %d byte float data as float\n",typesize);
		 ret = H5LTread_dataset(H5_INF,setdata,H5T_NATIVE_FLOAT,floatarr);
	    }
           
            inarr=malloc(values*bytes);
            if(RATE)
	    { 
	        if(RATE_TO_DBZ)
	        {
		    sprintf(ODIMquantity,"DBZH");
                    ODIM_quantcode = ODIM_DBZH;
		} else 
                { 
		   sprintf(ODIMquantity,"RATE"); 
                   ODIM_quantcode = ODIM_RATE; 
                   Ggain=0.01;
                   Goffset=0.0;
                   Gnodata=65535;
                   Gundetect=0;
                }
	    }

            if(RATE_TO_DBZ) { Ggain=fingain; Goffset=finoffset; }

            for(i=0;i<values;i++)
	    {
	       if(typesize==8) memcpy(&D,&floatarr[8*i],8);
               else
	       {
                   memcpy(&F,&floatarr[4*i],4);
                   D=(double)F;
	       }
    
               if(RATE || ACC)
	       {
		 if(RATE_TO_DBZ)
		 {
                   bytes=1;
                   while(1)
		   {
		     if(D==inUndetect) { N=0; break; }
		     if(D==inNodata)   { N=255; break; }
                     RR=inGain*D+inOffset;
                     if(DBZ) dBZ=RR;
                     else
		     {
		        if(RR<0.0)   { N=0; break; }
                        dBZ = la + ZRB10*log10(RR);
		     }
                     fN=(dBZ-finoffset)/fingain;
                     if(fN<0.0) { N=0; break; }
                     if(fN>=254.5) { N=254; break; }
                     N=(unsigned char)(fN+0.5);
                     break;
		   }
                   inarr[i]=N;
		 } else
		 {
                   while(1)
		   {
		     if(D==inUndetect) { N2=0; break; }
		     if(D==inNodata)   { N2=65535; break; }
                     RR=inGain*D+inOffset;
		     if(RR<0.0)   { N2=0; break; }
                     if(IRISRAIN) N2=IrisRain(RR); 
                     else 
		     {
		        if(RR>655.33)   { N2=65534; break; }
                        N2=(uint16_t)(RR*100.0);
		     }
                     break;
		   }
                   /* inarr[i]=N2; */
                   memcpy(&inarr[2*i],&N2,2);
		 }
	       }
	    }
            free(floatarr);
	    /*
            H5Fclose(H5_INF);
            return(1);
	    */
	  }    

    	  H5Fclose(H5_INF);
          sprintf(lutname,"%s/HDF5Cart_%s_%dx%d_%.0f,%.0f",LUTdir,sourceStr,indim[XD],indim[YD],inres[XD],inres[YD]);
          sprintf(projdef[IN],"%s",real_projdef);
	  projref[IN]=pj_init_plus(projdef[IN]);
          for(Ci=SW;Ci<=NE;Ci++) coords[CART][IN][Ci] = pj_fwd(coords[GEO][IN][Ci],projref[IN]);
      }  

      if(strcmp(objectStr,"PVOL") == 0 || strcmp(objectStr,"SCAN") == 0)  /* Polar HDF5 data */	    
      {
	  int quantI,LUT=0;
          char seekset[200];
          char quantity[30]={0};
          char dataset[200],datawhere[200],datahow[200],wanted_dataset[200];
          double cartres,cartradius,halfres,R,P,Rbin,tmplon,tmplat;
          double azres,halfbeam,zcartres,cartsub,U,V,fnodata,fundetect;
          long cartdim,XI,YI,N,AzI,binI,pN,Adim,Rdim,bins,maxbins,OY,arrsize;
          /* long insize; */
          unsigned char *bytearr;

          POLAR_HDF5=1;
    	  sprintf(indatatype,"Polar HDF5, quantity %s\n",ODIMquantity);
    	  if(VERB) printf("Polar HDF5, quantity %s\n",ODIMquantity);
	  ret=H5LTget_attribute_string(H5_INF, "what", "source", sourceStr);

	  if(ret>=0)
	  {
	    char *nod;
	    int i;

	    nod=strstr(sourceStr,"NOD:");
            if(nod)
	    {
               nod+=6;
	       for(i=0;i<3;i++) sitecode[i] = *(nod+i)-32;
	    }
	  }

	  for(quantI=1;;quantI++)
	  {
	       sprintf(wanted_dataset,"/dataset%d/data%d",sweepnumber,quantI);
	       sprintf(seekset,"%s/what",wanted_dataset);
	       ret=H5LTget_attribute_string(H5_INF, seekset, "quantity", quantity);
	       /*  printf("ret:%d %s\n",ret,ODIMquantity); */
	       if(!strcmp(ODIMquantity,quantity)) break;
	       if(ret<0) break;
	  }
	  if(ret<0) { if(!QUIET) printf("Quantity %s of sweep %d not found!\n",ODIMquantity,sweepnumber); return(1); }
          ret=H5LTget_attribute_double(H5_INF, seekset, "gain", &Ggain);
          ret=H5LTget_attribute_double(H5_INF, seekset, "offset", &Goffset);
          ret=H5LTget_attribute_double(H5_INF, seekset, "nodata", &fnodata);
          ret=H5LTget_attribute_double(H5_INF, seekset, "undetect", &fundetect);
          Gnodata=(int)fnodata;
          Gundetect=(int)fundetect;

	  ret=H5LTget_attribute_double(H5_INF, "where", "lon", &tmplon);   
	  ret=H5LTget_attribute_double(H5_INF, "where", "lat", &tmplat);
	  if(VERB) printf("HDF5 center %.6f %.6f\n",tmplon,tmplat);
	  IRISgeo.u=tmplon*DEG_TO_RAD;   
	  IRISgeo.v=tmplat*DEG_TO_RAD;   

	  sprintf(projdef[IN],"+proj=aeqd +ellps=WGS84 +lon_0=%.7f +lat_0=%.7f",
		  IRISgeo.u/DEG_TO_RAD,IRISgeo.v/DEG_TO_RAD);

	  projref[IN]=pj_init_plus(projdef[IN]);


	  sprintf(datawhere,"/dataset%d/where",sweepnumber);
	  sprintf(datahow,"/dataset%d/how",sweepnumber);
	  H5LTget_attribute_double(H5_INF, datawhere, "elangle", &elangle);
	  elangle*=DEG_TO_RAD;
	  H5LTget_attribute_double(H5_INF, datawhere, "rscale", &binlen);
#if __SIZEOF_LONG__==4
	  H5LTget_attribute_long_long(H5_INF, datawhere, "nbins", &raybins); 
	  H5LTget_attribute_long_long(H5_INF, datawhere, "nrays", &rays);
#else
	  H5LTget_attribute_long(H5_INF, datawhere, "nbins", &raybins); 
	  H5LTget_attribute_long(H5_INF, datawhere, "nrays", &rays);
#endif
	  azres=(double)rays/360.0;
          halfbeam=0.5/azres;


	  if(VERB) printf("Elangle %.1f, bins %d, binlen %.0f, rays %d\n",elangle/DEG_TO_RAD,(int)raybins,binlen,(int)rays);
          if(polarrange < 0) maxbins=raybins; else maxbins=polarrange*1000.0/binlen; 
	  cartres = binlen*cos(elangle);
	  cartradius = (double)maxbins*cartres;
	  inrange[XD]=inrange[YD] = cartradius*2;
	  cartdim=indim[XD]=indim[YD] = maxbins*2*polarzoom;
	  halfres=cartres/2;
	  sprintf(lutname,"%s/PolCar_Dim=%ld_Elev=%.1f_Bins=%d_Binlen=%.0f_Rays=%d_Zoom=%.0f.lut",
		  LUTdir,cartdim,elangle/DEG_TO_RAD,(int)maxbins,binlen,(int)rays,polarzoom);

	  inres[XD]=inres[YD] = cartres/polarzoom;   
	  coords[GEO][IN][CE] = IRISgeo;

	  coords[CART][IN][CE]=pj_fwd(coords[GEO][IN][CE],projref[IN]);

	  coords[CART][IN][SW].u=coords[CART][IN][CE].u - cartradius;
	  coords[CART][IN][SW].v=coords[CART][IN][CE].v - cartradius;
	  coords[CART][IN][NE].u=coords[CART][IN][CE].u + cartradius;
	  coords[CART][IN][NE].v=coords[CART][IN][CE].v + cartradius;

	  coords[GEO][IN][SW]=pj_inv(coords[CART][IN][SW],projref[IN]);
	  coords[GEO][IN][NE]=pj_inv(coords[CART][IN][NE],projref[IN]);

	  sprintf(dataset,"%s/data",wanted_dataset); 
	  ret = H5LTget_dataset_info(H5_INF,dataset,dims,NULL,&bytes);
	  if(bytes==1) datatype=H5T_NATIVE_UCHAR;
	  if(bytes==2) datatype=H5T_NATIVE_USHORT; 
	  Adim=(long int)dims[0];
	  Rdim=(long int)dims[1];
	  bins=Adim*Rdim;
	  bytearr=malloc(bins*bytes);
	  ret = H5LTread_dataset(H5_INF,dataset,datatype,bytearr);
	  shades=pow(2,8*bytes)-1;

          if(ODIM_quantcode==ODIM_VRAD && bytes==1) 
	  {
             H5LTget_attribute_double(H5_INF, datahow, "NI", &NI);
	     /*             Ggain *= NI;
			    Goffset *= NI; */
	  }
          
	  /*
	  if(bytes==2)
	  { 
	    if(!CONVERT_TO_8)
	    {
	      shades=65535;
	    } else
	      {
		 H5LTget_attribute_double(H5_INF, wanted_dataset, "gain", &gain);
		 H5LTget_attribute_double(H5_INF, wanted_dataset, "offset", &offset);
		 H5LTget_attribute_double(H5_INF, wanted_dataset, "nodata", &nodata);
		 H5LTget_attribute_double(H5_INF, wanted_dataset, "undetect", &undetect);
	      }
	  } else shades=255;
	  */
	  arrsize=cartdim*cartdim;
	  lutsize=arrsize*sizeof(long);
	  insize=arrsize*bytes;
	  inarr=malloc(insize);
	  memset(inarr,255,insize);

	  zcartres=cartres/polarzoom;
	  cartsub=-cartradius; /* +halfres; */ 


	  /* READING / GENERATING OF POLAR TO CARTESIAN LUT */

	  LUTarr=malloc(lutsize);
	  memset(LUTarr,255,lutsize);    
	  LUTF=fopen(lutname,"r");
	  if(LUTF) 
	  {
	     if(VERB) printf("Reading polar to cartesian LUT %s\n",lutname);
	     readsize=fread(LUTarr,lutsize,1,LUTF);
	     LUT = 1;
	  }   
	  else 
	  {  
	     LUTF=fopen(lutname,"w");
	     if(VERB) printf("Generating polar to cartesian LUT %s\n",lutname);
	     LUT = 0;
	  }
    
	  for(YI=0;YI<cartdim;YI++)
	  { 
	    /* OY=(cartdim-YI-1)*cartdim; */
	    OY=YI*cartdim;
	    V = (double)YI*zcartres+cartsub; 
	    for(XI=0;XI<cartdim;XI++)
	    {
	       N=OY+XI;
	       if(!LUT)
	       {
		  U = (double)XI*zcartres+cartsub; 
		  R=sqrt(V*V + U*U);
		  Rbin=R/cartres;
		  binI=(int)Rbin;
		  if(binI>=raybins) continue;
		  P=atan2(U,V)/DEG_TO_RAD;
		  if(P<0) P+=360;
		  AzI=(int)((P+halfbeam)*azres);
		  if(AzI >= rays) AzI-=rays;
		  pN=AzI*raybins+binI;
		  LUTarr[N]=pN;
	       } else pN=LUTarr[N];
	       if(pN<0) continue;

	       memcpy(&inarr[N*bytes],&bytearr[pN*bytes],bytes); 
	    }
	  }
         
          if(bytes==1 && (strcmp(quantity,"DBZH")==0 || strcmp(quantity,"TH")==0 || strcmp(quantity,"DBZH_QC")==0))
	  {
		    if(Ggain!=fingain || Goffset!=finoffset) 
		    { 
		       unsigned char N,M,undetect,nodata;
		       long i;             
		       double gain,offset,fM;

		       undetect=(unsigned char)fundetect;          
		       nodata=(unsigned char)fnodata;
		       gain=Ggain/fingain;
		       offset=(Goffset-finoffset)/fingain;
		       for(i=0;i<insize;i++)
		       {
			  N=inarr[i];
			  while(1)
			  {
			    if(N==undetect) { M=0; break; }
			    if(N==nodata) 
                            { 
			       if(FILL_NODATA) M=0; else M=255; 
                               break; 
                            }
			    fM=gain*(double)N+offset;
			    if(fM<0) { M=0; break; }
			    if(fM>=255.0) { M=254; break; }
			    M=(unsigned char)fM;
			    break;
			  }
			  inarr[i]=M;
		       }
                       Ggain=fingain;
                       Goffset=finoffset;
                       Gnodata=255;
                       Gundetect=0;
		    }                 
	  }
      

	  /*	  if(bytes==2 && HOST_LE) swab(inarr,inarr,insize); */
          inGain=Ggain;
          inOffset=Goffset;
          inUndetect=Gundetect;
          inNodata=Gnodata;

	  H5Fclose(H5_INF);
	  free(bytearr);
	  if(!IGN_INCE) USE_INCE=1;
	  if(!LUT) fwrite(LUTarr,lutsize,1,LUTF);
	  fclose(LUTF);
	  free(LUTarr);
          lutname[0]=0;
          if(VERB) printf("Polar HDF5 data %s read\n",infile);    
      }
      if(!(POLAR_HDF5 | CARTESIAN_HDF5)) { if(!QUIET) printf("HDF5 input data is not polar nor cartesian, exiting...\n"); return(1); }
  }

  /*================================== INPUT GEOMETRY ====================================================================*/

  /* Set the input geometrical properties for conversion */
  if(PGM || PPM)
  {
    int IS_CORNER=0;
    int STERE=0;
    char lats[20]={0},lons[20]={0},tslats[20]={0},*p=NULL,*qdp=NULL;

    qdp=strstr(QDparam,"Precipitation");
    for(Co=SW;Co<=NE;Co++) if(coords[GEO][IN][Co].u < 4) IS_CORNER++; 
    sprintf(outpnmcode,"%s",inpnmcode);
    for(i=0;i<3;i++)
    {
         retp=fgets(buf,299,INF);
         p=NULL;
         if(p=strstr(buf,"bottomleft"))
         { 
            sscanf(p,"%*s %lf %lf",&coords[GEO][IN][SW].u,&coords[GEO][IN][SW].v);
            coords[GEO][IN][SW].u *= DEG_TO_RAD;
            coords[GEO][IN][SW].v *= DEG_TO_RAD;
            IS_CORNER=1;
         }
             
         if(p=strstr(buf,"topright")) 
         {
            sscanf(p,"%*s %lf %lf",&coords[GEO][IN][NE].u,&coords[GEO][IN][NE].v);
            coords[GEO][IN][NE].u *= DEG_TO_RAD;
            coords[GEO][IN][NE].v *= DEG_TO_RAD;
            IS_CORNER=1;
         }

         p=strstr(buf,"Precipitation");
         if(p || qdp)
         { 
	     char *q;            

             if(qdp) p=qdp;
             p+=13;
             if(q=strstr(p,"Rate")) 
	     sprintf(ODIMquantity,"RATE"); 
             ODIM_quantcode = ODIM_RATE; 
             
             if(q=strchr(p,'h'))
	     {
	       char shours[3]="\0\0\0";
               
               strncpy(shours,p,q-p);
               Acchours=atoi(shours);
               sprintf(ODIMquantity,"ACRR"); 
               ODIM_quantcode = ODIM_ACRR; 
	     }
             Ggain=0.01;
             Goffset=0.0;
             Gnodata=65535;
             Gundetect=0;
         }
         if((p=strstr(buf,"CorrectedReflectivity")) || ODIM_quantcode==ODIM_DBZH || ODIM_quantcode==ODIM_DBZH_QC) 
         { 
             sprintf(ODIMquantity,"DBZH"); 
             ODIM_quantcode = ODIM_DBZH; 
             Ggain=0.5;
             Goffset=-32.0;
             Gnodata=255;
             Gundetect=0;
         }
         if(p=strstr(buf,"BackgroundPrecipitationType")) 
         { 
             sprintf(ODIMquantity,"HCLASSBG"); 
             ODIM_quantcode = ODIM_HCLASSBG; 
             Ggain=1.0;
             Goffset=0.0;
             Gnodata=255;
             Gundetect=0;
         }
         if(p=strstr(buf,"obstime")) sprintf(obstime,"%.12s",strchr(p,' ')+1); 
         if(p=strstr(buf,"fortime")) { Forecast=1; sprintf(fortime_start,"%.12s",strchr(p,' ')+1); } 
	 
         if(p=strstr(buf,"stereographic")) STERE=1;

         if(STERE)
	 {
	    if(p=strstr(buf,"centrallongitude")) sprintf(lons,"%s",strchr(p,' ')+1);
            if(p=strstr(buf,"centrallatitude")) sprintf(lats,"%s",strchr(p,' ')+1);
            if(p=strstr(buf,"truelatitude")) sprintf(tslats,"%s",strchr(p,' ')+1);
	 }

         if(buf[0]=='#') { i--; continue; }

         if(i==1) sscanf(buf,"%d %d\n",&indim[XD],&indim[YD]);
         if(i==2) sscanf(buf,"%d\n",&shades);
    }
    if(shades > 255) bytes=2; else bytes=1;
    if(Acchours) printf("Input PGM %s %d hours\n",ODIMquantity,Acchours);

    if(lons[0] && lats[0] && tslats[0]) 
      if(!projdef[IN][0]) sprintf(projdef[IN],"+proj=stere +a=6371288 +lon_0=%f +lat_0=%f +lat_ts=%f",atof(lons),atof(lats),atof(tslats)); 
    if(!projdef[IN][0]) { if(!QUIET) printf("Give --inproj= PROJ4 projection string for PNM data!\n"); return(1); }
    projref[IN]=pj_init_plus(projdef[IN]);

    if(! IS_CORNER)
    {
       if(coords[IN][GEO][CE].u > 4) { if(!QUIET) printf("Not even center input geo coords!\n"); return(1); }
       coords[CART][IN][CE]=pj_fwd(coords[GEO][IN][CE],projref[IN]);
       for(Co=XD;Co<=YD;Co++)
       { 
         while(1)
	 {
	   if(inres[Co] && indim[Co]) { inrange[Co]=inres[Co]*indim[Co]; break; }
	   if(inrange[Co] && indim[Co]) { inres[Co]=inrange[Co]/indim[Co]; break; }
           break;
	 }
       }
       coords[CART][IN][SW].u=coords[CART][IN][CE].u-inrange[XD]/2.0;
       coords[CART][IN][SW].v=coords[CART][IN][CE].v-inrange[YD]/2.0;
       coords[CART][IN][NE].u=coords[CART][IN][CE].u+inrange[XD]/2.0;
       coords[CART][IN][NE].v=coords[CART][IN][CE].v+inrange[YD]/2.0;
       for(Ci=SW;Ci<=NE;Ci++) coords[GEO][IN][Ci] = pj_inv(coords[CART][IN][Ci],projref[IN]);
    }
    else
    {

       /* Ci used both index for corner and after corner conversion to cartesian index for X and Y */
       for(Ci=SW;Ci<=NE;Ci++) coords[CART][IN][Ci] = pj_fwd(coords[GEO][IN][Ci],projref[IN]);

       inrange[XD]=coords[CART][IN][NE].u-coords[CART][IN][SW].u;
       inrange[YD]=coords[CART][IN][NE].v-coords[CART][IN][SW].v;
       for(Co=XD;Co<=YD;Co++) inres[Co]=inrange[Co]/(double)indim[Co];
    }

    if(PPM) { Chans=3; sprintf(indatatype,"PPM %d bytes",bytes); } else sprintf(indatatype,"PGM %d bytes",bytes);
    PNM=1;
  }

  /* Reading geometry and projection metadata from IRIS input data */  
  if(IRISDATA)
  {

	  int32_t IRISres[2],IRISdim[2],IRPar1,IRPar2;
	  int32_t IRIScentlon,IRIScentlat,Area_reflat,Area_reflon;
	  uint32_t ERad_cm,IRFlatten;
          int16_t TOPb;
	  uint8_t projcode;
	  double std_par1,std_par2;
	  double IRISrange[2];
	  projUV IRIScorners[2];
	  int Base=332,OPTCENTER=0;
	  char IRIS_projdef[300],cornersname[300];
	  FILE *OFFSFILE;  

	  bytes=1;
	  coords[CART][IN][CE].u=0.0;     
	  coords[CART][IN][CE].v=0.0;       
          if(coords[GEO][IN][CE].u < 4) OPTCENTER = 1;

          memcpy(&IRIS_prodtype,irishdr+24,2);
          if(IRIS_prodtype == IRIS_CAPPI)
	  {
	     ODIM_quantcode = ODIM_DBZH;
             Ggain=0.5;
             Goffset=-32.0;
          }   
          if(IRIS_prodtype == IRIS_TOPS) 
          { 
             memcpy(&TOPb,irishdr+180,2); 
             TOPthr=(double)TOPb/16.0; 
	     ODIM_quantcode = ODIM_ETOP;
             Ggain=0.1;
             Goffset=-0.1;
          }
	  memcpy(&IRISres[XD],irishdr+100,4);
	  memcpy(&IRISres[YD],irishdr+104,4);
	  memcpy(&IRISdim[XD],irishdr+112,4);
	  memcpy(&IRISdim[YD],irishdr+116,4);
	  memcpy(&IRIScentlat,irishdr+440,4);
	  memcpy(&IRIScentlon,irishdr+444,4);
	  sprintf(sitecode,"%.3s",irishdr+406);

	  indim[XD]=IRISdim[XD];
	  indim[YD]=IRISdim[YD];
          if(!OPTCENTER)
	  {
	     IRISgeo.u=M_PI*((double)IRIScentlon/((double)INT32_MAX + 1));
	     IRISgeo.v=M_PI*((double)IRIScentlat/((double)INT32_MAX + 1));
	  } else IRISgeo=coords[GEO][IN][CE];

          if(inrange[XD] || inrange[YD])
	  { 
	    if(inrange[XD]) IRISrange[XD]=inrange[XD]/2.0;
            if(inrange[YD]) IRISrange[YD]=inrange[YD]/2.0;
	  }
	  else for(Co=XD;Co<=YD;Co++) IRISrange[Co] = (double)IRISres[Co]*(double)IRISdim[Co]/200.0; /* range from center, IRISres [cm] */

	  if(VERB) printf("IRIS %s center lon %f, lat %f\n",sitecode,IRISgeo.u/DEG_TO_RAD,IRISgeo.v/DEG_TO_RAD);
	  if(VERB) printf("IRIS Xres %.2f Yres %.2f, Xdim %d Ydim %d\n",IRISres[XD]/100.0,IRISres[YD]/100.0,IRISdim[XD],IRISdim[YD]);


	  sprintf(indatatype,"IRIS area");
    /*          if(VERB) printf("IRIS area center lon %f lat %f\n",IRISgeo.u/DEG_TO_RAD,IRISgeo.v/DEG_TO_RAD); */
	  projcode=irishdr[158];
	  memcpy(&IRPar1,irishdr+Base+212,4);
	  memcpy(&IRPar2,irishdr+Base+216,4);
	  memcpy(&ERad_cm,irishdr+Base+220,4);
	  if(!ERad_cm) ERad=6371000.0; else ERad=(double)ERad_cm/100.0;
	  memcpy(&IRFlatten,irishdr+Base+224,4);
	  if(IRFlatten) ERflatten=1.0/((double)IRFlatten/1e+6); else ERflatten=0.0;
	  memcpy(&Area_reflat,irishdr+Base+240,4);
	  memcpy(&Area_reflon,irishdr+Base+244,4);
	  projtype=projtypes[projcode];

	  lon_0 = 180.0*((double)Area_reflon/((double)INT32_MAX + 1));
	  lat_0 = 180.0*((double)Area_reflat/((double)INT32_MAX + 1));
	  std_par1 = 180.0*((double)IRPar1/((double)INT32_MAX + 1));
	  std_par2 = 180.0*((double)IRPar2/((double)INT32_MAX + 1));
      
          if(OPTCENTER) { lon_0 = IRISgeo.u/DEG_TO_RAD; lat_0 = IRISgeo.v/DEG_TO_RAD; }
	  sprintf(InLUT,"IRIS:%s_%dx%d_%s_P%d,%d_R%d,%d_C%d,%d_a%u_f%u",sitecode,IRISdim[0],IRISdim[1],
		  projtype,IRISres[0],IRISres[1],Area_reflon,Area_reflat,IRIScentlon,IRIScentlat,ERad_cm,IRFlatten);

	  /* printf("%s\n",InLUT); */


	  /* set the input projection according to IRIS projection definition */
	  sprintf(IRIS_projdef,"+proj=%s +a=%.6f +f=%.16f ",projtype,ERad,ERflatten);
	  switch(projcode)
	  {
	     case(P_aeqd): case(P_gnom):
	       sprintf(IRIS_projdef,"%s +lon_0=%.6f +lat_0=%.6f",
		       IRIS_projdef,lon_0,lat_0); break;

	     case(P_merc): case(P_eqc):
	       sprintf(IRIS_projdef,"%s +lon_0=%.6f +lat_ts=%.6f",
		       IRIS_projdef,lon_0,lat_0); break;

	     case(P_stere):
	       if(lat_0 == 0) lat_0=90;
	         sprintf(IRIS_projdef,"%s +lon_0=%.6f +lat_0=90N +lat_ts=%.6f",
		 IRIS_projdef,lon_0,lat_0); break;

	     case(P_utm):
	       sprintf(IRIS_projdef,"%s +lon_0=%.6f",
		       IRIS_projdef,lon_0); break;

	     case(P_nsper):
	       sprintf(IRIS_projdef,"%s +lon_0=%.6f +h=%.0f",
		       IRIS_projdef,lon_0,35786000.0); break;

	     case(P_lcc):
	       sprintf(IRIS_projdef,"%s +lon_0=%.6f +lat_0=%.6f +lat_1=%.6f +lat_2=%.6f ",
		       IRIS_projdef,lon_0,lat_0,std_par1,std_par2); break;

	     default: 
	       if(!QUIET) printf("Unknown IRIS projection code: %d\n",projcode); return(1);
	  }

       if(VERB) printf("IRIS projdef: %s\n",IRIS_projdef);
       if(!projdef[IN][0]) { sprintf(projdef[IN],"%s",IRIS_projdef); 
       if(!QUIET) printf("IRIS projection used\n"); }
       projref[IN]=pj_init_plus(projdef[IN]);
       coords[GEO][IN][CE] = IRISgeo; /* IRIS projection center trusted, 
                                         cartesian corners recalculated after geographical corner iteration */
       /*       coords[CART][IN][CE] = pj_fwd(IRISgeo, projref[IN]); */    

       sprintf(cornersname,"%s/%s.geocorners",LUTdir,InLUT);
       OFFSFILE = fopen(cornersname,"r");
       if(OFFSFILE) reti=fscanf(OFFSFILE,"%lf %lf %lf %lf",&coords[GEO][IN][SW].u,&coords[GEO][IN][SW].v,
                    &coords[GEO][IN][NE].u,&coords[GEO][IN][NE].v);
       else
       {
       /* IRIS ranges are expressed along surface of 6371 km sphere beginning from rectangular image
	  centre towards Y and X in cartesian projection coordinates having reference lon_0 and lat_0. 
          Using these, the geographical corners are iterated by give_IRIS_geocorners function.
	  Cartesian coordinates of corners are then calculated using proj with actual projection.
          The center may be recalculated, because it differs from IRIS definition except with
          gnomonic projection. */

          give_IRIS_geocorners(IRIScorners,lon_0,lat_0,ERad,ERflatten,projtype,IRISgeo,IRISrange);
          coords[GEO][IN][SW]=IRIScorners[SW];
          coords[GEO][IN][NE]=IRIScorners[NE];
	  OFFSFILE = fopen(cornersname,"w");
	  if(!QUIET) printf("Writing IRIS area corner geocoords to %s\n",cornersname);
	  fprintf(OFFSFILE,"%.12f %.12f %.12f %.12f",IRIScorners[SW].u,IRIScorners[SW].v,IRIScorners[NE].u,IRIScorners[NE].v);
       }
       fclose(OFFSFILE);

       coords[CART][IN][SW] = pj_fwd(coords[GEO][IN][SW],projref[IN]);
       coords[CART][IN][NE] = pj_fwd(coords[GEO][IN][NE],projref[IN]);
       inrange[XD] = coords[CART][IN][NE].u - coords[CART][IN][SW].u;
       inrange[YD] = coords[CART][IN][NE].v - coords[CART][IN][SW].v;
       for(Co=XD;Co<=YD;Co++) inres[Co]=inrange[Co]/indim[Co];
       coords[CART][IN][CE].u=coords[CART][IN][SW].u+inrange[XD]/2.0;  
       coords[CART][IN][CE].v=coords[CART][IN][SW].v+inrange[YD]/2.0;
 
       if(!IGN_INCE) USE_INCE = 1; /* In case no output georeference given */
  }

  /* Calculating the center geocoords if not IRIS data */
  if(!USE_INCE)
  {
     coords[CART][IN][CE].u=coords[CART][IN][SW].u+inrange[XD]/2.0;  
     coords[CART][IN][CE].v=coords[CART][IN][SW].v+inrange[YD]/2.0;
     coords[GEO][IN][CE]=pj_inv(coords[CART][IN][CE],projref[IN]);  
  }

  coords[CART][IN][NW].u=coords[CART][IN][SW].u;  
  coords[CART][IN][NW].v=coords[CART][IN][SW].v+inrange[YD];  
  coords[CART][IN][SE].u=coords[CART][IN][NE].u;  
  coords[CART][IN][SE].v=coords[CART][IN][NE].v-inrange[YD];  
  coords[GEO][IN][NW]=pj_inv(coords[CART][IN][NW],projref[IN]);
  coords[GEO][IN][SE]=pj_inv(coords[CART][IN][SE],projref[IN]);


  /* Setting the cartesian offsets for input data SW corner (= common origo for inout and output) */ 
  inXoff=-coords[CART][IN][SW].u;
  inYoff=-coords[CART][IN][SW].v;

  if(indim[XD] == 0) { if(!QUIET) printf("Give inXdim for input data\n"); return(1); }
  if(indim[YD] == 0) { if(!QUIET) printf("Give inYdim for input data\n"); return(1); }

  if(VERB)
  {
     for(Co=SW;Co<=SE;Co++)
     {
        DEG[Co].u = coords[GEO][IN][Co].u/DEG_TO_RAD;
        DEG[Co].v = coords[GEO][IN][Co].v/DEG_TO_RAD;
     }

     printf("Input datatype: %s\n\n",indatatype);
     printf("inXdim %11d, inYdim %11d\n",indim[XD],indim[YD]);
     printf("inXrange %9.3f, inYrange %9.3f\n",inrange[XD]/1000.0,inrange[YD]/1000.0);
     printf("inXres %11.3f, inYres %11.3f\n\n",inres[XD],inres[YD]);

     printf("inNWlon %10.6f (%s), ",DEG[NW].u,toDMS(DEG[NW].u)); 
     printf("inNWlat %10.6f (%s)\n",DEG[NW].v,toDMS(DEG[NW].v)); 
     printf("inSWlon %10.6f (%s), ",DEG[SW].u,toDMS(DEG[SW].u)); 
     printf("inSWlat %10.6f (%s)\n",DEG[SW].v,toDMS(DEG[SW].v)); 
     printf("inNElon %10.6f (%s), ",DEG[NE].u,toDMS(DEG[NE].u)); 
     printf("inNElat %10.6f (%s)\n",DEG[NE].v,toDMS(DEG[NE].v)); 
     printf("inSElon %10.6f (%s), ",DEG[SE].u,toDMS(DEG[SE].u)); 
     printf("inSElat %10.6f (%s)\n",DEG[SE].v,toDMS(DEG[SE].v)); 
     printf("inCElon %10.6f (%s), ",DEG[CE].u,toDMS(DEG[CE].u)); 
     printf("inCElat %10.6f (%s)\n\n",DEG[CE].v,toDMS(DEG[CE].v)); 

     printf("inNWx %12.3f, inNWy %12.3f\n",coords[CART][IN][NW].u,coords[CART][IN][NW].v);
     printf("inSWx %12.3f, inSWy %12.3f\n",coords[CART][IN][SW].u,coords[CART][IN][SW].v);
     printf("inNEx %12.3f, inNEy %12.3f\n",coords[CART][IN][NE].u,coords[CART][IN][NE].v);
     printf("inSEx %12.3f, inSEy %12.3f\n",coords[CART][IN][SE].u,coords[CART][IN][SE].v);
     printf("inCEx %12.3f, inCEy %12.3f\n",coords[CART][IN][CE].u,coords[CART][IN][CE].v);

     printf("\nInproj %s\n\n",projdef[IN]); 
  }      

  /* Read input data */ 
  if(!HDF5) /* If input data is HDF5, the inarr is already filled when cartesian read or during polar to cartesian conversion */
  {
     insize=indim[XD]*indim[YD]*bytes*Chans;
     inarr=malloc(insize);
     readsize=fread(inarr,insize,1,INF);
     fclose(INF);
  }
  /* Swab byte order if 2 byte PNM data and little endian host, because PNM is always big endian */
  if(bytes==2 && PNM && HOST_LE) swab(inarr,inarr,insize);

  /*================================== OUTPUT GEOMETRY ====================================================================*/

  if(VERB) printf("OUTPUT GEOMETRY\n==============\n\nOutput initial geometry: ");
  /* set output projection reference */
  if(!projdef[OUT][0])
  {
    if(!QUIET) printf("Give option --outproj=PROJ4 definition string!\n");
    return(1);
  }
  if(USE_INCE) 
  {
    char tproj[300];

    sprintf(tproj,"+lon_0=%f +lat_0=%f %s",coords[GEO][IN][CE].u/DEG_TO_RAD,coords[GEO][IN][CE].v/DEG_TO_RAD,projdef[OUT]);
    sprintf(projdef[OUT],"%s",tproj);
    if(VERB) printf("Because using input center, added center to outproj: %s\n",projdef[OUT]);
  }
  projref[OUT]=pj_init_plus(projdef[OUT]);

  /* check the georeferences first */

  if(OUT_BBOX_CART)
  {
     coords[GEO][OUT][SW] = pj_inv(coords[CART][OUT][SW],projref[OUT]);  
     coords[GEO][OUT][NE] = pj_inv(coords[CART][OUT][NE],projref[OUT]);
  }

  if(USE_INSW) coords[GEO][OUT][SW]=coords[GEO][IN][SW];
  if(USE_INNE) coords[GEO][OUT][NE]=coords[GEO][IN][NE]; 
  
  if(coords[GEO][OUT][SW].u<4 && coords[GEO][OUT][SW].v<4) { coords[CART][OUT][SW]=pj_fwd(coords[GEO][OUT][SW],projref[OUT]); SW_CORNER=1; } 
  if(coords[GEO][OUT][NE].u<4 && coords[GEO][OUT][NE].v<4) { coords[CART][OUT][NE]=pj_fwd(coords[GEO][OUT][NE],projref[OUT]); NE_CORNER=1; } 
  if(coords[GEO][OUT][CE].u < 4 && coords[GEO][OUT][CE].v < 4) OUTCENTER = 1;

  if(!SW_CORNER && !NE_CORNER)
  {
     if(!OUTCENTER)
     { 
       coords[GEO][OUT][CE]=coords[GEO][IN][CE]; 
       if(!QUIET) printf("No output georeference given, using input center coordinates as reference!\n");
     }
     coords[CART][OUT][CE]=pj_fwd(coords[GEO][OUT][CE],projref[OUT]); 
     OUTCENTER = 1; 
  }

  if(SW_CORNER && NE_CORNER)
  {
     if(!QUIET) printf("Both output corners given, overriding possible ranges given as options!\n");
     outrange[XD]=coords[CART][OUT][NE].u-coords[CART][OUT][SW].u;
     outrange[YD]=coords[CART][OUT][NE].v-coords[CART][OUT][SW].v;
     /* Center calculated just for information */
     coords[CART][OUT][CE].u = coords[CART][OUT][SW].u+outrange[XD]/2;
     coords[CART][OUT][CE].v = coords[CART][OUT][SW].v+outrange[YD]/2;
     coords[GEO][OUT][CE] = pj_inv(coords[CART][OUT][CE],projref[OUT]);
  }

  /* Check the consistency of output geometry parameters */ 
  for(Co=XD;Co<=YD;Co++)
  {
      int range,dim,res;

      range=(outrange[Co]>0);
      dim=(outdim[Co]>0);
      res=(outres[Co]>0);

      while(1)
      {
        /* case of all three outbut geometry parameters are given  */
        if(dim && range && res)
        {
            if(!QUIET) printf("All output geometry parameters given. Give only two or one!\n");
            return(1);
        }
        /* cases of two output geometry parameters are given (minimum to calculate the third one) */
        if(dim && range && !res)
          { outres[Co]=outrange[Co]/(double)outdim[Co]; break; }

        if(dim && !range && res)
          { outrange[Co]=(double)outdim[Co]*outres[Co]; break; }

        if(!dim && range && res)
          { outdim[Co]=outrange[Co]/outres[Co];
            outres[Co]=outrange[Co]/(double)outdim[Co]; break; }

        /* fall to cases when only one parameter is given. The corresponding input parameter will be forced to be used */
        if(dim)
          { outrange[Co]=inrange[Co];
     	    outres[Co]=outrange[Co]/(double)outdim[Co]; break; }

        if(res)
          { outrange[Co]=inrange[Co];
            outdim[Co]=outrange[Co]/outres[Co]; break; }

        if(range)
          { outdim[Co]=outrange[Co]/inres[Co];
            outres[Co]=outrange[Co]/(double)outdim[Co]; break; }

        /* fall to case when output geometry parameters not given at all */ 
        if(!QUIET) printf("No output %c-coordinate geometry parameters given, using input range!\n",coord[Co]);
        if(!QUIET) printf("Output dimension will be calculated, and resolution recalculated!\n");
        outres[Co]=inres[Co];
        outrange[Co]=inrange[Co];
        outdim[Co]=(int)(outrange[Co]/outres[Co]);
        outres[Co]=outrange[Co]/(double)outdim[Co];
        break;
      }
  }
  /*        printf("OUTRANGE XY %f %f\n",outrange[XD],outrange[YD]); */

  /* Calculating the remaining output coordinates according to output geometry. 
     Corners has the priority over center */

  if(SW_CORNER ^ NE_CORNER) /* Both defined handled before */
  {
     if(SW_CORNER) 
     {
        if(!QUIET) printf("SW corner only\n");
        coords[CART][OUT][NE].u = coords[CART][OUT][SW].u + outrange[XD];
        coords[CART][OUT][NE].v = coords[CART][OUT][SW].v + outrange[YD];
        coords[GEO][OUT][NE] = pj_inv(coords[CART][OUT][NE],projref[OUT]);
     }

     if(NE_CORNER) 
     {
        if(!QUIET) printf("NE corner only\n");
        coords[CART][OUT][SW].u = coords[CART][OUT][NE].u - outrange[XD];
        coords[CART][OUT][SW].v = coords[CART][OUT][NE].v - outrange[YD];
        coords[GEO][OUT][SW] = pj_inv(coords[CART][OUT][SW],projref[OUT]);
     }

     coords[CART][OUT][CE].u = coords[CART][OUT][SW].u + outrange[XD]/2.0;
     coords[CART][OUT][CE].v = coords[CART][OUT][SW].v + outrange[YD]/2.0;

  } else if(!SW_CORNER && !NE_CORNER) /* Only center defined */
    {    
       coords[CART][OUT][SW].u = coords[CART][OUT][CE].u - outrange[XD]/2;
       coords[CART][OUT][SW].v = coords[CART][OUT][CE].v - outrange[YD]/2;
       coords[CART][OUT][NE].u = coords[CART][OUT][CE].u + outrange[XD]/2;
       coords[CART][OUT][NE].v = coords[CART][OUT][CE].v + outrange[YD]/2;

       coords[GEO][OUT][SW] = pj_inv(coords[CART][OUT][SW],projref[OUT]);
       coords[GEO][OUT][NE] = pj_inv(coords[CART][OUT][NE],projref[OUT]);
    }

  if(!OUTCENTER) coords[GEO][OUT][CE] = pj_inv(coords[CART][OUT][CE],projref[OUT]);

  coords[CART][OUT][NW].u = coords[CART][OUT][SW].u;
  coords[CART][OUT][NW].v = coords[CART][OUT][SW].v + outrange[YD];
  coords[CART][OUT][SE].u = coords[CART][OUT][SW].u + outrange[XD];
  coords[CART][OUT][SE].v = coords[CART][OUT][SW].v;
  
  coords[GEO][OUT][NW] = pj_inv(coords[CART][OUT][NW],projref[OUT]);  
  coords[GEO][OUT][SE] = pj_inv(coords[CART][OUT][SE],projref[OUT]);

  /* List the recalculated output geometry parameters */
  if(VERB)
  {
     for(Co=SW;Co<=SE;Co++)
     {
        DEG[Co].u = coords[GEO][OUT][Co].u/DEG_TO_RAD;
        DEG[Co].v = coords[GEO][OUT][Co].v/DEG_TO_RAD;
     }

     printf("\noutXdim %10d, outYdim %10d\n",outdim[XD],outdim[YD]);
     printf("outXrange %8.3f, outYrange %8.3f\n",outrange[XD]/1000.0,outrange[YD]/1000.0);
     printf("outXres %10.3f, outYres %10.3f\n\n",outres[XD],outres[YD]);

     printf("outNWlon %10.6f (%s), ",DEG[NW].u,toDMS(DEG[NW].u)); 
     printf("outNWlat %10.6f (%s)\n",DEG[NW].v,toDMS(DEG[NW].v)); 
     printf("outSWlon %10.6f (%s), ",DEG[SW].u,toDMS(DEG[SW].u)); 
     printf("outSWlat %10.6f (%s)\n",DEG[SW].v,toDMS(DEG[SW].v)); 
     printf("outNElon %10.6f (%s), ",DEG[NE].u,toDMS(DEG[NE].u)); 
     printf("outNElat %10.6f (%s)\n",DEG[NE].v,toDMS(DEG[NE].v)); 
     printf("outSElon %10.6f (%s), ",DEG[SE].u,toDMS(DEG[SE].u)); 
     printf("outSElat %10.6f (%s)\n",DEG[SE].v,toDMS(DEG[SE].v)); 
     printf("outCElon %10.6f (%s), ",DEG[CE].u,toDMS(DEG[CE].u)); 
     printf("outCElat %10.6f (%s)\n\n",DEG[CE].v,toDMS(DEG[CE].v)); 

     printf("outNWx %12.3f, outNWy %12.3f\n",coords[CART][OUT][NW].u , coords[CART][OUT][NW].v ); 
     printf("outSWx %12.3f, outSWy %12.3f\n",coords[CART][OUT][SW].u , coords[CART][OUT][SW].v ); 
     printf("outNEx %12.3f, outNEy %12.3f\n",coords[CART][OUT][NE].u , coords[CART][OUT][NE].v ); 
     printf("outSEx %12.3f, outSEy %12.3f\n",coords[CART][OUT][SE].u , coords[CART][OUT][SE].v ); 
     printf("outCEx %12.3f, outCEy %12.3f\n\n",coords[CART][OUT][CE].u , coords[CART][OUT][CE].v ); 
     printf("Outproj\t%s +units=m\n",projdef[OUT]);
     printf("ullr %.2f %.2f %.2f %.2f\n",coords[CART][OUT][NW].u, coords[CART][OUT][NW].v,coords[CART][OUT][SE].u, coords[CART][OUT][SE].v);  
  }

  /* Define half pixel widths to get sample to center of output pixel */
  halfXres=outres[XD]/2.0;
  halfYres=outres[YD]/2.0;

  /* Define output data size in bytes */ 
  outsize=outdim[XD]*outdim[YD]*bytes*Chans;
  datarr=malloc(outsize);
  memset(datarr,255,outsize);

  /*========================  READING / GENERATING PROJECTION CONVERSION LUT =========================================*/

  lutsize=outdim[YD]*outdim[XD]*sizeof(long);
  LUTarr=malloc(lutsize);
  memset(LUTarr,255,lutsize);    
  LUTF=NULL;
  LUT=0;
  for(i=j=0;i<strlen(projdef[IN]);i++) if(projdef[IN][i]!=' ') inlutproj[j++]=projdef[IN][i];
  for(i=j=0;i<strlen(projdef[OUT]);i++) if(projdef[OUT][i]!=' ') outlutproj[j++]=projdef[OUT][i];
  
  if(lutname[0]==0)
  {
     if(IRISDATA) sprintf(lutname,"%s/%s",LUTdir,InLUT); 
     else sprintf(lutname,"%s/I_%s:%dx%d",LUTdir,sitecode,indim[XD],indim[YD]);
     /*     else sprintf(lutname,"%s/I_%s:%dx%d_%.3f,%.3f_%.3f,%.3f%s",LUTdir,sitecode,
          indim[XD],indim[YD],coords[GEO][IN][SW].u/DEG_TO_RAD,coords[GEO][IN][SW].v/DEG_TO_RAD,
          coords[GEO][IN][NE].u/DEG_TO_RAD,coords[GEO][IN][NE].v/DEG_TO_RAD,inlutproj); */
  }

  { /* Using tmpname because it depends on gcc version if catenating to same variable is done by sprintf */
    char tmpname[500];
    sprintf(tmpname,"%s",lutname);

    sprintf(lutname,"%s_O:%dx%d_%.3f,%.3f_%.3f,%.3f%s_B%dC%ld.lut",tmpname,
          outdim[XD],outdim[YD],coords[GEO][OUT][SW].u/DEG_TO_RAD,coords[GEO][OUT][SW].v/DEG_TO_RAD,
          coords[GEO][OUT][NE].u/DEG_TO_RAD,coords[GEO][OUT][NE].v/DEG_TO_RAD,outlutproj,bytes,Chans);
  }


  LUTF=fopen(lutname,"r");
  if(LUTF) 
  {
     if(VERB) printf("\nReading projection conversion LUT %s\n\n",lutname);
     readsize=fread(LUTarr,lutsize,1,LUTF);
     LUT = 1;
  }   
  else 
  {  
     LUTF=fopen(lutname,"w");
     if(VERB) printf("\nGenerating projection conversion LUT %s\n\n",lutname);
     LUT = 0;
  }

  /*================================== CONVERSION LOOP ====================================================================*/

  maxoutN=indim[XD]*indim[YD]*Chans;
  for(oY=0;oY<outdim[YD];oY++)
  {
     if(!LUT) tocoord.v = outrange[YD] - ((double)oY*outres[YD]+halfYres) + coords[CART][OUT][SW].v; 
     for(oX=0;oX<outdim[XD];oX++)
     {
        if(FLIP_OUTPUT_DATA) oN = (outdim[YD]-oY-1)*outdim[XD] + oX;
        else oN = oY*outdim[XD] + oX;

        if(!LUT)
        {
           tocoord.u = (double)oX*outres[XD] + halfXres + coords[CART][OUT][SW].u;
           lonlat = pj_inv(tocoord,projref[OUT]);
           fromcoord = pj_fwd(lonlat,projref[IN]);
           X = (int)((fromcoord.u + inXoff) / inres[XD]);
           relY = (int)((fromcoord.v + inYoff) / inres[YD]);
           if(PNM || CARTESIAN_HDF5) Y = indim[YD] - relY; else Y = relY;

        /* if outside input image */
           if(X<0 || X>=indim[XD] || Y<0 || Y>=indim[YD]) continue;
     
           N = (Y*indim[XD] + X)*Chans;
           LUTarr[oN]=N;
	}  else N=LUTarr[oN];

        if(N<0) continue;
        if(N>=maxoutN) continue;

        cN = oN*Chans;

	/* datarr is array between input data and output data. Byte order is always host byte order 
           datarr is used directly for output if there is no further conversion needed */
        for(C=0;C<Chans;C++)
        {
           if(bytes==1) datarr[cN+C] = inarr[N+C];
           else memcpy(&datarr[(cN+C)*2],&inarr[(N+C)*2],2);
        }
     }  
  }

  /* ------------------------- Output data definitions and final conversions ---------------------------------- */


  if((ODIM_quantcode==ODIM_DBZH || ODIM_quantcode==ODIM_TH || ODIM_quantcode==ODIM_DBZHC || ODIM_quantcode==ODIM_DBZH_QC ) &&  DBZ_TO_RATE)
  {
    uint16_t ival,*ratearr,RR;
    size_t ratesize,Oi,Ri;
    double dBZ,fRR,effUndetect,effNodata,effGain,effOffset;

    if(Chans!=1) { if(!QUIET) printf("Conversion of dBZ to rain rate for colour images is not possible!\n"); return(1); }
    ratesize=outsize;
    if(bytes==1)
    { 
       ratesize *= 2;
       effUndetect=Gundetect;
       effNodata=Gnodata;
       effGain=Ggain;
       effOffset=Goffset;
    }

    if(bytes==2)
    { 
       effUndetect=inUndetect;
       effNodata=inNodata;
       effGain=inGain;
       effOffset=inOffset;
    }


    /*    printf("Bytes %d, outsize %ld, Gain %f Offset %f\n",bytes,outsize,effGain,effOffset); */

    ratearr=calloc(ratesize,1);
    for(Ri=Oi=0;Oi<outsize;Oi+=bytes,Ri++)
    {
       if(bytes==1) ival=datarr[Oi];
       if(bytes==2) memcpy(&ival,&datarr[Oi],bytes);
       /*       if(ival==effUndetect) { ratearr[Ri]=0; continue; } */
       if(ival==effUndetect) continue;
       if(ival==effNodata) 
       {
	  if(!FILL_NODATA) ratearr[Ri]=UINT16_MAX; 
          continue; 
       }
       dBZ = effGain*(double)ival + effOffset;
       fRR = Rfactor * dBZtoR(dBZ);
       /*       printf("%d\t%.2f\t%.2f\n",ival,dBZ,fRR); */
       ratearr[Ri]=(uint16_t)(100.0*fRR);
    } 

    if(bytes==1) datarr=realloc(datarr,ratesize);
    memcpy(datarr,ratearr,ratesize);
    free(ratearr);

    bytes=2;
    outsize=ratesize;
    Gnodata=UINT16_MAX;
    Gundetect=0;
    Ggain=0.01;
    Goffset=0.0;
    sprintf(ODIMquantity,"RATE"); 
    ODIM_quantcode = ODIM_RATE; 
  }

  /* PNM output  ---------------------------------------------------------------------------------------------*/ 
  if(PNMOUT)
  {
     uint8_t *outarr=NULL;

     outarr=malloc(outsize);
     memset(outarr,255,outsize);
     if(bytes>1 && HOST_LE) 
     { 
        /* Swabbing because PGM is always big endian */ 
       swab(datarr,outarr,outsize);
     } else memcpy(outarr,datarr,outsize);

     OUTF=fopen(pnmoutfile,"w");
     if(QDHEADERS && strstr(projdef[OUT],"=stere ")) 
     {
       char *p;

       p=strstr(projdef[OUT],"lon_0=");
       if(p) { p=strchr(p,'=')+1; lon_0=atof(p); }

       p=strstr(projdef[OUT],"lat_ts=");
       if(p) { p=strchr(p,'=')+1; lat_0=atof(p); }

       if(lon_0 > 180 || lat_0 > 180) { if(!QUIET) printf("Central coordinates of PNM data for QueryData not defined!\n"); return(1); }
       /* add # param PrecipitationRate etc */
       if(!QDparam[0])
       {
          if(ODIM_quantcode==ODIM_RATE || DBZ_TO_RATE) sprintf(QDparam,"PrecipitationRate");
          if(ODIM_quantcode==ODIM_DBZH || RATE_TO_DBZ) sprintf(QDparam,"CorrectedReflectivity");
          if(ODIM_quantcode==ODIM_HCLASS) sprintf(QDparam,"HydroClass");
       }

       fprintf(OUTF,"%s\n# obstime %s\n# param %s\n# projection radar {\n# type stereographic\n",outpnmcode,obstime,QDparam); 
       fprintf(OUTF,"# centrallongitude %.6f\n# centrallatitude 90.0\n# truelatitude %.6f\n",lon_0,lat_0); 
       fprintf(OUTF,"# bottomleft %.6f %.6f\n# topright %.6f %.6f\n",coords[GEO][OUT][SW].u/DEG_TO_RAD, coords[GEO][OUT][SW].v/DEG_TO_RAD, 
               coords[GEO][OUT][NE].u/DEG_TO_RAD, coords[GEO][OUT][NE].v/DEG_TO_RAD); 
       fprintf(OUTF,"# }\n%d %d\n%d\n",outdim[XD],outdim[YD],(int)(pow(256,bytes))-1);       
     } 
     else
     {
          fprintf(OUTF,"%s\n# SWlon=%.7f SWlat=%.7f\n# NElon=%.7f NElat=%.7f\n",
          outpnmcode,coords[GEO][OUT][SW].u/DEG_TO_RAD,coords[GEO][OUT][SW].v/DEG_TO_RAD,coords[GEO][OUT][NE].u/DEG_TO_RAD,
          coords[GEO][OUT][NE].v/DEG_TO_RAD);
          fprintf(OUTF,"# bottomleft %.6f %.6f\n# topright %.6f %.6f\n",coords[GEO][OUT][SW].u/DEG_TO_RAD, coords[GEO][OUT][SW].v/DEG_TO_RAD, 
               coords[GEO][OUT][NE].u/DEG_TO_RAD, coords[GEO][OUT][NE].v/DEG_TO_RAD); 
          fprintf(OUTF,"%d %d\n%d\n",outdim[XD],outdim[YD],(int)(pow(256,bytes))-1);
     }

     fwrite(outarr,outsize,1,OUTF);   
     fclose(OUTF);
     free(outarr);
  }


  /* GeoTIFF output  --------------------------------------------------------------------------------*/ 
  if(GTIFF) 
  { 
# ifndef NO_GEOTIFF
     GDALDriverH hDriver;
     GDALDatasetH hDstDS;        
     char **papszOptions = NULL;
     double adfGeoTransform[6];
     OGRSpatialReferenceH hSRS;
     char *pszSRS_WKT = NULL;
     GDALRasterBandH hBand;
     int GDT_datatype;
     char projstr[400]={0};
     char **papszMetaData=NULL;
     char tmptag[500]={0},metatag[100000]={0};
     QuantityStrings Qstr;
     struct stat gbuf;
     int nonexist,ii;
     CPLErr OVerr;
     char temp_path[300];

     /* Define geometry, projection and compression */
     sprintf(projstr,"%s +units=m +no_defs",projdef[OUT]); 

     adfGeoTransform[0] = coords[CART][OUT][NW].u;
     adfGeoTransform[3] = coords[CART][OUT][NW].v;
     adfGeoTransform[1] = outres[XD];
     adfGeoTransform[5] = -outres[YD];
     adfGeoTransform[2] = 0;
     adfGeoTransform[4] = 0;

     GDALAllRegister();

     GDT_datatype=bytes;
     hDriver = GDALGetDriverByName("GTiff");

     if(strcmp(TIFFcompression,"NONE")) 
     {
         papszOptions = CSLSetNameValue( papszOptions,"COMPRESS",TIFFcompression);
         if(Zlevel) papszOptions = CSLSetNameValue( papszOptions,"ZLEVEL",sZlevel);
     }
     if(tiling)
     {
         papszOptions = CSLSetNameValue( papszOptions,"TILED","YES");
         papszOptions = CSLSetNameValue( papszOptions,"BLOCKXSIZE",Xtile);
         papszOptions = CSLSetNameValue( papszOptions,"BLOCKYSIZE",Ytile);
     }

     /* Creating temporary filename */
     {
        char *tpath,estr[100],*tname;

        tpath=tmpnam(estr);
        tname=strrchr(tpath,'/');
        if(!tname) tname=tpath; else tname+=1;
        sprintf(temp_path,"%s.%s",GTiffoutfile,tname);
     }


     /* Creating the temporary GeoTIFF file */

     hDstDS = GDALCreate( hDriver, temp_path, outdim[XD], outdim[YD], 1, GDT_datatype, papszOptions );

     /* Insert metadata as XML-tags */
     if(INSERT_METAFILE || ODIM_quantcode>=0)
     {
         FILE *TMPF; /* Had to use this, because sprintf(str,"%s%s",str,...) doesn's work in every system */
        
         TMPF=tmpfile();
         fprintf(TMPF,"<GDALMetadata>\n");

	 if(INSERT_METAFILE)
	 {
	     FILE *METAF;
             char linebuf[500];
    
             METAF=fopen(metafile,"r");
             while(1)
	     {
	        memset(linebuf,0,500);
                fgets(linebuf,499,METAF);
                if(feof(METAF)) break;
                fprintf(TMPF,"%s",linebuf);
	     }
             fclose(METAF);
    	 }
          	 
         if(ODIM_quantcode>=0 && DEFAULTMETA)
	 {

	     Qstr=QuantStr[ODIM_quantcode];
	     if(obstime[0])
	     {
	       fprintf(TMPF,"<Item name=\"Observation time\" format=\"YYYYMMDDhhmm\">%s</Item>\n",obstime);
	       if(Forecast && (ODIM_quantcode==ODIM_DBZH || ODIM_quantcode==ODIM_RATE))  
		  fprintf(TMPF,"<Item name=\"Forecast time\" format=\"YYYYMMDDhhmm\">%s</Item>\n",fortime_start);
	     }
	     if(ODIM_quantcode==ODIM_ACCPROB)
	     {
	       fprintf(TMPF,"<Item name=\"Forecast start time\" format=\"YYYYMMDDhhmm\">%s</Item>\n",fortime_start);
	       fprintf(TMPF,"<Item name=\"Forecast end time\" format=\"YYYYMMDDhhmm\">%s</Item>\n",fortime_end);
	       fprintf(TMPF,"<Item name=\"Time zone\">UTC</Item>\n");
	       fprintf(TMPF,"<Item name=\"Accumulation time\" unit=\"min\">%d</Item>\n",accmins);
	       fprintf(TMPF,"<Item name=\"Accumulation threshold\" unit=\"mm\">%.2f</Item>\n",acc_threshold);
	     }

	     fprintf(TMPF,"<Item name=\"Quantity\" unit=\"%s\">%s</Item>\n",Qstr.unit,Qstr.GTIFF_qstr);
	     fprintf(TMPF,"<Item name=\"Gain\">%f</Item>\n",Ggain);
	     fprintf(TMPF,"<Item name=\"Offset\">%f</Item>\n",Goffset);
	     fprintf(TMPF,"<Item name=\"Nodata\">%d</Item>\n",Gnodata);
	     fprintf(TMPF,"<Item name=\"Undetect\">%d</Item>\n",Gundetect);

	     if(ODIM_quantcode == ODIM_VRAD) 
		fprintf(TMPF,"<Item name=\"Maximum unambiguous velocity\" unit=\"m/s\">%.2f</Item>\n",NI);

	     if(ODIM_quantcode == ODIM_ETOP) 
		fprintf(TMPF,"<Item name=\"Echo top reflectivity threshold\" unit=\"dBZ\">%.1f</Item>\n",TOPthr);

	     if(ODIM_quantcode == ODIM_HCLASS || ODIM_quantcode == ODIM_HCLASSBG)
	     { 
		fprintf(TMPF,"<Item name=\"Non-meteorological echo class\">1</Item>\n");
		fprintf(TMPF,"<Item name=\"Rain class\">2</Item>\n");
		fprintf(TMPF,"<Item name=\"Wet snow class\">3</Item>\n");
		fprintf(TMPF,"<Item name=\"Dry snow class\">4</Item>\n");
		fprintf(TMPF,"<Item name=\"Graupel class\">5</Item>\n");
		fprintf(TMPF,"<Item name=\"Hail class\">6</Item>\n");
	     }

	     if(ODIM_quantcode == ODIM_ACRR)
	     {
		fprintf(TMPF,"<Item name=\"Accumulation time\" unit=\"h\">%d</Item>\n",Acchours);
		if(Forecast) fprintf(TMPF,"<Item name=\"Temporal type\">Forecast</Item>\n"); 
		else fprintf(TMPF,"<Item name=\"Temporal type\">Past</Item>\n");
	     }

	     if(ODIM_quantcode == ODIM_RATE)
	     {
		fprintf(TMPF,"<Item name=\"Adjustment factor of precipitation intensity\">%f</Item>\n",Rfactor);
	     }

	     if(ODIM_quantcode == ODIM_ACRR || ODIM_quantcode == ODIM_RATE)
	     {
		fprintf(TMPF,"<Item name=\"Z=AR^B relation parameter A\">%f</Item>\n",ZRA);
		fprintf(TMPF,"<Item name=\"Z=AR^B relation parameter B\">%f</Item>\n",ZRB);
                if(ZSA>0)
	        {
   		   fprintf(TMPF,"<Item name=\"Z=AS^B relation parameter A for snow\">%f</Item>\n",ZSA);
		   fprintf(TMPF,"<Item name=\"Z=AS^B relation parameter B for snow\">%f</Item>\n",ZSB);
		}
		if(conv_dBZlim < 100.0)
		{
		   fprintf(TMPF,"<Item name=\"Z=AR^B relation parameter A convective\">%f</Item>\n",ZRAc);
		   fprintf(TMPF,"<Item name=\"Z=AR^B relation parameter B convective\">%f</Item>\n",ZRBc);
		   fprintf(TMPF,"<Item name=\"Reflectivity threshold between frontal and convective rain\">%f</Item>\n",conv_dBZlim);
		   fprintf(TMPF,"<Item name=\"Intensity threshold between frontal and convective rain\">%f</Item>\n",conv_Rlim);
		}
	     }

	     if(POLAR_HDF5) 
	     {
		fprintf(TMPF,"<Item name=\"Elevation angle\" unit=\"deg\">%.1f</Item>\n",elangle/DEG_TO_RAD);
		fprintf(TMPF,"<Item name=\"Bin length\" unit=\"m\">%.0f</Item>\n",binlen);
		fprintf(TMPF,"<Item name=\"Bins per ray\">%ld</Item>\n",(long)raybins);
		fprintf(TMPF,"<Item name=\"Rays per scan\">%ld</Item>\n",(long)rays);
	     }

	 }
    
    	 fprintf(TMPF,"</GDALMetadata>\n");
	 rewind(TMPF);
	 while(1)
	 {
	    fgets(tmptag,499,TMPF);
	    if(feof(TMPF)) break;
	    strcat(metatag,tmptag);
	 }
	 
	 fclose(TMPF);
	 papszMetaData = CSLSetNameValue(papszMetaData,"GDAL_METADATA",metatag);
	 GDALSetMetadata(hDstDS,papszMetaData,NULL);
     }

     /* Insert geolocation information */
     GDALSetGeoTransform( hDstDS, adfGeoTransform );
     hSRS = OSRNewSpatialReference( NULL );
     if(EPSG) OSRImportFromEPSG(hSRS,EPSG);	
     else OSRImportFromProj4(hSRS,projstr);	
     OSRExportToWkt( hSRS, &pszSRS_WKT );
     OSRDestroySpatialReference( hSRS );
     GDALSetProjection( hDstDS, pszSRS_WKT );
     CPLFree( pszSRS_WKT );

     /* Build overviews */
     if(eff_OverviewList)
     { 
        OVerr = GDALBuildOverviews(hDstDS,overview_method,overviews,eff_OverviewList,0,NULL,NULL,NULL );
        if(OVerr == CE_Failure) GDALBuildOverviews(hDstDS,"NONE",overviews,eff_OverviewList,0,NULL,NULL,NULL );
     }

     /* Insert the data */
     hBand = GDALGetRasterBand( hDstDS, 1 );
     GDALRasterIO( hBand, GF_Write, 0, 0, outdim[XD], outdim[YD], 
                  datarr, outdim[XD], outdim[YD], GDT_datatype, 0, 0 );    

     GDALClose( hDstDS );
     CSLDestroy( papszOptions );

     /* Checking the GeoTIFF existence and renaming temporary file / overwriting old one */
     nonexist=stat(GTiffoutfile,&gbuf);
     if(nonexist || OVERWRITE) rename(temp_path, GTiffoutfile);
     else 
     {
         RetVal = 4;
         if(!QUIET) printf("GeoTIFF %s not created, exists already!\n",GTiffoutfile);
	 remove(temp_path);
     }
# else
     printf("GeoTIFF output not available! Compile with -DWITH_GEOTIFF if GDAL devel is installed.\n");
# endif
  }

  /* HDF5 output ---------------------------------------------------------------------------------------------*/
  if(HDFOUT)
  {
        hid_t OUTH5,whatg,whereg,howg,datasetg,datag,datawhatg,setwhatg,setwhereg,sethowg,dataset,plist,dataspace,datatype;
        hsize_t dims[2],chunkdims[2]={100,100};
        double fundetect,fnodata,startepochs,endepochs;
        int64_t xsize,ysize;
        time_t startsecs=0,endsecs=0;
        int ii;
        herr_t ret;
        char odate[9]="",otime[7]="";

        OUTH5 = H5Fcreate(HDFoutfile,H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
        H5LTset_attribute_string(OUTH5,"/","Conventions","ODIM_H5/V2_1");
	
        whatg=H5Gcreate2(OUTH5,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        whereg=H5Gcreate2(OUTH5,"where",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        howg=H5Gcreate2(OUTH5,"how",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        datasetg=H5Gcreate2(OUTH5,"dataset1",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        setwhatg=H5Gcreate2(datasetg,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        sethowg=H5Gcreate2(datasetg,"how",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        datag=H5Gcreate2(datasetg,"data1",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
        /* create new empty dataset and attributes to desination file */ 
        if(bytes == 1) datatype= H5T_STD_U8LE;
        if(bytes == 2) datatype= H5T_STD_U16LE;

        dims[0]=(hsize_t)outdim[YD];
        dims[1]=(hsize_t)outdim[XD];
        dataspace=H5Screate_simple(2, dims, NULL);
        plist = H5Pcreate(H5P_DATASET_CREATE);
        if(dims[0]>chunkdims[0] && dims[1]>chunkdims[1])
	{ 
            H5Pset_chunk(plist, 2, chunkdims);
            H5Pset_deflate( plist, 6);
	}
        dataset = H5Dcreate2(datag,"data", datatype, dataspace,
                  H5P_DEFAULT, plist, H5P_DEFAULT);
        H5Dwrite(dataset, datatype, H5S_ALL, dataspace, H5P_DEFAULT,datarr);
        H5LTset_attribute_string( datag, "data", "CLASS", "IMAGE");
        H5LTset_attribute_string( datag, "data", "IMAGE_VERSION", "1.2");

        datawhatg=H5Gcreate2(datag,"what",H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);	

        H5LTset_attribute_string(OUTH5,"dataset1/data1","IMAGE_VERSION","1.2");
        H5LTset_attribute_string(OUTH5,"dataset1/data1/what","quantity",QuantStr[ODIM_quantcode].ODIM_qstr);
        H5LTset_attribute_double(OUTH5,"dataset1/data1/what","gain",&Ggain,1);
        H5LTset_attribute_double(OUTH5,"dataset1/data1/what","offset",&Goffset,1);
        fnodata=(double)Gnodata;
        H5LTset_attribute_double(OUTH5,"dataset1/data1/what","nodata",&fnodata,1);
        fundetect=(double)Gundetect;
        H5LTset_attribute_double(OUTH5,"dataset1/data1/what","undetect",&fundetect,1);
	/* prodpar attribute for specific products */
        if(ODIM_quantcode==ODIM_ETOP) 
	   H5LTset_attribute_double(OUTH5,"dataset1/data1/what","prodpar",&TOPthr,1);
        if(ODIM_quantcode==ODIM_THEIGHT) 
	   H5LTset_attribute_double(OUTH5,"dataset1/data1/what","prodpar",&isotherm,1);
        if(ODIM_quantcode==ODIM_PHEIGHT) 
	   H5LTset_attribute_double(OUTH5,"dataset1/data1/what","prodpar",&isobar,1);

        if(obstime[0])
	{
	  sprintf(odate,"%.8s",obstime);
	  sprintf(otime,"%.4s00",obstime+8);
	}

        H5LTset_attribute_string(OUTH5,"what","date",odate);
        H5LTset_attribute_string(OUTH5,"what","time",otime);
        H5LTset_attribute_string(OUTH5,"what","object",ODIMobject);
        H5LTset_attribute_string(OUTH5,"what","source",ODIMsource);
        H5LTset_attribute_string(OUTH5,"what","version","H5rad 2.1");

        if(endtimes[0][0] || Acchours)
	{
          if(Acchours) sprintf(endtimes[0],"%s",obstime);
	  endsecs=sec_from_date(endtimes[0]);
          endepochs=(double)endsecs;
          if(Acchours)
	  {
             startsecs=endsecs-3600*Acchours;
             startepochs=(double)startsecs;
             date_from_sec(starttimes[0],startsecs);
	  }
	  sprintf(odate,"%.8s",endtimes[0]);
	  sprintf(otime,"%.4s00",endtimes[0]+8);
          H5LTset_attribute_string(OUTH5,"dataset1/what","enddate",odate);
          H5LTset_attribute_string(OUTH5,"dataset1/what","endtime",otime);
          H5LTset_attribute_double(OUTH5,"dataset1/how","endepochs",&endepochs,1);
	}

        if(starttimes[0][0])
	{
          if(!startsecs)
	  {
             startsecs=sec_from_date(starttimes[0]);
             startepochs=(double)startsecs;
	  }
	  sprintf(odate,"%.8s",starttimes[0]);
	  sprintf(otime,"%.4s00",starttimes[0]+8);
          H5LTset_attribute_string(OUTH5,"dataset1/what","startdate",odate);
          H5LTset_attribute_string(OUTH5,"dataset1/what","starttime",otime);
          H5LTset_attribute_double(OUTH5,"dataset1/how","startepochs",&startepochs,1);
	}

        if(ODIMproduct[0])
	{
          H5LTset_attribute_string(OUTH5,"dataset1/what","product",ODIMproduct);
	}

        if(ODIMACCnum) H5LTset_attribute_long(OUTH5,"dataset1/how","ACCnum",&ODIMACCnum,1);
    
        for(Co=SW;Co<=SE;Co++)
        {
           DEG[Co].u = coords[GEO][OUT][Co].u/DEG_TO_RAD;
           DEG[Co].v = coords[GEO][OUT][Co].v/DEG_TO_RAD;
         }

        H5LTset_attribute_double(OUTH5,"where","LL_lon",&DEG[SW].u,1);
        H5LTset_attribute_double(OUTH5,"where","LL_lat",&DEG[SW].v,1);
        H5LTset_attribute_double(OUTH5,"where","UL_lon",&DEG[NW].u,1);
        H5LTset_attribute_double(OUTH5,"where","UL_lat",&DEG[NW].v,1);
        H5LTset_attribute_double(OUTH5,"where","UR_lon",&DEG[NE].u,1);
        H5LTset_attribute_double(OUTH5,"where","UR_lat",&DEG[NE].v,1);
        H5LTset_attribute_double(OUTH5,"where","LR_lon",&DEG[SE].u,1);
        H5LTset_attribute_double(OUTH5,"where","LR_lat",&DEG[SE].v,1);

        xsize=(int64_t)outdim[XD];
        ysize=(int64_t)outdim[YD];
#if __SIZEOF_LONG__==4
        H5LTset_attribute_long_long(OUTH5,"where","xsize",&xsize,1);
        H5LTset_attribute_long_long(OUTH5,"where","ysize",&ysize,1);
#else
        H5LTset_attribute_long(OUTH5,"where","xsize",&xsize,1);
        H5LTset_attribute_long(OUTH5,"where","ysize",&ysize,1);
#endif
        H5LTset_attribute_double(OUTH5,"where","xscale",&outres[XD],1);
        H5LTset_attribute_double(OUTH5,"where","yscale",&outres[YD],1);

        H5LTset_attribute_string(OUTH5,"where","projdef",projdef[OUT]);

        if(ODIM_NODES) H5LTset_attribute_string(OUTH5,"how","nodes",nodes);
        else H5LTset_attribute_string(OUTH5,"how","nodes","all");
	if(ODIM_quantcode == ODIM_ACRR || ODIM_quantcode == ODIM_RATE)
	{
           H5LTset_attribute_double(OUTH5,"how","zr_a",&ZRA,1);
           H5LTset_attribute_double(OUTH5,"how","zr_b",&ZRB,1);
           if(conv_Rlim>0)
	   {
              H5LTset_attribute_double(OUTH5,"how","zr_a_conv",&ZRAc,1);
              H5LTset_attribute_double(OUTH5,"how","zr_b_conv",&ZRBc,1);
              H5LTset_attribute_double(OUTH5,"how","conv_limit_dBZ",&conv_dBZlim,1);
              H5LTset_attribute_double(OUTH5,"how","conv_limit_rate",&conv_Rlim,1);
	   }
           if(ZSA>0)
	   {
              H5LTset_attribute_double(OUTH5,"how","zs_a",&ZSA,1);
              H5LTset_attribute_double(OUTH5,"how","zs_b",&ZSB,1);
	   }
	}

        H5Dclose(dataset);
        H5Gclose(whatg);
        H5Gclose(whereg);
        H5Gclose(howg);
        H5Gclose(datawhatg);
        H5Gclose(datag);
        H5Gclose(datasetg);
	
        H5Fclose(OUTH5);
  }


  free(datarr);
  if(!LUT) fwrite(LUTarr,lutsize,1,LUTF);
  fclose(LUTF);
  free(LUTarr);
  free(inarr);
  return(RetVal);
}

/*######################################### F U N C T I O N S ###################################################*/

int configreader(char *cfgfile)
{

  FILE *CFGFILE = NULL;
  char cfgline[1000],optname[1000],optarg[1000];
  char *p;
  int namelen,ret;
  struct option option;

  CFGFILE=fopen(cfgfile,"r");
  if(!CFGFILE) return(0);
  while(1)
  {
    memset(cfgline,0,1000);
    retp=fgets(cfgline,999,CFGFILE);
    if(feof(CFGFILE)) break;

    if(cfgline[0]=='#') continue;
    p=strchr(cfgline,'=');
    if(p) namelen=p-cfgline; else namelen=strlen(cfgline)-1;
    memset(optname,0,1000);
    strncpy(optname,cfgline,namelen);
    option.name=optname;
    if(p) sprintf(optarg,"%s",p+1); else optarg[0]=0;
    if(VERB) printf("Config ");
    ret=option_solver(option,optarg);
    if(!ret) if(!QUIET) printf("Unknown config option %s\n",cfgline);
  }
  fclose(CFGFILE);
  return(1);
}

/*____________________________________________________________________________________________________________*/


int option_handler(int argc, char **argv)  /* option handler */
{
               static char short_options[500]={0};
               static struct option option; 

               int c,oi=0,li=0,options=0,old_optind=1;
               int option_index = 0;

               opterr=0;
               while(1)
               {
                  if(long_options[li].name==NULL) break;
                  short_options[oi++]=long_options[li].val; 
                  if(long_options[li++].has_arg)  short_options[oi++]=':';
               }
               options=li;

               while(1)
               { 
		  option_index=-1;
                  c = getopt_long (argc, argv, short_options,long_options, &option_index);
                  if (c == -1) break;
                  if(option_index==-1 && !QUIET) printf("Unknown option %s\n",argv[old_optind]); 
                  option=long_options[option_index];
                  option.val=c;
                  old_optind=optind;

                  if(option_index>=0)
		  {
                        if(VERB) printf("Option handler: #%d --%s=%s\n",option_index,option.name,optarg);
                        if(VERB) printf("Option ");
		  } else continue;
                  if(option_solver(option,optarg)) continue;

                  switch(c)
                  {
                     case '?': case 'h':
		        print_options();
                     default:
                         if(!QUIET) printf ("Unknown option --%s.\n", option.name);
                  }
               }
               return(0);
}

/*____________________________________________________________________________________________________________*/

void print_options(void)
{
   int oi;

   printf("\n\nConverts radar data projections, geometry and format.\n");
   printf("All options should be given in long format, e.g. --option=value or without value (flag options).\n");
   printf("Accepted input formats are IRIS, ODIM HDF5, PGM and PPM.\n");
   printf("Accepted output formats are GeoTIFF, ODIM HDF5, PPM and PGM.\n");
   printf("If verbose mode is wanted, give --verbose as the first option.\n");
   printf("\nAll options are possible to give also in the configuration file (--cfgfile=file),\n");
   printf("or as environment variables. In both cases use option name without leading -- .\n");
   printf("Handling order of options are: configuration file, environment variables, command line options.\n");

   printf("\nCurrently available ODIM quantities when using --ODIMquantity option for GeoTIFF output are:\n\n");
   
   for(oi=0;;oi++)
   {
     if(QuantStr[oi].ODIM_qstr[0]==0) break;
     printf("%-15s%s\n",QuantStr[oi].ODIM_qstr,QuantStr[oi].GTIFF_qstr);
   }

   printf("\n\nAvailable options:\n\n");

   for(oi=0;oi<num_options;oi++)
   {
      printf("--%-15s\t%s\n",optionlist[oi].name,optionlist[oi].expl);
   }
   printf("\n\nUsage: reprojection_radardata [--option=value | --flag_option] inputfile [outputfiles]\n\n");
   /*  printf("\n\nGive the input file as the first non-option argument, the rest of arguments are reserved for output files.\n"); */
   printf("\n\n");

#  ifdef NO_GEOTIFF
   printf("!!! The current executable is not compiled with libgdal, so GeoTIFF output is not possible.\nCompile without -DNO_GEOTIFF if GeoTIFF output is needed !!!\n\n");
#  endif   
}

/*____________________________________________________________________________________________________________*/


char *toDMS(double degs)
{
  int D,M;
  double S;

  D=(int)degs;
  S=3600.0*(degs-D);
  M=S/60;
  S=(S-60*M);
  if(D<0) { M=-M; S=-S; }
  sprintf(DMS,"%4dd %2d\'%7.4f\"",D,M,S);

  return(DMS);
}

/*____________________________________________________________________________________________________________*/

/* Iteration to get IRIS geographical SW and NE corners of IRIS area using cartesian 
   northing/easting as initial values for distance along spherical great circle between 
   area center and center of north/east border of IRIS area. 
   Gnomonic projection uses the actual ellipsoid for geographical/cartesian
   conversions. All other IRIS projections uses sphere of 6371 km radius. */

void give_IRIS_geocorners(projUV *IRIScorners, double reflon, double reflat,
     double a, double f, char *projtype, projUV geoC, double *IRISrange)
{
    int Co;
    double Sa=6371000.0;
    double cartoffs[2]={0.0},gd,dS;
    PJ *projref,*spheref,*effref;
    char projdef[500],sphedef[500];
    projUV cartC,cartT,geoT,cartNE,cartSW;

    sprintf(projdef,"+proj=%s +lon_0=%.8f +lat_0=%.8f +a=%.6f +f=%.12f",
            projtype,reflon,reflat,a,f);
    sprintf(sphedef,"+proj=%s +lon_0=%.8f +lat_0=%.8f +a=%.6f +f=0",
            projtype,reflon,reflat,Sa);
    if(VERB) printf("IRIS projection definition %s\nCenter lon = %.8f, lat =  %.8f\n",projdef,geoC.u/DEG_TO_RAD,geoC.v/DEG_TO_RAD); 

    projref=pj_init_plus(projdef);    
    spheref=pj_init_plus(sphedef);

    cartC=pj_fwd(geoC,projref);
    if(strcmp(projtype,"gnom")==0) effref=projref; else effref=spheref;

    for(Co=0;Co<=1;Co++)
    {
       cartT=cartC;
       Co ? (cartT.v+=IRISrange[YD]) : (cartT.u+=IRISrange[XD]); 

       while(1)
       {
          geoT=pj_inv(cartT,effref);
          gd = geodist(geoC,geoT,Sa,0);
          dS = gd-IRISrange[Co];
	  /*          printf("%d %.8f %.8f, %.8f %.8f\n",Co,gd,dS,cartT.u-cartC.u,cartT.v-cartC.v); */
          if(fabs(dS) <= 0.1) break;
          Co ? (cartT.v-=dS) : (cartT.u-=dS);
       }
       Co ? (cartoffs[Co] = cartT.v-cartC.v) : (cartoffs[Co] = cartT.u-cartC.u);
    }
   
    cartSW.u = cartC.u - cartoffs[XD];
    cartSW.v = cartC.v - cartoffs[YD];
    cartNE.u = cartC.u + cartoffs[XD];
    cartNE.v = cartC.v + cartoffs[YD];

    IRIScorners[0]=pj_inv(cartSW,projref);
    IRIScorners[1]=pj_inv(cartNE,projref);
 
    if(VERB)
    {
      printf("IRIS SW corner %.8f %.8f\n",IRIScorners[0].u/DEG_TO_RAD,IRIScorners[0].v/DEG_TO_RAD);
      printf("IRIS NE corner %.8f %.8f\n",IRIScorners[1].u/DEG_TO_RAD,IRIScorners[1].v/DEG_TO_RAD);
    }

    return;
}  

/*____________________________________________________________________________________________________________*/

/* Geodesic distance between two points of geographic coordinates in given ellipsoid of a-radius (a) and flattening (f)
   (Laberts's approximation formulae, error < 5 m / 50 degree) */

double geodist(projUV geo1, projUV geo2,double a,double f)
{
  double lon1=geo1.u,lat1=geo1.v,lon2=geo2.u,lat2=geo2.v;
  double psi1,psi2,X,Y,P,Q,d,r1,sigma,sinlats,sinlons;
  double r,sinP,sinQ,cosP,cosQ,sinsigma,sinhsig,coshsig;

  if(f<1e-8) d = a*acos(cos(lat1)*cos(lat2)*cos(lon1 - lon2)+sin(lat1)*sin(lat2)); /* Sphere*/
  else /* Ellipsoid with Lambert's approximation */
  {
     r=1.0/f;
     r1=(r-1)/r;

     psi1 = atan(r1 * tan(lat1)); 
     psi2 = atan(r1 * tan(lat2)); 
 
     sinlats = sin((psi2 - psi1)/2.0);
     sinlons = sin((lon2 - lon1)/2.0);

     sigma = 2.0 * asin(sqrt(sinlats*sinlats + cos(psi1)*cos(psi2)*sinlons*sinlons)); 

     P = (psi1 + psi2)/2.0;  
     Q = (psi2 - psi2)/2.0;  

     sinP = sin(P);
     sinQ = sin(Q);
     cosP = cos(P);
     cosQ = cos(Q);

     sinsigma = sin(sigma);
     sinhsig = sin(sigma/2.0);
     coshsig = sin(sigma/2.0);

     X = (sigma - sinsigma)*sinP*sinP*cosQ*cosQ/coshsig*coshsig;  
     Y = (sigma + sinsigma)*cosP*cosP*sinQ*sinQ/sinhsig*sinhsig;  

     d = a*(sigma - 0.5*(X+Y)/r);
  }
  return(d);
}

/* _________________________________  dBZ to R _________________________________________________________________________ */

double dBZtoR(double dBZ)
{
  double RR;

  if(dBZ<conv_dBZlim) RR=pow(10.0,dBZ/ZRB10 + ZRC);
  else RR=pow(10.0,dBZ/ZRB10c + ZRCc);

  return(RR);
}

/*___________________________________get CSV item _________________________________________________________________________*/

double get_csv_item(char *csv, int it)
{
   char *sp,*ep,sval[100]={0};
   int i;
   double val;

   sp=csv;
   for(i=1;i<it;i++) sp=strchr(sp,',')+1;
   ep=strchr(sp,',');
   if(ep==NULL) sprintf(sval,"%s",sp);
   else strncpy(sval,sp,ep-sp);
   val=atof(sval);

   return(val);
}
/*--------------------------------------------------------------------------------------------------*/

uint16_t IrisRain(double F)
{
  uint32_t E,Eb=0,D;
  uint16_t R=0;
  
  while(1)
  {
    D=(uint32_t)(1000.0*F + 1e-6);

    if(D <= 8192)      { R=D;     break; }
    if(D >= 134184960) { R=65534; break; }

    E = D >> 12;
    while(E)
    {
       E >>= 1;
       Eb++;
    }
    R = (Eb << 12) | ((D >> (Eb-1)) & 0xFFF); 

    break;
  }

  return(R);
}

/*--------------------------------------------------------------------------------------------------*/

void date_from_sec(char *date,time_t secs)
{
   struct tm *Sdd;

   Sdd=gmtime(&secs);

   sprintf(date,"%d%02d%02d%02d%02d",Sdd->tm_year+1900,Sdd->tm_mon+1,
                                     Sdd->tm_mday,Sdd->tm_hour,
                                     Sdd->tm_min);
   return;
}

/*--------------------------------------------------------------------------------------------------*/

time_t sec_from_date(char *date)
{
   struct tm Sdd;
   time_t secs;

   sscanf(date,"%4d%2d%2d%2d%2d",&Sdd.tm_year,&Sdd.tm_mon,
                                 &Sdd.tm_mday,&Sdd.tm_hour,
                                 &Sdd.tm_min);
   Sdd.tm_year-=1900;
   Sdd.tm_mon--;
   Sdd.tm_sec=0;
   secs=mktime(&Sdd);
   return(secs);
}


/*____________________________________________________________________________________________________________*/

int option_solver(struct option option, char *optarg)
{
          int c,len;      

          c=option.val;

          if(optarg)
          {
             len=strlen(optarg);
             if(optarg[len-1]=='\n') optarg[len-1]=0;
          }          
          if(VERB) printf("solver  --%s",option.name);
          if(optarg && VERB) printf("=%s\n",optarg); else if(VERB) printf("\n");

          if(c=='h' || c=='?') return(0); 

          if(!strcmp(option.name,"cfgfile"))
            {
               strcpy(cfgfile,optarg);
               /* set options according to possible config file. If --cfgfile option in command line
               preceeds other options, those options overrides config file options! */
               if(!configreader(cfgfile)) printf("Config file %s don't exist!\n",cfgfile);
               return(1);
            }

          if(!strcmp(option.name,"LUTdir"))
            {
  	       sprintf(LUTdir,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"inSWlon"))
            {
               coords[GEO][IN][SW].u = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"inSWlat"))
            {
               coords[GEO][IN][SW].v = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"inNElon"))
            {
               coords[GEO][IN][NE].u = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"inNElat"))
            {
               coords[GEO][IN][NE].v = atof(optarg)*DEG_TO_RAD;
               return(1);
            }
          if(!strcmp(option.name,"inCElon"))
            {
               coords[GEO][IN][CE].u = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"inCElat"))
            {
               coords[GEO][IN][CE].v = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"inBBOXgeo"))
            {
	       coords[GEO][IN][SW].u = get_csv_item(optarg,1)*DEG_TO_RAD;
	       coords[GEO][IN][SW].v = get_csv_item(optarg,2)*DEG_TO_RAD;
	       coords[GEO][IN][NE].u = get_csv_item(optarg,3)*DEG_TO_RAD;
	       coords[GEO][IN][NE].v = get_csv_item(optarg,4)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"indim"))
            {
               indim[XD] = indim[YD] = atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"inXdim"))
            {
               indim[XD] = atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"inYdim"))
            {
               indim[YD] = atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"inres"))
            {
               inres[XD] = inres[YD] = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"inXres"))
            {
               inres[XD] = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"inYres"))
            {
               inres[YD] = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"inrange"))
            {
               inrange[XD] = inrange[YD] = atof(optarg)*1000.0;
               return(1);
            }

          if(!strcmp(option.name,"inXrange"))
            {
               inrange[XD] = atof(optarg)*1000.0;
               return(1);
            }
          if(!strcmp(option.name,"inYrange"))
            {
               inrange[YD] = atof(optarg)*1000.0;
               return(1);
            }

          if(!strcmp(option.name,"inproj"))
            {
               sprintf(projdef[IN],"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"infile"))
            {
               sprintf(infile,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"outSWlon"))
            {
               coords[GEO][OUT][SW].u = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"outSWlat"))
            {
               coords[GEO][OUT][SW].v = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"outNElon"))
            {
               coords[GEO][OUT][NE].u = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"outNElat"))
            {
               coords[GEO][OUT][NE].v = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"outCElon"))
            {
               coords[GEO][OUT][CE].u = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"outCElat"))
            {
               coords[GEO][OUT][CE].v = atof(optarg)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"outBBOXgeo"))
            {
	       coords[GEO][OUT][SW].u = get_csv_item(optarg,1)*DEG_TO_RAD;
	       coords[GEO][OUT][SW].v = get_csv_item(optarg,2)*DEG_TO_RAD;
	       coords[GEO][OUT][NE].u = get_csv_item(optarg,3)*DEG_TO_RAD;
	       coords[GEO][OUT][NE].v = get_csv_item(optarg,4)*DEG_TO_RAD;
               return(1);
            }

          if(!strcmp(option.name,"outBBOXcart"))
            {
	       coords[CART][OUT][SW].u = get_csv_item(optarg,1);
	       coords[CART][OUT][SW].v = get_csv_item(optarg,2);
	       coords[CART][OUT][NE].u = get_csv_item(optarg,3);
	       coords[CART][OUT][NE].v = get_csv_item(optarg,4);
               OUT_BBOX_CART=1;
               return(1);
            }

          if(!strcmp(option.name,"lon_0"))
            {
               lon_0 = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"lat_0"))
            {
               lat_0 = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"outdim"))
            {
               outdim[XD] = outdim[YD] = atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"outXdim"))
            {
               outdim[XD] = atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"outYdim"))
            {
               outdim[YD] = atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"outres"))
            {
               outres[XD] = outres[YD] = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"outXres"))
            {
               outres[XD] = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"outYres"))
            {
               outres[YD] = atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"outrange"))
            {
               outrange[XD] = outrange[YD] = atof(optarg)*1000.0;
               return(1);
            }

          if(!strcmp(option.name,"outXrange"))
            {
               outrange[XD] = atof(optarg)*1000.0;
               return(1);
            }

          if(!strcmp(option.name,"outYrange"))
            {
               outrange[YD] = atof(optarg)*1000.0;
               return(1);
            }

          if(!strcmp(option.name,"outproj"))
            {
               sprintf(projdef[OUT],"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"outformat"))
            {
              sprintf(outformat,"%s",optarg);
              return(1);
            }

          if(!strcmp(option.name,"outfile"))
            {
               sprintf(outfile,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"outbytes"))
            {
               OUTBYTES=atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ODIMsource"))
            {
	       sprintf(ODIMsource,"%s",optarg);;
               return(1);
            }

          if(!strcmp(option.name,"ODIMquantity"))
            {
	       int q;

               sprintf(ODIMquantity,"%s",optarg);
               if(strcmp(ODIMquantity,"ACCR")==0) sprintf(ODIMquantity,"ACRR");
               /* ODIM string search here */
             
               for(q=0;;q++)
	       {
                 if(QuantStr[q].ODIM_qstr[0]==0) return(0);
		 if(strcmp(ODIMquantity,QuantStr[q].ODIM_qstr)==0)
		 {
		   ODIM_quantcode=q;
                   break;
		 }
	       }
               return(1);
            }

          if(!strcmp(option.name,"sweepnumber"))
            {
	       sweepnumber=atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ODIMdatasetnum"))
            {
	       datasetnumber=atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ODIMdatanum"))
            {
	       datanumber=atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ODIMobject"))
            {
	       sprintf(ODIMobject,"%s",optarg);;
               return(1);
            }

          if(!strcmp(option.name,"ODIMproduct"))
            {
	       sprintf(ODIMproduct,"%s",optarg);;
               return(1);
            }

          if(!strcmp(option.name,"ODIMACCnum"))
            {
	       ODIMACCnum=atol(optarg);
               return(1);
            }

          if(!strcmp(option.name,"polarzoom"))
            {
	       polarzoom=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"polarrange"))
            {
	       polarrange=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"use_incorners"))
            {
               USE_INSW=1;
               USE_INNE=1;
               return(1);
            }

          if(!strcmp(option.name,"EPSG"))
            {
 	       EPSG=atoi(optarg);
               return(1);
            }

          if(!strcmp(option.name,"TIFFcompression"))
            {
 	       sprintf(TIFFcompression,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"Zlevel"))
            {
	       Zlevel=atoi(optarg);
	       if(Zlevel>0 && Zlevel<10) sprintf(sZlevel,"%s",optarg); else Zlevel=0;
               return(1);
            }

          if(!strcmp(option.name,"tiling"))
            {
	       tiling=atoi(optarg);
               if(tiling>1) { sprintf(Xtile,"%d",tiling); sprintf(Ytile,"%d",tiling); }
               return(1);
            }

          if(!strcmp(option.name,"Xtile"))
            {
	       sprintf(Xtile,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"Ytile"))
            {
	       sprintf(Ytile,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"overviews"))
            {
	       int i,list[20];
               char *sp,*ep,sval[50];
           
               sp=ep=optarg;
               for(i=0;ep;i++) 
 	       {
		  memset(sval,0,50);
                  ep=strchr(sp,',');
                  if(ep==NULL) sprintf(sval,"%s",sp);
                  else strncpy(sval,sp,ep-sp);
                  if(i==0) sprintf(overview_method,"%s",sval);
                  else list[i-1]=atoi(sval);
                  sp=ep+1;
               }
               overviews=i-1;
               if(!overviews) eff_OverviewList=NULL;
               else
	       {
                  OverviewList=calloc(overviews,sizeof(int));
                  for(i=0;i<overviews;i++) OverviewList[i]=list[i];
                  eff_OverviewList=OverviewList;
	       }
               return(1);
            }

          if(!strcmp(option.name,"obstime"))
            {
 	       sprintf(obstime,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"starttimes"))
            {
	       int i,len;
               char *p=NULL,*timep=optarg;

               for(i=0;;i++)
	       {
                 p=strchr(timep,',');
 	         if(p) len=p-timep; else len=strlen(timep); 
                 strncpy(starttimes[i],timep,len);
                 if(!p) break;
                 timep=p+1;
	       }

               return(1);
            }

          if(!strcmp(option.name,"endtimes"))
            {
	       int i,len;
               char *p=NULL,*timep=optarg;

               for(i=0;;i++)
	       {
                 p=strchr(timep,',');
 	         if(p) len=p-timep; else len=strlen(timep); 
                 strncpy(endtimes[i],timep,len);
                 if(!p) break;
                 timep=p+1;
	       }

               return(1);
            }

          if(!strcmp(option.name,"nodefaultmeta"))
            {
               DEFAULTMETA=0;     
               return(1);
            }

          if(!strcmp(option.name,"metafile"))
            {
 	       sprintf(metafile,"%s",optarg);
               INSERT_METAFILE=1;     
               return(1);
            }

          if(!strcmp(option.name,"QDparam"))
            {
 	       sprintf(QDparam,"%s",optarg);
               return(1);
            }

          if(!strcmp(option.name,"use_inSW"))
            {
               USE_INSW=1;
               return(1);
            }

          if(!strcmp(option.name,"use_inNE"))
            {
               USE_INNE=1;
               return(1);
            }

          if(!strcmp(option.name,"use_incenter"))
            {
               USE_INCE=1;
               return(1);
            }

          if(!strcmp(option.name,"ignore_incenter"))
            {
               IGN_INCE=1;
               return(1);
            }

          if(!strcmp(option.name,"QDheaders"))
            {
               QDHEADERS=1;
               return(1);
            }

          if(!strcmp(option.name,"RtodBZ"))
            {
               RATE_TO_DBZ=1;
               return(1);
            }

          if(!strcmp(option.name,"dBZtoR"))
            {
	       DBZ_TO_RATE=1;
               return(1);
            }

          if(!strcmp(option.name,"ZRA"))
            {
	       ZRA=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ZRB"))
            {
	       ZRB=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ZSA"))
            {
	       ZSA=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ZSB"))
            {
	       ZSB=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ZRAc"))
            {
	       ZRAc=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"ZRBc"))
            {
	       ZRBc=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"Rfactor"))
            {
	       Rfactor=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"isotherm"))
            {
	       isotherm=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"isobar"))
            {
	       isobar=atof(optarg);
               return(1);
            }

          if(!strcmp(option.name,"flip_output"))
            {
	       FLIP_OUTPUT_DATA=1;
               return(1);
            }

          if(!strcmp(option.name,"fill_nodata"))
            {
	       FILL_NODATA=1;
               return(1);
            }

          if(!strcmp(option.name,"overwrite"))
            {
	       OVERWRITE=1;
               return(1);
            }

          if(!strcmp(option.name,"IRISrain"))
            {
	       IRISRAIN=1;
               return(1);
            }

          if(!strcmp(option.name,"verbose"))
            {
               VERB=1;
               return(1);
            }

          if(!strcmp(option.name,"quiet"))
            {
               QUIET=1;
               return(1);
            }

          return(0);
}  

