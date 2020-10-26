#!/bin/bash
# Configuration settings of fmippn
# if PROCTIME is exported, the input data check is not performed, and TIMESTAMP is set to PROCTIME
# Using this script as source terminates the calling script if exit is called (like with option --list) !

if [ $# == 0 ] && [ ! $DOMAIN ]; then
   HELP=True
else
   OPTS=`echo $* | tr ' ' "\n" | grep ^'\-\-'`
   HELP=`echo $OPTS | grep help`
   LIST=`echo $OPTS | grep list`
   POSTPROC=`echo $OPTS | grep postproc`
   ARGDOMAIN=`echo $* | tr ' ' "\n" | grep -v '\-\-' | head -1`
fi

# echo $OPTS 
# echo $LIST

if [ $ARGDOMAIN ]; then export DOMAIN=${ARGDOMAIN} ; fi
if [ ! $DOMAIN ]; then export DOMAIN=test ; fi
export DOMAIN=${DOMAIN:-'test'}

if [ $HELP ]; then
   echo -e "\nUsage: set_common_conf.sh [--help --list] DOMAIN\n"
   echo -e "This script sould be executed as source in the beginning of each script running PPN"
   echo -e "with argument DOMAIN, e.g. \"source $HOME/ppn/config/set_common_conf.sh europe\""
   echo -e "The DOMAIN argument should be used when sourcing this script," 
   echo -e "otherwise the PPN will run with DOMAIN set as environment variable"
   echo -e "or \"test\" if not specified anywhere. In this case the PPN will run using data from test directory\n"
   echo -e "The current domain is \"${DOMAIN}\"\n"
   echo -e "With option --list, the environment and configuration files of specified domain will be listed.\n"
   exit 1
fi

export TZ=UTC
export LOGDAY=`date -u +%d`
export LOGHOUR=`date -u +%H`
export RECONFIG=66

export CONDAENV=${CONDAENV:-'fmippn'}
export DOMAIN=${DOMAIN:-'test'}
export OPERDIR=${OPERDIR:-"$HOME/fmippn-oper"}
export CONFDIR=$OPERDIR/config
export PPNDIR=$OPERDIR/fmippn
export DOMAINDIR=$PPNDIR/config
export RUNDIR=$OPERDIR/run-and-distribution
export PREPROCDIR=$OPERDIR/preprocess
export POSTPROCDIR=$OPERDIR/postprocess
export LOGDIR=/var/tmp/ppn/log
 if [ ! -e $LOGDIR ]; then mkdir -p $LOGDIR ; fi

# define common functions
source $CONFDIR/common_functions.sh

export RUNLOG=$LOGDIR/run_DOMAIN=${DOMAIN}_day=${LOGDAY}.log
export PPNLOG=$LOGDIR/ppn_DOMAIN=${DOMAIN}_hour=${LOGHOUR}.log
if [ ! $LIST ]; then
   ln -sf $RUNLOG $LOGDIR/run_DOMAIN=${DOMAIN}.log
   ln -sf $PPNLOG $LOGDIR/ppn_DOMAIN=${DOMAIN}.log
fi
export POSTPROCLOG=$LOGDIR/postproc_DOMAIN=${DOMAIN}.log
export PREPROCLOG=$LOGDIR/preproc_DOMAIN=${DOMAIN}.log
export DISTRIBLOG=$LOGDIR/distrib_DOMAIN=${DOMAIN}.log

export LD_LIBRARY_PATH=`cat $CONFDIR/ld_library_path`
export PYTHON=$HOME/miniconda3/envs/fmippn/bin/python
export START_CONDA=$HOME/miniconda3/etc/profile.d/conda.sh

# pysteps configuration file (data input of all domains are configured here)
export PYSTEPSRC=$PPNDIR/pystepsrc

# Domain specific configuration
export DOMAINCONF=$DOMAINDIR/${DOMAIN}.json
if [ ! -e $DOMAINCONF ]; then
   echo "Create $DOMAINCONF for domain \"$DOMAIN\" first!"
   exit 2
fi

# Default output directories
export PPN_OUTPUT_DIR=`grep OUTPUT_PATH $DOMAINCONF | cut -d: -f2 | tr -d \"\,' '`
export INTERPDIR=$PPN_OUTPUT_DIR/interp
export OBSDIR=$PPN_OUTPUT_DIR/obs # Contains the latest observed radar field
export ENSMEANDIR=$PPN_OUTPUT_DIR/ensmean
export PROBDIR=$PPN_OUTPUT_DIR/prob
export VISUALDIR=$PPN_OUTPUT_DIR/visual
export DISTRIBDIR=$PPN_OUTPUT_DIR/distrib

if [ ! -d $PPN_OUTPUT_DIR ]; then mkdir -p $PPN_OUTPUT_DIR ; fi
if [ ! -d $OBSDIR ]; then mkdir -p $OBSDIR ; fi
if [ ! -d $INTERPDIR ]; then mkdir -p $INTERPDIR ; fi
if [ ! -d $ENSMEANDIR ]; then mkdir -p $ENSMEANDIR ; fi
if [ ! -d $PROBDIR ]; then mkdir -p $PROBDIR ; fi
if [ ! -d $VISUALDIR ]; then mkdir -p $VISUALDIR ; fi
if [ ! -d $DISTRIBDIR ]; then mkdir -p $DISTRIBDIR ; fi

# Set possible environment variables for postprocessing of the domain data
export POSTPROC_CONF=$CONFDIR/postproc_DOMAIN=${DOMAIN}.conf
if [ -e $POSTPROC_CONF ]; then source $POSTPROC_CONF ; fi


# Input data from pysteps config
DATAROOT=`grep -A100 \"$DOMAIN\": $PYSTEPSRC | grep root_path | head -1 | cut -d: -f2 | tr -d \"\,' '`
if [ ! $DATAROOT ]; then
   echo "No input data defined for domain \"$DOMAIN\" in $PYSTEPSRC"
   exit 3
fi

export INPUT_DATAROOT=`eval echo $DATAROOT`
export INPUT_DATADIR=`grep -A100 \"$DOMAIN\": $PYSTEPSRC | grep path_fmt | head -1 | cut -d: -f2 | tr -d \"\,' '`
EXT=`grep -A100 \"$DOMAIN\": $PYSTEPSRC | grep fn_ext | head -1 | cut -d: -f2 | tr -d \"\,' '`
export FILEPATTERN='*_'`grep -A100 \"$DOMAIN\": $PYSTEPSRC | grep fn_pattern | head -1 | cut -d: -f2 | tr -d \"\, | cut -d_ -f2-`.${EXT}
export INPUT_PATH=${INPUT_DATAROOT}/${INPUT_DATADIR}

if [ $PROCTIME ]; then
   export TIMESTAMP=$PROCTIME
   NAMETAIL=`echo $FILEPATTERN | cut -d'*' -f2`
   export LATEST_OBSFILE=$OBSDIR/"$TIMESTAMP""$NAMETAIL"
   echo $LATEST_OBSFILE
   if [ ! -e $LATEST_OBSFILE ]; then
      cp $INPUT_PATH/`basename $LATEST_OBSFILE` $LATEST_OBSFILE
   fi
else
   # Last observed radar data (time attributes of the input files must be in real time order!)
   # This is used as argument to PPN model run, and as first field of interpolation
   export OBSFILE=`ls -1t ${INPUT_PATH}/${FILEPATTERN} | head -1`
   if [ ! $OBSFILE ]; then
      echo "No input files found!"
      exit 4
   fi

   # Store the latest obs file to temporary obs directory
   cp $OBSFILE $OBSDIR
   OBSNAME=`basename $OBSFILE`
   export LATEST_OBSFILE=$OBSDIR/$OBSNAME

   # Setting the default timestamp requires the input filename starting with YYYYmmddHHMM
   # If not so, please preprocess the filenames, otherwise parse the timestamp from OBSFILE
   export TIMESTAMP=`basename $OBSFILE | cut -c-12`
   # ADD HERE VALIDITY CHECk OF THE TIMESTAMP !
fi

# CPU cores per member for parallel execution in pysteps
if [ -e $CONFDIR/membercores_DOMAIN=${DOMAIN}.conf ]; then
   OMPTHREADS=`cat $CONFDIR/membercores_DOMAIN=${DOMAIN}.conf`
fi
if [ ! $OMP_NUM_THREADS ]; then
   export OMP_NUM_THREADS=${OMPTHREADS:-"1"}
fi

# Finally, the output HDF5 file of PPN model run
export PPN_OUTPUT_FILE=$PPN_OUTPUT_DIR/nc_${TIMESTAMP}.h5

# List the environment and configurations
if [ $LIST ]; then
       echo '_________________________________'
       echo -e "Environment\n"
       echo '_________________________________'
       echo CONDAENV=${CONDAENV}
       echo DOMAIN=${DOMAIN}
       echo OPERDIR=${OPERDIR}
       echo CONFDIR=${CONFDIR}
       echo PPNDIR=${PPNDIR}
       echo DOMAINDIR=${DOMAINDIR}
       echo RUNDIR=${RUNDIR}
       echo PREPROCDIR=${PREPROCDIR}
       echo POSTPROCDIR=${POSTPROCDIR}
       echo DOMAINCONF=${DOMAINCONF}
       echo POSTPROC_CONF=${POSTPROC_CONF}
       ls -1 ${POSTPROC_CONF}
       echo LOGDIR=${LOGDIR}
       echo PPNLOG=${PPNLOG}
       echo POSTPROCLOG=${POSTPROCLOG}
       echo PREPROCLOG=${PREPROCLOG}
       echo DISTRIBLOG=${DISTRIBLOG}

       echo PYTHON=${PYTHON}
       echo START_CONDA=${START_CONDA}

       echo PYSTEPSRC=${PYSTEPSRC}

       echo PPN_OUTPUT_DIR=${PPN_OUTPUT_DIR}
       echo OBSDIR=${OBSDIR}
       echo LATEST_OBSFILE=${LATEST_OBSFILE}
       echo INTERPDIR=${INTERPDIR}
       echo ENSMEANDIR=${ENSMEANDIR}
       echo PROBDIR=${PROBDIR}
       echo VISUALDIR=${VISUALDIR}
       echo DISTRIBDIR=${DISTRIBDIR}
       echo INPUT_DATAROOT=${INPUT_DATAROOT}
       echo INPUT_DATADIR=${INPUT_DATADIR}
       echo INPUT_PATH=${INPUT_PATH}
       echo FILEPATTERN=${FILEPATTERN}
       echo OBSFILE=${OBSFILE}
       echo TIMESTAMP=${TIMESTAMP}
       echo PPN_OUTPUT_FILE=${PPN_OUTPUT_FILE}
       echo -n "CPU cores per member for parallel execution in pySTEPS: "
       echo OMP_NUM_THREADS=${OMP_NUM_THREADS}

       echo '_________________________________'
       echo "Configuration files"
       echo '_________________________________'
       echo -e "\n"PYSTEPSRC=${PYSTEPSRC}", configuration of \"${DOMAIN}\" data source \n"
       cat $PYSTEPSRC | grep  -A100 -B1 \"${DOMAIN}\": | grep -m2 -B100 \/\/ | head --lines=-1
       echo '_________________________________'
       echo -e "\n"DOMAINCONF=${DOMAINCONF}" :\n"    
       cat $DOMAINCONF

       if [ -e $POSTPROC_CONF ]; then
       echo '_________________________________'
          echo "Postprocessing environment for domain \"$DOMAIN\""
          echo -e "\n"POSTPROC_CONF=${POSTPROC_CONF}" :\n"
          grep -v ^# $POSTPROC_CONF | grep -v ^$ | cut -d' ' -f2-
       fi
       echo '_________________________________'

       exit 1
fi

# Starting anaconda/miniconda
source $START_CONDA
if [ $? != 0 ]; then
   echo "No anaconda available!"  
   exit 5
fi

# Activate fmippn in conda
conda activate $CONDAENV

