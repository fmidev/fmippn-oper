# include <hdf5.h>
# include <hdf5_hl.h>
# include <stdlib.h>

int main(int argc, char *argv[])
{
  hid_t srcfile, dstfile,fcplist;
  herr_t ret;

  H5Eset_auto(H5E_DEFAULT,NULL,NULL); /*suppress HDF error messages */
  srcfile=H5Fopen(argv[1],H5F_ACC_RDONLY,H5P_DEFAULT);
  fcplist=H5Fget_create_plist(srcfile); /* store creation parameter list of source file */

  dstfile=H5Fcreate(argv[2],H5F_ACC_TRUNC,fcplist,H5P_DEFAULT);
  H5Ocopy(srcfile,"meta",dstfile,"meta",H5P_DEFAULT,H5P_DEFAULT);
  H5Ocopy(srcfile,"motion",dstfile,"motion",H5P_DEFAULT,H5P_DEFAULT);
  H5Fclose(dstfile);
  H5Fclose(srcfile);

  return(0);
}

