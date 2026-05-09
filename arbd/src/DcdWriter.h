///////////////////////////////////////////////////////////////////////  
// Modified dcd reader from NAMD.
// Author: Jeff Comer <jcomer2@illinois.edu>

/**
***  Copyright (c) 1995, 1996, 1997, 1998, 1999, 2000 by
***  The Board of Trustees of the University of Illinois.
***  All rights reserved.
**/

/*
  dcdlib contains C routines for reading and writing binary DCD
  files.  The output format of these files is based on binary FORTRAN
  output, so its pretty ugly.  If you are squeamish, don't look!
*/

#pragma once
#include <unistd.h>
#include <cstring>
#include <fstream>
#include <fcntl.h>
#include <sys/stat.h>

// using namespace std;

#define NFILE_POS (off_t) 8
#define NPRIV_POS (off_t) 12
#define NSAVC_POS (off_t) 16
#define NSTEP_POS (off_t) 20

#ifndef O_LARGEFILE
#define O_LARGEFILE 0x0
#endif

/*  DEFINE ERROR CODES THAT MAY BE RETURNED BY DCD ROUTINES		*/
#define DCD_DNE		-2	/*  DCD file does not exist		*/
#define DCD_OPENFAILED	-3	/*  Open of DCD file failed		*/
#define DCD_BADREAD 	-4	/*  read call on DCD file failed	*/
#define DCD_BADEOF	-5	/*  premature EOF found in DCD file	*/
#define DCD_BADFORMAT	-6	/*  format of DCD file is wrong		*/
#define DCD_FILEEXISTS  -7	/*  output file already exists		*/
#define DCD_BADMALLOC   -8	/*  malloc failed			*/

// Just use write instead of NAMD_write --JRC
#define NAMD_write write

class DcdWriter {
public:
  DcdWriter(const char* fileName) {
    fd = openDcd(fileName);    
    
    if (fd == DCD_OPENFAILED) {
      printf("DcdWriter::DcdWriter Failed to open dcd file %s.", fileName);
      exit(-1);
    }
  }

  ~DcdWriter() {
    closeDcd();
  }
private:
  int fd;

private:

void pad(char *s, int len)
{
	int curlen;
	int i;

	curlen=strlen(s);

	if (curlen>len)
	{
		s[len]='\0';
		return;
	}

	for (i=curlen; i<len; i++)
	{
		s[i]=' ';
	}

	s[i]='\0';
}


  /*********************************************************************/
  /*								     */
  /*			FUNCTION open_dcd_write			     */
  /*								     */
  /*   INPUTS:							     */
  /*	dcdfile - Name of the dcd file				     */
  /*								     */
  /*   OUTPUTS:							     */
  /*	returns an open file descriptor for writing		     */
  /*								     */
  /*	This function will open a dcd file for writing.  It takes    */
  /*   the filename to open as its only argument.	 It will return a    */
  /*   valid file descriptor if successful or DCD_OPENFAILED if the    */
  /*   open fails for some reason.  If the file specifed already       */
  /*   exists, it is renamed by appending .BAK to it.		     */
  /*								     */
  /*********************************************************************/
  int openDcd(const char* dcdname)
  {
    struct stat sbuf;
    int dcdfd;
    char *newdcdname = 0;

    if (stat(dcdname, &sbuf) == 0) 
      {
	newdcdname = new char[strlen(dcdname)+5];
	if(newdcdname == (char *) 0)
	  return DCD_OPENFAILED;
	strcpy(newdcdname, dcdname);
	strcat(newdcdname, ".BAK");
	if(rename(dcdname, newdcdname))
	  return(DCD_OPENFAILED);
	delete [] newdcdname;
      } 


    if ( (dcdfd = open(dcdname, O_RDWR|O_CREAT|O_EXCL|O_LARGEFILE,
		       S_IRUSR|S_IWUSR|S_IRGRP|S_IROTH)) < 0)
      {
	return(DCD_OPENFAILED);
      }

    return dcdfd;
  }

  /****************************************************************/
  /*								*/
  /*			FUNCTION close_dcd_write		*/
  /*								*/
  /*   INPUTS:							*/
  /*	fd - file descriptor to close				*/
  /*								*/
  /*   OUTPUTS:							*/
  /*	the file pointed to by fd				*/
  /*								*/
  /*	close_dcd_write close a dcd file that was opened for    */
  /*   writing							*/
  /*								*/
  /****************************************************************/

  void closeDcd()

  {	
    close(fd);
  }

public:
  /*****************************************************************************/
  /*									     */
  /*				FUNCTION write_dcdheader		     */
  /*									     */
  /*   INPUTS:								     */
  /*	fd - file descriptor for the dcd file				     */
  /*	filename - filename for output					     */
  /*	N - Number of atoms						     */
  /*	NFILE - Number of sets of coordinates				     */
  /*	NPRIV - Starting timestep of DCD file - NOT ZERO		     */
  /*	NSAVC - Timesteps between DCD saves				     */
  /*	NSTEP - Number of timesteps					     */
  /*	DELTA - length of a timestep					     */
  /*									     */
  /*   OUTPUTS:								     */
  /*	none								     */
  /*									     */
  /*	This function prints the "header" information to the DCD file.  Since*/
  /*   this is duplicating an unformatted binary output from FORTRAN, its ugly.*/
  /*   So if you're squeamish, don't look.					     */
  /*									     */
  /*****************************************************************************/
  int writeHeader(const char *filename, int N, int NFILE, int NPRIV, 
		  int NSAVC, int NSTEP, float DELTA, int with_unitcell)
  {
    int	out_integer;
    float   out_float;
    char	title_string[200];
    //int	user_id;
    time_t 	cur_time;
    struct  tm *tmbuf;
    char    time_str[11];

    out_integer = 84;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    strcpy(title_string, "CORD");
    NAMD_write(fd, title_string, 4);
    out_integer = NFILE;  /* located at fpos 8 */
    out_integer = 0;  /* ignore the lies */
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    out_integer = NPRIV;  /* located at fpos 12 */
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    out_integer = NSAVC;  /* located at fpos 16 */
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    out_integer = NSTEP;  /* located at fpos 20 */
    out_integer = NPRIV - NSAVC;  /* ignore the lies */
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    out_integer=0;
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    out_float = DELTA;
    NAMD_write(fd, (char *) &out_float, sizeof(float));
    out_integer = with_unitcell ? 1 : 0;
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    out_integer = 0;
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    out_integer = 24;  // PRETEND TO BE CHARMM24 -JCP
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    out_integer = 84;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));

    out_integer = 164;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    out_integer = 2;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));

    sprintf(title_string, "REMARKS FILENAME=%s CREATED BY NAMD", filename);
    pad(title_string, 80);
    NAMD_write(fd, title_string, 80);

    char username[100];
    //user_id= (int) getuid();
    //pwbuf=getpwuid(user_id);
    //if ( pwbuf ) sprintf(username,"%s",pwbuf->pw_name);
    //else sprintf(username,"%d",user_id);
    sprintf(username,"%s", "BrownTown");

    cur_time=time(NULL);
    tmbuf=localtime(&cur_time);
    strftime(time_str, 10, "%m/%d/%y", tmbuf);

    sprintf(title_string, "REMARKS DATE: %s CREATED BY USER: %s",
	    time_str, username);
    pad(title_string, 80);
    NAMD_write(fd, title_string, 80);
    out_integer = 164;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    out_integer = 4;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    out_integer = N;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));
    out_integer = 4;
    NAMD_write(fd, (char *) & out_integer, sizeof(int));

    return(0);
  }

  
  /************************************************************************/
  /*									*/
  /*				FUNCTION write_dcdstep			*/
  /*									*/
  /*   INPUTS:								*/
  /*	fd - file descriptor for the DCD file to write to		*/
  /*	N - Number of atoms						*/
  /*	X - X coordinates						*/
  /*	Y - Y coordinates						*/
  /*	Z - Z coordinates						*/
  /*  unitcell - a, b, c, alpha, beta, gamma of unit cell */
  /*									*/
  /*   OUTPUTS:								*/
  /*	none								*/
  /*									*/
  /*	write_dcdstep writes the coordinates out for a given timestep   */
  /*   to the specified DCD file.						*/
  /*                                                                      */
  /************************************************************************/
  int writeStep(int N, const float *X, const float *Y, const float *Z, const double *cell)

  {
    int NSAVC,NSTEP,NFILE;
    int out_integer;

    /* Unit cell */
    if (cell) {
      out_integer = 6*8;
      NAMD_write(fd, (char *) &out_integer, sizeof(int));
      NAMD_write(fd, (char *) cell, out_integer);
      NAMD_write(fd, (char *) &out_integer, sizeof(int));
    }

    /* Coordinates */
    out_integer = N*4;
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) X, out_integer);
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) Y, out_integer);
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) &out_integer, sizeof(int));
    NAMD_write(fd, (char *) Z, out_integer);
    NAMD_write(fd, (char *) &out_integer, sizeof(int));

    /* don't update header until after write succeeds */
    lseek(fd,NSAVC_POS,SEEK_SET);
    read(fd,(void*) &NSAVC,sizeof(int));
    lseek(fd,NSTEP_POS,SEEK_SET);
    read(fd,(void*) &NSTEP,sizeof(int));
    lseek(fd,NFILE_POS,SEEK_SET);
    read(fd,(void*) &NFILE,sizeof(int));
    NSTEP += NSAVC;
    NFILE += 1;
    lseek(fd,NSTEP_POS,SEEK_SET);
    NAMD_write(fd,(char*) &NSTEP,sizeof(int));
    lseek(fd,NFILE_POS,SEEK_SET);
    NAMD_write(fd,(char*) &NFILE,sizeof(int));
    lseek(fd,0,SEEK_END);

    return(0);
  }
};

