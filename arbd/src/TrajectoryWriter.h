///////////////////////////////////////////////////////////////////////  
// Author: Jeff Comer <jcomer2@illinois.edu>
#ifndef TRAJECTORYWRITER_H
#define TRAJECTORYWRITER_H

#define PDB_TEMPLATE_LINE "ATOM      1  CA  MOL S   1      -6.210  -9.711   3.288  0.00  0.00      ION"

#include <cstdio>
#include "useful.h"
#include "DcdWriter.h"

class TrajectoryWriter {
public:
  const static int formatDcd = 0;
  const static int formatPdb = 1;
  const static int formatTraj = 2;

  TrajectoryWriter(const char* filePrefix, const char* formatName, Matrix3 box0, int num0, float timestep0, int outputPeriod0) 
    : box(box0), num(num0), timestep(timestep0), outputPeriod(outputPeriod0)  {
    pdbTemplate = PDB_TEMPLATE_LINE;
    format = getFormatCode(String(formatName));
    makeUnitCell();

    fileName = filePrefix;
    fileName.add(".");
    fileName.add(getFormatName(format));
    
    if (format == formatDcd) {
      dcd = new DcdWriter(fileName);
      dcd->writeHeader(fileName.val(), num, 1, 0, outputPeriod, 0, timestep, 1);
    }
  }

  ~TrajectoryWriter() {
    if (format == formatDcd) delete dcd;
  }

private:
  Matrix3 box;
  String fileName;
  int format;
  String pdbTemplate;
  double unitCell[6];		/* use double for dcd format */
  int num;
  float timestep;
  int outputPeriod;
  DcdWriter* dcd;

  void makeUnitCell() {
    float pi = 4.0f*atan(1.0f);

    unitCell[0] = box.ex().length();
    unitCell[2] = box.ey().length();
    unitCell[5] = box.ez().length();
    
    float bc = box.ey().dot(box.ez());
    float ac = box.ex().dot(box.ez());
    float ab = box.ex().dot(box.ey());

    unitCell[1] = bc/unitCell[0]/unitCell[2]/pi*180.0f;
    unitCell[3] = ac/unitCell[0]/unitCell[5]/pi*180.0f;
    unitCell[4] = ab/unitCell[0]/unitCell[1]/pi*180.0f;
  }

public:
  static int getFormatCode(String format) {
    format.lower();
    if (format == String("dcd")) return formatDcd;
    if (format == String("pdb")) return formatPdb;
    if (format == String("traj")) return formatTraj;
    return formatDcd;
  }

  static String getFormatName(int formatCode) {
    switch(formatCode) {
    case formatPdb:
      return String("pdb");
    case formatTraj:
      return String("traj");
    case formatDcd:
      return String("dcd");
    default:
      return String("dcd");
    }
  }

  void newFile(const Vector3* pos, const String* name, float t, int n) const {
    switch(format) {
    case formatPdb:
      newPdb(fileName, pos, name);
      break;
    case formatTraj:
      newTraj(pos, name, t, n);
      break;
    case formatDcd:
    default:
      newDcd(pos, name);
      break;
    }
  }

  void newFile(const Vector3* pos, const String* name, const int* id, float t, int n) const {
    switch(format) {
    case formatPdb:
      newPdb(fileName, pos, name);
      break;
    case formatTraj:
      newTraj(pos, name, id, t, n);
      break;
    case formatDcd:
    default:
      newDcd(pos, name);
      break;
    }
  }

  void append(const Vector3* pos, const String* name, float t, int n) const {
    switch(format) {
    case formatPdb:
      appendPdb(pos, name);
      break;
    case formatTraj:
      appendTraj(pos, name, t, n);
      break;
    case formatDcd:
    default:
      appendDcd(pos);
      break;
    }
  }
  
  void append(const Vector3* pos, const String* name, const int* id, float t, int n) const {
    switch(format) {
    case formatPdb:
      appendPdb(pos, name);
      break;
    case formatTraj:
      appendTraj(pos, name, id, t, n);
      break;
    case formatDcd:
    default:
      appendDcd(pos);
      break;
    }
  }

  void newPdb(const char* outFile, const Vector3* pos, const String* name) const {
    char s[128];

    sprintf(s, "CRYST1   %.3f   %.3f   %.3f  90.00  90.00  90.00 P 1           1\n", box.exx, box.eyy, box.ezz);
    String sysLine(s);

    sprintf(s, "REMARK   frameTime %.10g ns\n", outputPeriod*timestep);
    String remarkLine(s);
    
    String line;

    FILE* out = fopen(outFile, "w");
    fprintf(out, "%s", sysLine.val());
    fprintf(out, "%s", remarkLine.val());

    for (int i = 0; i < num; i++) {
      line = makePdbLine(pdbTemplate, i, name[i], i, name[i], pos[i], 0.0);
      fprintf(out, "%s",  line.val());
      fprintf(out, "\n");
    }
    fprintf(out, "END\n");
    fclose(out);
  }

  void appendPdb(const Vector3* pos, const String* name) const {
    String line;

    FILE* out = fopen(fileName, "a");
    for (int i = 0; i < num; i++) {
      line = makePdbLine(pdbTemplate, i, name[i], i, name[i], pos[i], 0.0);
      fprintf(out, "%s", line.val());
      fprintf(out, "\n");
    }
    fprintf(out, "END\n");
    fclose(out);
  }
 
  void newTraj(const Vector3* pos, const String* name, float t, int n) const {
    FILE* out = fopen(fileName, "w");
    for (int i = 0; i < n; i++)
      fprintf(out, "%s %.10g %.10g %.10g %.10g\n", name[i].val(), t, pos[i].x, pos[i].y, pos[i].z);
    fprintf(out, "END\n");
    fclose(out);
  }

  void newTraj(const Vector3* pos, const String* name, const int* id, float t, int n) const {
    FILE* out = fopen(fileName, "w");
    for (int i = 0; i < n; i++)
      fprintf(out, "%s %.10g %.10g %.10g %.10g %d\n", name[i].val(), t, pos[i].x, pos[i].y, pos[i].z, id[i]);
    fprintf(out, "END\n");
    fclose(out);
  }

  void appendTraj(const Vector3* pos, const String* name, float t, int n) const {
    FILE* out = fopen(fileName, "a");
    for (int i = 0; i < n; i++)
      fprintf(out, "%s %.10g %.10g %.10g %.10g\n", name[i].val(), t, pos[i].x, pos[i].y, pos[i].z);
    fprintf(out, "END\n");
    fclose(out);
  }

  void appendTraj(const Vector3* pos, const String* name, const int* id, float t, int n) const {
    FILE* out = fopen(fileName, "a");
    for (int i = 0; i < n; i++)
      fprintf(out, "%s %.10g %.10g %.10g %.10g %d\n", name[i].val(), t, pos[i].x, pos[i].y, pos[i].z, id[i]);
    fprintf(out, "END\n");
    fclose(out);
  }

  void newDcd(const Vector3* pos, const String* name) const {
    /*  
    // Write a new pdb to store the atom names and such.
    char pdbFile[128];
    sprintf(pdbFile, "%s.pdb", fileName.val());
    newPdb(pdbFile, pos, name);
    */

    // Write first frame.
    appendDcd(pos);
  }

  void appendDcd(const Vector3* pos) const {
    float* x = new float[num];
    float* y = new float[num];
    float* z = new float[num];

    for (int i = 0; i < num; i++) {
      x[i] = pos[i].x;
      y[i] = pos[i].y;
      z[i] = pos[i].z;
    }
    dcd->writeStep(num, x, y, z, unitCell);

    delete[] x;
    delete[] y;
    delete[] z;
  }

  static String makePdbLine(const String& tempLine, int index, const String& segName, int resId, 
			    const String& name, Vector3 r, float beta) {
    char s[128];

    String record("ATOM  ");
    sprintf(s, "     %5i ", index);
    String si = String(s).range(-6,-1);
    if (name.length() == 4) sprintf(s, "%s   ", name.val());
    else sprintf(s, " %s   ", name.val());
    String nam = String(s).range(0,3);
    String temp0 = tempLine.range(16,21);
  
    sprintf(s, "    %d", resId);
    String res = String(s).range(-4,-1);
    String temp1 = tempLine.range(26,29);
  
    sprintf(s,"       %.3f", r.x);
    String sx = String(s).range(-8,-1);
    sprintf(s,"       %.3f", r.y);
    String sy = String(s).range(-8,-1);
    sprintf(s,"       %.3f", r.z);
    String sz = String(s).range(-8,-1);

    String temp2 = tempLine.range(54,59);
    sprintf(s,"    %.2f", beta);
    String bet = String(s).range(-6,-1);
    String temp3 = tempLine.range(66,71);

    sprintf(s, "%s    ", segName.val());
    String seg = String(s).range(0,3);

    String ret(record);
    ret.add(si);
    ret.add(nam);
    ret.add(temp0);
    ret.add(res);
    ret.add(temp1);
    ret.add(sx);
    ret.add(sy);
    ret.add(sz);
    ret.add(temp2);
    ret.add(bet);
    ret.add(temp3);
    ret.add(seg);
  
    return ret;
  }
};
#endif
