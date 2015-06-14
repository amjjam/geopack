/******************************************************************************
 * This is geopack.H - it defines C++ interface functions to most, but not    *
 * Fortran subroutines in the Fortran geopack library.                        *
 ******************************************************************************/

#ifndef _GEOPACK_H_
#define _GEOPACK_H_

#include <iostream>
#include <vector>
#include <math.h>

extern "C"{
  void igrf_gsw_08_(double *xGsw, double *yGsw, double *zGsw,
		    double *hxgsw,double *hygsw,double *hzgsw);
  void igrf_geo_08_(double *r, double *theta, double *phi, 
		    double *br, double *btheta, double *bphi);  
  void dip_08_(double *xgsw, double *ygsx, double *zgsw,
	       double *bxgsw, double *bygsw, double *bzgsw);
  void sun_08_(int *yr, int *dy, int *hr, int *mn, int *se,
	       double *gst, double *slong, double *srasn, double *sdec);
  void sphcar_08_(double *r, double *theta, double *phi,
		  double *x, double *y, double *z, int *j);
  void bspcar_08_(double *theta, double *phi, double *br, double *btheta, 
		  double *bphi,double *x, double *y, double *z);
  void bcarsp_08_(double *x, double *y, double *z, double *bx, double *by,
		  double *bz, double *br, double *btheta, double *bphi);
  void recalc_08_(int *yr, int *dy, int *hr, int *mn, int *se,
		  double *vgsex, double *vgsey, double *vgsez);
  void gswgse_08_(double *xgsw, double *ygsw, double *zgsw,
		  double *xgse, double *ygse, double *zgse, int *j);
  void geomag_08_(double *xgsm, double *ygsm, double *zgsm,
		  double *xmag, double *ymag, double *zmag, int *j);
  void geigeo_08_(double *xgei, double *ygei, double *zgei,
		  double *xgeo, double *ygeo, double *zgeo, int *j);
  void magsm_08_(double *xmag, double *ymag, double *zmag,
		 double *xsm, double *ysm, double *zsm, int *j);
  void smgsw_08_(double *xsm, double *ysm, double *zsm,
		 double *xgsw, double *ygsw, double *zgsw, int *j);
  void geogsw_(double *xgeo, double *ygeo, double *zgeo,
	       double *xgsw, double *ygsw, double *zgsw, int *j);
  void geodgeo_(double *h, double *xmum, double *r, double *theta, int *j);
}

void igrfGSW(double xGsw, double yGsw, double zGsw, 
	     double &hxgsw, double &gygsw, double &hzgsw);
void igrfGEO(double r, double theta, double phi, 
	     double &br, double &btheta, double &bphi);
void dip(double xgsw, double ygsw, double zgsw,
	 double &bxgsw, double &bygsw, double &bzgsw);
void sun(int yr, int dy, int hr, int mn, int se, double &gst, double &slong,
	 double &srasn, double &sdec);
void sphcar(double r, double theta, double phi, 
	    double &x, double &y, double &z);
void sphcar(double x, double y, double z,
	    double &r, double &theta, double &phi); 
void bspcar(double theta, double phi, double br, double btheta, double bphi,
	    double &bx, double &by, double &bz);
void bcarsp(double x, double y, double z, double bx, double by, double bz,
	    double &br, double &btheta, double &bphi);
void recalc(int yr, int dy, int hr, int mn, int se, 
	    double vxgse, double vygse, double vzgse);
void gswgse(double xgsw, double ygsw, double zgsw,
	    double &xgse, double &ygse, double &zgse);
void gsegsw(double xgse, double ygse, double zgse,
	    double &xgsw, double &ygsw, double &zgsw);
void geomag(double xgeo, double ygeo, double zgeo, 
	    double &xmag, double &ymag, double &zmag);
void maggeo(double xmag, double ymag, double zmag,
	    double &xgeo, double &ygeo, double &zgeo);
void geigeo(double xgei, double ygei, double zgei,
	    double &xgeo, double &ygeo, double &zgeo);
void geogei(double xgeo, double ygeo, double zgeo,
	    double &xgei, double &ygei, double &zgei);
void magsm(double xmag, double ymag, double zmag,
	   double &xsm, double &ysm, double &zsm);
void smmag(double xsm, double ysm, double zsm,
	   double &xmag, double &ymag, double &zmag);
void smgsw(double xsm, double ysm, double zsm,
	   double &xgsw, double &ygsw, double &zgsw);
void gswsm(double xgsw, double ygsw, double zgsw,
	   double &xsm, double &ysm, double &zsm);
void geogsw(double xgeo, double ygeo, double zgeo,
	    double &xgsw, double &ygsw, double &zgsw);
void gswgeo(double xgsw, double ygsw, double zgsw,
	    double &xgeo, double &ygeo, double &zgeo);
void geodgeo(double h, double xmu, double &r, double &theta);
void geogeod(double r, double theta, double &h, double &xmu);
void rhand(double x, double y, double z, double ds3, 
	   double &r1, double &r2, double &r3,
	   void (*field)(double x, double y, double z, 
			 double &bx, double &by, double &bz));
void step(double &x, double &y, double &z, double &ds, double dsmax, 
	  double errin, 
	  void (*field)(double x, double y, double z, 
			double &bx, double &by, double &bz));
void trace(double xi, double yi, double zi, double dir, double dsmax, 
	   double err, double rlim, double r0, 
	   void (*field)(double x, double y, double z,
			 double &bx, double &by, double &bz),
	   double &xf, double &yf, double &zf,
	   std::vector<double> &xx, std::vector<double> &yy, 
	   std::vector<double> &zz, int &l, int lmax);

// rhand_08
// step_08
// trace_08
// shuetal_mgnp_08
// t96_mgnp_08

#endif
