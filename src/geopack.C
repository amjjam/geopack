#include "../include/geopack.H"

void igrfGSW(double xGsw, double yGsw, double zGsw, 
	     double &hxgsw, double &hygsw, double &hzgsw){
  igrf_gsw_08_(&xGsw,&yGsw,&zGsw,&hxgsw,&hygsw,&hzgsw);
}

void igrfGEO(double r, double theta, double phi, 
	     double &br, double &btheta, double &bphi){
  igrf_geo_08_(&r,&theta,&phi,&br,&btheta,&bphi);
}

void dip(double xgsw, double ygsw, double zgsw,
	 double &bxgsw, double &bygsw, double &bzgsw){
  dip_08_(&xgsw,&ygsw,&zgsw,&bxgsw,&bygsw,&bzgsw);
}

void sun(int yr, int dy, int hr, int mn, int se, double &gst, double &slong,
	 double &srasn, double &sdec){
  sun_08_(&yr,&dy,&hr,&mn,&se,&gst,&slong,&srasn,&sdec);
}

void sphcar(double r, double theta, double phi, 
	    double &x, double &y, double &z){
  int j=1;
  sphcar_08_(&r,&theta,&phi,&x,&y,&z,&j);
}

void carsph(double x, double y, double z,
	    double &r, double &theta, double &phi){
  int j=-1;
  sphcar_08_(&r,&theta,&phi,&x,&y,&z,&j);
}

void bspcar(double theta, double phi, double br, double btheta, double bphi,
	    double &bx, double &by, double &bz){
  bspcar_08_(&theta,&phi,&br,&btheta,&bphi,&bx,&by,&bz);
}

void bcarsp(double x, double y, double z, double bx, double by, double bz,
	    double &br, double &btheta, double &bphi){
  bcarsp_08_(&x,&y,&z,&bx,&by,&bz,&br,&btheta,&bphi);
}

void recalc(int yr, int dy, int hr, int mn, int se, 
	    double vxgse, double vygse, double vzgse){
  recalc_08_(&yr,&dy,&hr,&mn,&se,&vxgse,&vygse,&vzgse);
}

void gswgse(double xgsw, double ygsw, double zgsw,
	    double &xgse, double &ygse, double &zgse){
  int j=1;
  gswgse_08_(&xgsw,&ygsw,&zgsw,&xgse,&ygse,&zgse,&j);
}

void gsegsw(double xgse, double ygse, double zgse,
	    double &xgsw, double &ygsw, double &zgsw){
  int j=-1;
  gswgse_08_(&xgsw,&ygsw,&zgsw,&xgse,&ygse,&zgse,&j);
}

void geomag(double xgeo, double ygeo, double zgeo, 
	    double &xmag, double &ymag, double &zmag){
  int j=1;
  geomag_08_(&xgeo,&ygeo,&zgeo,&xmag,&ymag,&zmag,&j);
}

void maggeo(double xmag, double ymag, double zmag,
	    double &xgeo, double &ygeo, double &zgeo){
  int j=-1;
  geomag_08_(&xgeo,&ygeo,&zgeo,&xmag,&ymag,&zmag,&j);
}

void geigeo(double xgei, double ygei, double zgei,
	    double &xgeo, double &ygeo, double &zgeo){
  int j=1;
  geigeo_08_(&xgei,&ygei,&zgei,&xgeo,&ygeo,&zgeo,&j);
}

void geogei(double xgeo, double ygeo, double zgeo,
	    double &xgei, double &ygei, double &zgei){
  int j=1;
  geigeo_08_(&xgei,&ygei,&zgei,&xgeo,&ygeo,&zgeo,&j);
}

void magsm(double xmag, double ymag, double zmag,
	   double &xsm, double &ysm, double &zsm){
  int j=1;
  magsm_08_(&xmag,&ymag,&zmag,&xsm,&ysm,&zsm,&j);
}

void smmag(double xsm, double ysm, double zsm,
	   double &xmag, double &ymag, double &zmag){
  int j=-1;
  magsm_08_(&xmag,&ymag,&zmag,&xsm,&ysm,&zsm,&j);
}

void smgsw(double xsm, double ysm, double zsm,
	   double &xgsw, double &ygsw, double &zgsw){
  int j=1;
  smgsw_08_(&xsm,&ysm,&zsm,&xgsw,&ygsw,&zgsw,&j);
}

void gswsm(double xgsw, double ygsw, double zgsw,
	   double &xsm, double &ysm, double &zsm){
  int j=-1;
  smgsw_08_(&xsm,&ysm,&zsm,&xgsw,&ygsw,&zgsw,&j);
}

void geogsw(double xgeo, double ygeo, double zgeo,
	    double &xgsw, double &ygsw, double &zgsw){
  int j=1;
  geogsw_08_(&xgeo,&ygeo,&zgeo,&xgsw,&ygsw,&zgsw,&j);
}

void gswgeo(double xgsw, double ygsw, double zgsw,
	    double &xgeo, double &ygeo, double &zgeo){
  int j=-1;
  geogsw_08_(&xgeo,&ygeo,&zgeo,&xgsw,&ygsw,&zgsw,&j);
}

void geodgeo(double h, double xmu, double &r, double &theta){
  int j=1;
  geodgeo_08_(&h,&xmu,&r,&theta,&j);
}

void geogeod(double r, double theta, double &h, double &xmu){
  int j=-1;
  geodgeo_08_(&h,&xmu,&r,&theta,&j);
}


/*=============================================================================
  void rhand(double x, double y, double z, double ds3, 
             double &r1, double &r2, double &r3,
             void (*field)(double x, double y, double z, 
                           double &bx, double &by, double &bz)) -
  calculates the right-hand side of the field tracing equation. It is
  a little different from the Fortran routine in that ds3 must be
  passed in and in that only one field model is passed int. It is up
  to the user to define a field model which incorporates the internal
  and external fields, and in arranging to pass appropriate parameters
  tot he external field model, if necessary.  
  ============================================================================*/
void rhand(double x, double y, double z, double ds3, 
	   double &r1, double &r2, double &r3,
	   void (*field)(double x, double y, double z, 
			 double &bx, double &by, double &bz)){
  double bx,by,bz,b;
  field(x,y,z,bx,by,bz);
  b=ds3/sqrt(bx*bx+by*by+bz*bz);
  r1=bx*b;
  r2=by*b;
  r3=bz*b;
}


/*=============================================================================
  void step(double &x, double &y, double &z, double &ds, double dsmax, 
            double errin, 
            void (*field)(double x, double y, double z, 
	                  double &bx, double &by, double &bz)) - step update

  This function takes a step along the field line. The position x,y,z
  is updated to move ds along the field line (or a smaller step if
  needed, in which case ds is updated to the smaller size), keeping
  the estimated error below some value.
  
  double x,y,z - input start position, output end position
  doubld ds - input step length, output adjusted step length
  double dsmax - largest step size permitted.
  double erriin - largest position error permitted
  ============================================================================*/
void step(double &x, double &y, double &z, double &ds, double dsmax, 
	  double errin, 
	  void (*field)(double x, double y, double z, 
			double &bx, double &by, double &bz)){

  // Don't make any step larger than size dsmax
  if(fabs(ds)>dsmax){
    if(ds>0)
      ds=fabs(dsmax);
    else
      ds=-fabs(dsmax);
  }

  // Loop until precision is good enough
  double r11,r12,r13;
  double r41,r42,r43;
  double r51,r52,r53;
  double errcur;
  bool done=false;
  do{
    double ds3=-ds/3;
    rhand(x,y,z,ds3,r11,r12,r13,field);
    double r21,r22,r23;
    rhand(x+r11,y+r12,z+r13,ds3,r21,r22,r23,field);
    double r31,r32,r33;
    rhand(x+0.5*(r11+r21),y+0.5*(r12+r22),z+0.5*(r13+r23),ds3,
	  r31,r32,r33,field);
    rhand(x+0.375*(r11+3*r31),y+0.375*(r12+3*r32),z+0.375*(r13+3*r33),ds3,
	  r41,r42,r43,field);
    rhand(x+1.5*(r11-3*r31+4*r41),y+1.5*(r12-3*r32+4*r42),
	  z+1.5*(r13-3*r33+4*r43),ds3,r51,r52,r53,field);
    if((errcur=fabs(r11-4.5*r31+4*r41-0.5*r51)+fabs(r12-4.5*r32+4*r42-0.5*r52)+
	fabs(r13-4.5*r33+4*r43-0.5*r53))>errin)
      ds*=0.5;
    else
      done=true;
  }while(!done);

  x+=0.5*(r11+4*r41+r51);
  y+=0.5*(r12+4*r42+r52);
  z+=0.5*(r13+4*r43+r53);

  // If the error is really small then increase the step size
  if(errcur<0.04*errin&&ds<dsmax/1.5)
    ds*=1.5;
}

  
/*=============================================================================
  void trace(double xi, double yi, double zi, double dsmax, double err, 
	   double rlim, double r0, 
	   void (*field)(double x, double y, double z,
			 double &bx, double &by, double &bz),
	   double &xf, double &yf, double &zf,
	   std::vector<double> &xx, std::vector<double> &yy, 
	   std::vector<double> &zz, int l, int lmax) - trace a field line
  and fill vectors with the field line positions
  double xi,yi,zi - starting position
  double dir - the direction of field line tracing. If dir=1.0 then the 
    tracing direction is made antiparallel to the total field vector (e.g. 
    from Northern to Southern conjugate poin). If dir=-1.0 then the tracing 
    proceeds in the opposite direction, that is, parallel to the total field 
    vector. 
  double dsmax - largest step size permitted. This is the largest desired 
    spacing between points along the field line.
  double err - largest error permitted. According to the original Fortran 
    code on which this is based a value of err=0.0001 is sufficient for 
    most applications. 
  double rlim - radius of a spere in Earth radii which defines the outer edge 
    of the tracing volume. Field line tracing stops if the field line 
    traces outside of this radius.
  double r0 - radius of spere in Earth radii which defines the inner 
    boundary of the tracing volume. Field line tracing stops if the 
    field line traces inside this radius. 
  field - the magnetic field function to trace. 
  std::vector<double> xx, yy, zz - vectors in which the traced points 
    along the field line will be returned. 
  int &l - number of field line points traced.
  int lmax - the maximum number of field line points to trace. 

  NOTE: 

  The trace function will resize the vectors xx, yy, zz to lmax before
  starting, just in case they were not already set, and it will resize
  them to l before returning.

  presumable if l==lmax the tracing stopped because it ran out of
  points, or, by coincidence the number of points given, lmax, was
  exactly right.
  ============================================================================*/
void trace(double xi, double yi, double zi, double dir, double dsmax, 
	   double err, double rlim, double r0, 
	   void (*field)(double x, double y, double z,
			 double &bx, double &by, double &bz),
	   double &xf, double &yf, double &zf,
	   std::vector<double> &xx, std::vector<double> &yy, 
	   std::vector<double> &zz, int &l, int lmax){
  
  xx.resize(lmax);
  yy.resize(lmax);
  zz.resize(lmax);

  l=0;
  int nrev=0;
  double ds=0.5*dir;
  double x=xi,y=yi,z=zi,r1,r2,r3,ad;

  rhand(x,y,z,dir,r1,r2,r3,field);
  if(x*r1+y*r2+z*r3<0)
    ad=-0.01;
  else 
    ad=0.01;

  double rr=sqrt(x*x+y*y+z*z)+ad;
  double ryz,r,xr,yr,zr;
  double drp,dr;

  for(l=0;l<lmax;l++){
    xx[l]=x;
    yy[l]=y;
    zz[l]=z;
    ryz=y*y+z*z;
    r2=x*x+ryz;
    r=sqrt(r2);
    
    if(r>rlim||ryz>1600||x>20)
      break;

    if(r<r0||rr>r){
      r1=(r0-r)/(rr-r);
      x=x-(x-xr)*r1;
      y=y-(y-yr)*r1;
      z=z-(z-zr)*r1;
      break;
    }
    
    if(!(r>=rr||r>=3)){
      double fc=0.2;
      if(r-r0<0.05)
	fc=0.05;
      double al=fc*(r-r0+0.2);
      ds=dir*al;
    }
    
    xr=x;
    yr=y;
    zr=z;
    
    drp=r-rr;
    rr=r;
    
    step(x,y,z,ds,dsmax,err,field);
    
    r=sqrt(x*x+y*y+z*z);
    dr=r-rr;
    if(drp*dr<0)
      nrev++;
    if(nrev>4)
      break;
  }
  if(l==lmax)
    std::cout << "trace: loop ended early because l reached lmax" << std::endl;
  
  xf=x;
  yf=y;
  zf=z;
  
  l--;
  xx[l]=xf;
  yy[l]=yf;
  zz[l]=zf;

  l++;
  xx.resize(l);
  yy.resize(l);
  zz.resize(l);
}
