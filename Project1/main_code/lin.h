#ifndef LIN_H
#define LIN_H

class Lin
{
private:
  double x,y,z ;
  double norm ;
  double * pos [3] ;
  double m_gamma ;
  double m_gamma_sq ;
public:
  Lin () ;
  Lin (double x0, double y0, double z0) ;
  Lin (double , double , double , double ) ;
  Lin& operator= (const Lin& c) ;
  ~Lin () {}
  double X () const ;
  double Y () const ;
  double Z () const ;
  double GAMMA () const ;
  double GAMMA_SQ () const ;
  double abs () const ;
  double R_sq () const ;
  double R_sq_gamma () const ;
  double R_sq_gamma_sq () const ;
  double get ( int ) const ;
  Lin unit () ;
  friend Lin operator+ (const Lin&  a, const Lin& b) ;
  friend Lin operator- (const Lin&  a, const Lin& b) ;
  friend double operator, (const Lin&  a, const Lin& b) ;
  friend Lin operator* (const double&  a, const Lin& b) ;
  friend Lin operator/ (const Lin&  r, const Lin& v) ;
  void read () const ;
  void tweak ( int , double ) ;
  void tweak_star ( int , double ) ;
  void set ( int , double ) ;
  void set_gamma ( double ) ;
};

inline Lin operator+ ( const Lin& a , const Lin& b )
{ return Lin (a.X()+b.X(), a.Y()+b.Y(), a.Z()+b.Z()) ; }

inline Lin operator- ( const Lin& a , const Lin& b )
{ return Lin (a.X()-b.X(), a.Y()-b.Y(), a.Z()-b.Z()) ; }

inline double operator, ( const Lin& a , const Lin& b )
{ return a.X()*b.X() + a.Y()*b.Y() + a.Z()*b.Z() ; }

inline Lin operator* ( const double& a , const Lin& b )
{ return Lin (a*b.X(), a*b.Y(), a*b.Z()) ; }

inline Lin operator*= ( const Lin& r , const Lin& v )
{ return Lin(r.Y()*v.Z()-r.Z()*v.Y(),-r.X()*v.Z()+r.Z()*v.X(),r.X()*v.Y()-r.Y()*v.X()) ; }

#endif
