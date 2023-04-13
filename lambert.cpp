
#include <math.h>
#include <cmath>
#include <tgmath.h>
#include <string.h>
#include "lambert.h"
#include "seebatt.h"
#include "seebattk.h"
#include "vector.h"

void lambert(double ro[], double r[], char dm[], double Dtsec, double vo[], double v[]) {
	double small = 0.000001;
	double mu = 3.986004418e14;   // m3/s2
	double y1 = 0;
	double magr = norm(r, 3);
	double magro = norm(ro, 3);
	double CosDeltaNu = dot(ro,r, 3) / (magro*magr);
	double rcrossr[3];
	cross(ro,r,rcrossr);
	double magrcrossr = norm(rcrossr, 3);

	double SinDeltaNu;

	if (strcmp(dm,"pro") == 0)
	    SinDeltaNu = magrcrossr/(magro*magr);
	else
	    SinDeltaNu = - magrcrossr/(magro*magr);

	double DNu = atan2(SinDeltaNu, CosDeltaNu);

	// the angle needs to be positive to work for the long way
	if (DNu < 0.0)
	    DNu = 2.0 * M_PI + DNu;

	double RoR   = magr / magro;
	double eps   = RoR - 1.0;
	double tan2w = 0.25*eps*eps / ( sqrt( RoR ) + RoR * ( 2.0 + sqrt( RoR ) ) );
	double rp    = sqrt( magro*magr )*( (cos(DNu*0.25))*(cos(DNu*0.25)) + tan2w );
	
	double L;

	if ( DNu < M_PI )
		L = ( (sin(DNu*0.25))*(sin(DNu*0.25)) + tan2w ) / ( (sin(DNu*0.25))*(sin(DNu*0.25)) + tan2w + cos( DNu*0.5 ) );
	else
		L = ( (cos(DNu*0.25))*(cos(DNu*0.25)) + tan2w - cos( DNu*0.5 ) ) / ( (cos(DNu*0.25))*(cos(DNu*0.25)) + tan2w );

	double m     = mu*Dtsec*Dtsec / ( 8.0*rp*rp*rp );
	double x     = 10.0;
	double xn    = L;
	double chord = sqrt( magro*magro + magr*magr - 2.0*magro*magr*cos( DNu ) );
	double s     = ( magro + magr + chord )*0.5;
	double lim1  = sqrt(m/L);


	double tempx, Denom, h1, h2, b, u, k2, y;
	int Loops= 1;

	while (true) {
	    x     = xn;    
	    tempx = seebatt(x);
	    Denom = 1.0 / ( (1.0+2.0*x+L) * (4.0*x + tempx*(3.0+x) ) );
	    h1    = (L+x)*(L+x) * ( 1.0+ 3.0*x + tempx )*Denom;
	    h2    = m*( x - L + tempx )*Denom;
	    
	    // ----------------------- Evaluate CUBIC ------------------
	    b = 0.25*27.0*h2 / ((1.0+h1)*(1.0+h1)*(1.0+h1));
	    if (b < -1.0) // reset the initial condition
	        xn = 1.0 - 2.0*L;
	    else {
	        if (y1 > lim1)
	            xn = xn * (lim1 / y1);
	        else {
	            u = 0.5*b / ( 1.0 + sqrt( 1.0 + b ) );            
	            k2 = seebattk(u);
	            y = ( ( 1.0+h1 ) / 3.0 )*( 2.0 + sqrt( 1.0+b )/( 1.0+2.0e0*u*k2*k2 ) );
	            xn= sqrt( ((1.0e0-L)*0.5)*((1.0e0-L)*0.5) + m/(y*y) ) - ( 1.0e0+L )*0.5;
	        }
	    }
	    Loops ++;
	    y1 = sqrt(m/((L+x)*(1.0+x)) );
	    if ((abs(xn-x) < small) && (Loops > 30))
	        break;
	}

	double a = mu*Dtsec*Dtsec / (16.0*rp*rp*xn*y*y );

	double arg1, arg2, AlpH, BetH, DH, F, GDot, G, Sinv, Cosv, AlpE, BetE, am, ae, be, tm, DE;

	// ------------------ Find Eccentric anomalies -----------------
	// ------------------------ Hyperbolic -------------------------
	if ( a < -small ) {
	    arg1 = sqrt( s / ( -2.0*a ) );
	    arg2 = sqrt( ( s-chord ) / ( -2.0*a ) );
	    // ------- Evaluate f and g functions --------
	    AlpH = 2.0 * asinh( arg1 );
	    BetH = 2.0 * asinh( arg2 );
	    DH   = AlpH - BetH;
	    F    = 1.0 - (a/magro)*(1.0 - cosh(DH) );
	    GDot = 1.0 - (a/magr) *(1.0 - cosh(DH) );
	    G    = Dtsec - sqrt(-a*a*a/mu)*(sinh(DH)-DH);
	}
	else {
		// % ------------------------ Elliptical ---------------------
	    if ( a > small ) {
	        arg1 = sqrt( s / ( 2.0*a ) );
	        arg2 = sqrt( ( s-chord ) / ( 2.0*a ) );
	        Sinv = arg2;
	        Cosv = sqrt( 1.0 - (magro+magr-chord)/(4.0e0*a) );
	        BetE = 2.0*acos(Cosv);
	        BetE = 2.0*asin(Sinv);
	        if ( DNu > M_PI )
	            BetE= -BetE;

	        Cosv= sqrt( 1.0 - s/(2.0*a) );
	        Sinv= arg1;
	        am  = s*0.5;
	        ae  = M_PI;
	        be  = 2.0*asin( sqrt( (s-chord)/s ) );
	        tm  = sqrt(am*am*am/mu)*(ae - (be-sin(be)));
	        if ( Dtsec > tm )
	            AlpE= 2.0*M_PI-2.0*asin( Sinv );
	        else
	            AlpE= 2.0*asin( Sinv );

	        DE  = AlpE - BetE;
	        F   = 1.0 - (a/magro)*(1.0 - cos(DE) );
	        GDot= 1.0 - (a/magr)* (1.0 - cos(DE) );
	        G   = Dtsec - sqrt(a*a*a/mu)*(DE - sin(DE));
	    }
	    else {
	        // % --------------------- Parabolic ---------------------
	        arg1 = 0.0;
	        arg2 = 0.0;
			throw "a parabolic orbit";
	    }
	}

	int j;
	for(j = 0; j < 3; j++) {
		vo[j] = ( r[j] - F*ro[j] )/G;
    	v[j]  = ( GDot*r[j] - ro[j] )/G;
	}
}
