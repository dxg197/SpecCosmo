! MHD_Analysis
! lapse: Calculates the lapse for Gowdy Spacetime:
! Written by David Garrison
 
#include "cctk.h"  
#include "cctk_Arguments.h"  
#include "cctk_Parameters.h"
	
module gowdylapse_analysis
  implicit none
  private
  public glapse
  	!
  contains
  	!
	subroutine glapse (time, z, alp, bet, lz)
        use cactus_analysis
        DECLARE_CCTK_PARAMETERS
        !
        integer, parameter :: gowdy_expand   = 1
        integer, parameter :: gowdy_collapse = 2
        !
        CCTK_REAL  th0,gc,pi,t,lam,bj0,bj1,by0,by1
        CCTK_REAL, intent(in) :: time,z,lz  
        CCTK_REAL, intent(out) :: alp, bet(3)  
        integer gowdy_type
        !
        if (CCTK_EQUALS(gowdytyp, "expand")) then
           gowdy_type = gowdy_expand
        else if (CCTK_EQUALS(gowdytyp, "collapse")) then
           gowdy_type = gowdy_collapse
        end if
        !
        select case (gowdy_type)

          case (gowdy_expand)
          t = time
          case (gowdy_collapse)
          t = exp(-time)

        end select
   	
   	    bj0 = bessj0(t)
   	    bj1 = bessj1(t)
   	    by0 = bessy0(2.d0*t)
   	    by1 = bessy1(2.d0*t)
            pi  = 4.0d0*atan(1.0d0)
            gc  = nw*2.0d0*pi/lz  
            th0 = gc*z 
   	    
            lam = Amp**2*(-t*bj0*bj1*cos(z)**2 + 0.5d0*(t*bj0)**2 + &
                       0.5d0*(t*bj1)**2)      

	    alp = exp(0.25d0*lam)*t**(-0.25)
	    !alp_z = 2.0*Amp**2*t*bj0*bj1*cos(z)*sin(z)*exp(0.25d0*lam)*t**(-0.25)
        bet(:) = 0.0d0  
	!
	END Subroutine glapse 
    
  CCTK_REAL FUNCTION bessj0(x)
	CCTK_REAL x
	CCTK_REAL ax,xx,z
	DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4, &
 		r5,r6,s1,s2,s3,s4,s5,s6,y 
	SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
		s1,s2,s3,s4,s5,s6
	    
        p1=1.d0
        p2=-.1098628627d-2
        p3=.2734510407d-4
        p4=-.2073370639d-5
        p5=.2093887211d-6

        q1=-.1562499995d-1
        q2=.1430488765d-3
        q3=-.6911147651d-5
        q4=.7621095161d-6
        q5=-.934945152d-7
	
        r1=57568490574.d0
        r2=-13362590354.d0
        r3=651619640.7d0
        r4=-11214424.18d0
        r5=77392.33017d0
        r6=-184.9052456d0

	s1=57568490411.d0
        s2=1029532985.d0
        s3=9494680.718d0
        s4=59272.64853d0
        s5=267.8532712d0
        s6=1.d0
	
        if(abs(x).lt.8.)then 
		y=x**2
		bessj0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
		/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
	else
		ax=abs(x)
		z=8./ax
		y=z**2
		xx=ax-.785398164
		bessj0=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
		*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
	endif
	return
  END FUNCTION bessj0 
	
  CCTK_REAL FUNCTION bessy0(x)
	CCTK_REAL x
	CCTK_REAL xx,z
	DOUBLE PRECISION p1,p2,p3,p4,p5,q1, &
		q2,q3,q4,q5,r1,r2,r3,r4, &
		r5,r6,s1,s2,s3,s4,s5,s6,y
	SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4, &
		r5,r6,s1,s2,s3,s4,s5,s6
        
        p1=1.d0
        p2=-.1098628627d-2
        p3=.2734510407d-4
        p4=-.2073370639d-5
        p5=.2093887211d-6

        q1=-.1562499995d-1
        q2=.1430488765d-3
        q3=-.6911147651d-5
        q4=.7621095161d-6
        q5=-.934945152d-7

	r1=-2957821389.d0
        r2=7062834065.d0
        r3=-512359803.6d0
        r4=10879881.29d0
        r5=-86327.92757d0
        r6=228.4622733d0
		
        s1=40076544269.d0
        s2=745249964.8d0
        s3=7189466.438d0
        s4=47447.26470d0
        s5=226.1030244d0
        s6=1.d0
	
        if(x.lt.8.)then 
		y=x**2
		bessy0=(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y &
		*(s3+y*(s4+y*(s5+y*s6)))))+.636619772*bessj0(x)*log(x)
	else
		z=8./x
		y=z**2
		xx=x-.785398164
		bessy0=sqrt(.636619772/x)*(sin(xx)*(p1+y*(p2+y*(p3+y*(p4+y* &
		p5))))+z*cos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
	endif
	return
  END FUNCTION bessy0 
	
  CCTK_REAL FUNCTION bessj1(x)
	CCTK_REAL x
	CCTK_REAL ax,xx,z
	DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4, &
		r5,r6,s1,s2,s3,s4,s5,s6,y 
	SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6, &
		s1,s2,s3,s4,s5,s6
	    
	r1=72362614232.d0
        r2=-7895059235.d0
        r3=242396853.1d0
        r4=-2972611.439d0
        r5=15704.48260d0
        r6=-30.16036606d0
		
        s1=144725228442.d0
        s2=2300535178.d0
        s3=18583304.74d0
        s4=99447.43394d0
        s5=376.9991397d0
        s6=1.d0
	
        p1=1.d0
        p2=.183105d-2
        p3=-.3516396496d-4 
        p4=.2457520174d-5
        p5=-.240337019d-6
	
        q1=.04687499995d0
        q2=-.2002690873d-3
        q3=.8449199096d-5
        q4=-.88228987d-6
        q5=.105787412d-6

	if(abs(x).lt.8.)then 
		y=x**2
		bessj1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) &
		/(s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))))
	else 
		ax=abs(x)
		z=8./ax
		y=z**2
		xx=ax-2.356194491
		bessj1=sqrt(.636619772/ax)*(cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
		*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*(x/abs(x)) 
	endif
	return
  END FUNCTION bessj1 
	
  CCTK_REAL FUNCTION bessy1(x)
	CCTK_REAL x
	CCTK_REAL xx,z
	DOUBLE PRECISION p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4, &
		r5,r6,s1,s2,s3,s4,s5,s6,s7,y 
	SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4, &
		r5,r6,s1,s2,s3,s4,s5,s6,s7
	    
	p1=1.d0
        p2=.183105d-2
        p3=-.3516396496d-4
        p4=.2457520174d-5
        p5=-.240337019d-6
		
        q1=.04687499995d0
        q2=-.2002690873d-3
        q3=.8449199096d-5
        q4=-.88228987d-6
        q5=.105787412d-6

	r1=-.4900604943d13
        r2=.1275274390d13
        r3=-.5153438139d11
        r4=.7349264551d9
        r5=-.4237922726d7
        r6=.8511937935d4
     
	s1=.2499580570d14
        s2=.4244419664d12
        s3=.3733650367d10
        s4=.2245904002d8
        s5=.1020426050d6
        s6=.3549632885d3
        s7=1.d0

	if(x.lt.8.)then 
		y=x**2
		bessy1=x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6)))))/(s1+y*(s2+y* &
		(s3+y*(s4+y*(s5+y*(s6+y*s7))))))+.636619772 &
		*(bessj1(x)*log(x)-1./x)
	else
		z=8./x
		y=z**2
		xx=x-2.356194491
		bessy1=sqrt(.636619772/x)*(sin(xx)*(p1+y*(p2+y*(p3+y*(p4+y &
		*p5))))+z*cos(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))
	endif
	return
  END FUNCTION bessy1 

END MODULE gowdylapse_analysis
