!===============================================================================!
MODULE multifluid_riemann_mod
!===============================================================================!
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  !*It consists of two nonlinear acoustic waves (either shocks or rarefactions) and a contact discontinuity which is also a material interface
  !----------------------------------------------------------------------------

  INTERFACE multifluid_riemann
     MODULE PROCEDURE multifluid_riemann
  END INTERFACE multifluid_riemann



  !----------------------------------------------------------------------------
  PUBLIC  :: multifluid_riemann
  !----------------------------------------------------------------------------
CONTAINS
      SUBROUTINE multifluid_riemann(Gmm_L, pinf_L, density_L, velocity_L, pressure_L, Gmm_R, pinf_R, density_R, velocity_R, pressure_R, density_star_L, density_star_R, velocity_star, pressure_star)
                                     !*gl,  pinfl,        rl,         ul,         pl,    gr,  pinfr,        rr,         ur,         pr,             r1,             r2,         ustar,         pstar)
         IMPLICIT NONE
         !*---------------------------------------------
         REAL, INTENT(IN)        :: Gmm_L, pinf_L, density_L, velocity_L, pressure_L
         REAL, INTENT(IN)        :: Gmm_R, pinf_R, density_R, velocity_R, pressure_R
         REAL, INTENT(INOUT)     :: density_star_L
         REAL, INTENT(INOUT)     :: density_star_R
         REAL, INTENT(INOUT)     :: velocity_star, pressure_star
         !*---------------------------------------------
         REAL :: gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,cl,r1,r2,ustar,pstar
         REAL :: aaa, bbb
         REAL :: www1, www2


         !*======================
         !*LEFT STATE SINGLE FLUID
         !*======================
         rl    = density_L
         ul    = velocity_L
         pl    = pressure_L
         gl    = Gmm_L
         pinfl = pinf_L
         cl=SQRT(gl*(pl+pinfl)/rl)
         !*======================
         !*RIGHT STATE SINGLE FLUID
         !*======================
         rr    = density_R
         ur    = velocity_R
         pr    = pressure_R
         gr    = Gmm_R
         pinfr = pinf_R
         cr=SQRT(gr*(pr+pinfr)/rr)


         if(ul.ge.ur) then
            if(pl.lt.pr) then
               aaa=pl
               www1=pr+pinfr
               www2=2.0*(ul-ur)**2*gr*rr
               bbb=pr+MAX(www1,www2)
            else
               aaa=pr
               www1=pl+pinfl
               www2=2.0*(ul-ur)**2*gl*rl
               bbb=pl+MAX(www1,www2)
            endif
            call ss(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,r1,r2,ustar,pstar)
            if((pstar.gt.pr).and.(pstar.gt.pl)) then
               goto 128
            else
               if(pl.lt.pr) then
                  aaa=pl
                  bbb=pr
                  call sr(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,r1,r2,ustar,pstar)
               else
                  aaa=pr
                  bbb=pl
                  call rs(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cl,r1,r2,ustar,pstar)
               endif
            endif
         else
            aaa=1e-6
            bbb=MAX(pl,pr)
            call rrw(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,cl,r1,r2,ustar,pstar)
            if((pstar.lt.pr).and.(pstar.lt.pl)) then
               goto 128
            else
               if(pl.lt.pr) then
                  aaa=pl
                  bbb=pr
                  call sr(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,r1,r2,ustar,pstar)
               else
                  aaa=pr
                  bbb=pl
                  call rs(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cl,r1,r2,ustar,pstar)
               endif
            endif
         endif
         128    continue

         !*State star left
         ! cstarl=SQRT(gl*(pstar+pinfl)/r1)
         density_star_L  = r1

         !*State star right
         ! cstarr=SQRT(gr*(pstar+pinfr)/r2)
         density_star_R  = r2

         velocity_star   = ustar
         pressure_star   = pstar



      END SUBROUTINE multifluid_riemann


      !---------------------------------------------------------------------
      !
      ! Exact Shock-Shock solution of Riemann problems .....................
      !
      !---------------------------------------------------------------------
      subroutine ss(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,r1,r2,ust,pst)
         IMPLICIT NONE
         REAL  :: aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,cl,r1,r2,ust,pst
         INTEGER :: iter, i
         REAL  :: al, ar, bl, br, dmul, dmur, dmul1, dmur1, fal, far, fbl, fbr, fa, fb, rur, rul, fl, fr, f, ccc


         al=1.0/(gl-1.0)
         ar=1.0/(gr-1.0)

         iter=50
         dmul=2.0*gl*al
         dmur=2.0*gr*ar
         bl=(gl+1.0)*al
         br=(gr+1.0)*ar
         fal=SQRT(dmul)*(aaa-pl)/SQRT(gl*rl*(pl+bl*aaa+dmul*pinfl))
         far=SQRT(dmur)*(aaa-pr)/SQRT(gr*rr*(pr+br*aaa+dmur*pinfr))
         fbl=SQRT(dmul)*(bbb-pl)/SQRT(gl*rl*(pl+bl*bbb+dmul*pinfl))
         fbr=SQRT(dmur)*(bbb-pr)/SQRT(gr*rr*(pr+br*bbb+dmur*pinfr))
         fa=ur+far-ul+fal
         fb=ur+fbr-ul+fbl

         do i=1,iter
            ccc=0.5*(aaa+bbb)
            fl=SQRT(dmul)*(ccc-pl)/SQRT(gl*rl*(pl+bl*ccc+dmul*pinfl))
            fr=SQRT(dmur)*(ccc-pr)/SQRT(gr*rr*(pr+br*ccc+dmur*pinfr))
            f=ur-ul+fr+fl
            ust=0.5*(ur+ul+fr-fl)
            if(fa*f.gt.0.0) then
               fa=f
               aaa=ccc
            else
               fb=f
               bbb=ccc
            endif
         enddo
         pst=ccc
         r1=rl*(pl+bl*pst+dmul*pinfl)/(pst+bl*pl+dmul*pinfl)
         r2=rr*(pr+br*pst+dmur*pinfr)/(pst+br*pr+dmur*pinfr)
         return
      END SUBROUTINE
      !*---------------------------------------------------------------------
      !*
      !* Exact Rarefaction-Rarefaction solution of Riemann problems .........
      !*
      !*---------------------------------------------------------------------
      SUBROUTINE rrw(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,cl,r1,r2,ust,pst)
         IMPLICIT NONE
         REAL  :: aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,cl,r1,r2,ust,pst
         INTEGER :: iter, i
         REAL  :: al, ar, bl, br, dmul, dmur, dmul1, dmur1, fal, far, fbl, fbr, fa, fb, rur, rul, fl, fr, f, ccc

         al=1.0/(gl-1.0)
         ar=1.0/(gr-1.0)

         iter=50
         dmul=2.0*gl*al
         dmur=2.0*gr*ar
         dmul1=1.0/dmul
         dmur1=1.0/dmur
         fal=dmul*cl*(1.0-((aaa+pinfl)/(pl+pinfl))**dmul1)/gl
         far=dmur*cr*(1.0-((aaa+pinfr)/(pr+pinfr))**dmur1)/gr
         fbl=dmul*cl*(1.0-((bbb+pinfl)/(pl+pinfl))**dmul1)/gl
         fbr=dmur*cr*(1.0-((bbb+pinfr)/(pr+pinfr))**dmur1)/gr
         fa=ur-far-ul-fal
         fb=ur-fbr-ul-fbl

         do i=1,iter
            ccc=0.5*(aaa+bbb)
            fl=dmul*cl*(1.0-((ccc+pinfl)/(pl+pinfl))**dmul1)/gl
            fr=dmur*cr*(1.0-((ccc+pinfr)/(pr+pinfr))**dmur1)/gr
            f=ur-fr-ul-fl
            ust=0.5*(ul+ur-fr+fl)
            if(fa*f.gt.0.0) then
               fa=f
               aaa=ccc
            else
               fb=f
               bbb=ccc
            endif
         enddo
         pst=ccc
         r1=rl*((pst+pinfl)/(pl+pinfl))**(1.0/gl)
         r2=rr*((pst+pinfr)/(pr+pinfr))**(1.0/gr)
         return
      END SUBROUTINE
      !*---------------------------------------------------------------------
      !*
      !* Exact Shock-Rarefaction solution of Riemann problems ...............
      !*
      !*---------------------------------------------------------------------
      SUBROUTINE sr(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,r1,r2,ust,pst)
         IMPLICIT NONE
         REAL  :: aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,cl,r1,r2,ust,pst
         INTEGER :: iter, i
         REAL  :: al, ar, bl, br, dmul, dmur, dmul1, dmur1, fal, far, fbl, fbr, fa, fb, rur, rul, fl, fr, f, ccc

         al=1.0/(gl-1.0)
         ar=1.0/(gr-1.0)

         iter=50
         dmul=2.0*gl*al
         dmur=2.0*gr*ar
         bl=(gl+1.0)*al
         dmur1=1.0/dmur
         fal=SQRT(dmul)*(aaa-pl)/SQRT(gl*rl*(pl+bl*aaa+dmul*pinfl))
         far=dmur*cr*(1.0-((aaa+pinfr)/(pr+pinfr))**dmur1)/gr
         fbl=SQRT(dmul)*(bbb-pl)/SQRT(gl*rl*(pl+bl*bbb+dmul*pinfl))
         fbr=dmur*cr*(1.0-((bbb+pinfr)/(pr+pinfr))**dmur1)/gr
         fa=ur-far-ul+fal
         fb=ur-fbr-ul+fbl

         do i=1,iter
            ccc=0.5*(aaa+bbb)
            fl=SQRT(dmul)*(ccc-pl)/SQRT(gl*rl*(pl+bl*ccc+dmul*pinfl))
            fr=dmur*cr*(1.0-((ccc+pinfr)/(pr+pinfr))**dmur1)/gr
            f=ur-fr-ul+fl
            ust=0.5*(ur+ul-fr-fl)
            if(fa*f.gt.0.0) then
               fa=f
               aaa=ccc
            else
               fb=f
               bbb=ccc
            endif
         enddo
         pst=ccc
         r1=rl*(pl+bl*pst+dmul*pinfl)/(pst+bl*pl+dmul*pinfl)
         r2=rr*((pst+pinfr)/(pr+pinfr))**(1.0/gr)
         return
      END SUBROUTINE
      !*---------------------------------------------------------------------
      !*
      !* Exact Rarefaction-Shock solution of Riemann problems ...............
      !*
      !*---------------------------------------------------------------------
      SUBROUTINE rs(aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cl,r1,r2,ust,pst)
         IMPLICIT NONE
         REAL  :: aaa,bbb,gl,pinfl,rl,ul,pl,gr,pinfr,rr,ur,pr,cr,cl,r1,r2,ust,pst
         INTEGER :: iter, i
         REAL  :: al, ar, bl, br, dmul, dmur, dmul1, dmur1, fal, far, fbl, fbr, fa, fb, rur, rul, fl, fr, f, ccc

         al=1.0/(gl-1.0)
         ar=1.0/(gr-1.0)

         iter=50
         dmul=2.0*gl*al
         dmur=2.0*gr*ar
         br=(gr+1.0)*ar
         dmul1=1.0/dmul
         fal=dmul*cl*(1.0-((aaa+pinfl)/(pl+pinfl))**dmul1)/gl
         far=SQRT(dmur)*(aaa-pr)/SQRT(gr*rr*(pr+br*aaa+dmur*pinfr))
         fbl=dmul*cl*(1.0-((bbb+pinfl)/(pl+pinfl))**dmul1)/gl
         fbr=SQRT(dmur)*(bbb-pr)/SQRT(gr*rr*(pr+br*bbb+dmur*pinfr))
         fa=ur+far-ul-fal
         fb=ur+fbr-ul-fbl

         do i=1,iter
            ccc=0.5*(aaa+bbb)
            fl=dmul*cl*(1.0-((ccc+pinfl)/(pl+pinfl))**dmul1)/gl
            fr=SQRT(dmur)*(ccc-pr)/SQRT(gr*rr*(pr+br*ccc+dmur*pinfr))
            f=ur+fr-ul-fl
            ust=0.5*(ur+ul+fr+fl)
            if(fa*f.gt.0.0) then
               fa=f
               aaa=ccc
            else
               fb=f
               bbb=ccc
            endif
         enddo
         pst=ccc
         r1=rl*dabs((pst+pinfl)/(pl+pinfl))**(1.0/gl)
         r2=rr*(pr+br*pst+dmur*pinfr)/(pst+br*pr+dmur*pinfr)
         return
      END SUBROUTINE

!===============================================================================!
END MODULE multifluid_riemann_mod
!===============================================================================!
