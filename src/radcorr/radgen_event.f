
      subroutine radgen_event

      implicit none

       include "mc_set.inc"
       include "mconsp.inc"
       include "phiout.inc"
       include "tailcom.inc"
       include "cmpcom.inc"
       include "radgen.inc"
       include "mcRadCor.inc"

      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200) 
      INTEGER MSTP,MSTI
      REAL PARP,PARI
      SAVE /PYPARS/

      real PhRAD(4),q2true,nutrue,radweight
      real nu,q2,phi,yys,xxs
      DOUBLE PRECISION pbeamE, ebeamE, s

! calculate radiative corrections
      nu=sngl(gennu)
      q2=sngl(genq2)
      phi=sngl(genphi)
      yys=sngl(geny)
      xxs=sngl(genx)
      call RADGEN(mcSet_EneBeam,q2,nu,yys,xxs,phi,PhRAD,q2true,nutrue,
     +     radweight)

C      write(*,*)"PythiaQ2= ",genq2," RadgenQ2= ",q2true

C      write(*,*)"PythiaQ2= ",genq2," RadgenQ2= ",q2true
C      write(*,*)"PythiaNu= ",gennu," RadgenNu= ",nutrue
C      write(*,*)"Pythiay = ", gennu/mcSet_EneBeam, " Radgeny = ", nutrue/mcSet_EneBeam


! fill mcRadCor WCB with ...

      mcRadCor_ID=1

! ... true kinematics

* by definition we calculate xbj using the proton mass
* ---> xbj for elastic events is A
*     mcRadCor_XTrue=q2true/(2.0d0*amp*nutrue)
      mcRadCor_NuTrue=nutrue
      mcRadCor_Q2True=q2true
      mcRadCor_YTrue=nutrue/mcSet_EneBeam
      pbeamE=sqrt(pbeam**2+amp**2)
      ebeamE=sqrt(ebeam**2+masse**2)
      s = pbeamE**2 + 2.0*pbeamE*ebeamE + ebeamE**2 - (pbeam - ebeam)**2
      mcRadCor_XTrue= 
     + q2true/(2*nutrue*amp)

      mcRadCor_W2True=amp2 + (q2true*(1./mcRadCor_XTrue-1.))
C      
!     mcRadCor_XTrue=q2true/(2.0d0*0.938272d0*nutrue)
!      mcRadCor_W2True=amp2 - q2true + 2.*amp*nutrue

! ... kinematics of real photon
      mcRadCor_EBrems=phrad(4)
      mcRadCor_ThetaBrems=0.
      if(phrad(4).gt.0.) then
         mcRadCor_ThetaBrems = acos(phrad(3)/phrad(4))
C         write(*,*)"RealE=",deg,"RealTheta=",dthg,"realphi=",dphig
      endif
      mcRadCor_PhiBrems=0.
      
      if (.not.(phrad(1).eq.0..and.phrad(2).eq.0.)) then
        mcRadCor_PhiBrems = atan2(phrad(2),phrad(1))
        if (mcRadCor_PhiBrems.lt.0.)
     +       mcRadCor_PhiBrems = mcRadCor_PhiBrems + twopi
      endif

c...  if we would like to have the TSAI system angles for the real gamma than
C      mcRadCor_ThetaBrems = dthg
C      mcRadCor_PhiBrems = dphig
  
! ... radiative contributions

      mcRadCor_Sigrad=sigrad
      mcRadCor_Sigcor=sigcor
      mcRadCor_Sigcorerr=0.
      mcRadCor_TailIne=tine
      mcRadCor_TailEla=tpro
      mcRadCor_TailCoh=tnuc
      mcRadCor_Vacuum=vac
      mcRadCor_Vertex=vertex
      mcRadCor_Small=small
      mcRadCor_Redfac=redfac

! ... radiative correction type

      if (ita.eq.2) then
        mcRadCor_cType='elas'
      else if (ita.eq.3) then
        mcRadCor_cType='qela'
      else
        mcRadCor_cType='inel'
      endif
C      if(ita.eq.5) write(*,*)"ita=",ita
C      write(*,*)"mcRadCor_cType=",mcRadCor_cType
      end
