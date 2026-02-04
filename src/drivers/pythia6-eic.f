C...As an example, consider a main program of the form
C...Double precision and integer declarations.
*======================================================================
      program pythiaeic

      include 'pythia-radcorr/pythia.inc'              ! All PYTHIA commons blocks
      include "pythia-radcorr/mc_set.inc"
      include "pythia-radcorr/py6strf.inc"
      include "pythia-radcorr/mcRadCor.inc"
      include "pythia-radcorr/radgen.inc"
      include "pythia-radcorr/phiout.inc"

C...Added by liang 1/6/12
C...Switches for nuclear correction
      COMMON/PYNUCL/INUMOD,CHANUM,ORDER
      SAVE /PYNUCL/
      DOUBLE PRECISION INUMOD,CHANUM
      INTEGER ORDER

      integer NEV, NPRT, ievent, genevent, I, tracknr, ltype 
      integer lastgenevent, idum1, idum2, initseed, nrtrack
      REAL trueX, trueW2, trueNu
      DOUBLE PRECISION sqrts, radgamE, radgamp, radgamEnucl
      DOUBLE PRECISION pbeamE, pbeta, pgamma, ebeamE, epznucl
      CHARACTER PARAM*100
      CHARACTER PARAM2*100
      LOGICAL UseLut, GenLut

c ---------------------------------------------------------------------
c     Run parameter
c ---------------------------------------------------------------------
      integer*4 today(3), now(3)
c---------------------------------------------------------------------
c     ASCII output file
c ---------------------------------------------------------------------
      integer simcLun, pythiaLun
      parameter (simcLun=29)
      parameter (pythiaLun=30)
      CHARACTER*256 simcfilename
      CHARACTER*256 pythiafilename
      CHARACTER*256 baseout
      CHARACTER*256 outname
      real*8 deg2rad, maxOff

c---------------------------------------------------------------------
C...Function used for random seed generation
c---------------------------------------------------------------------
      integer get_random_seed
      REAL twoPi
      parameter(twoPi=6.28318)
      REAL pi2
      parameter(pi2=3.14159)

      LOGICAL do_simc_out, do_pythia_event_record
      LOGICAL do_proton,do_pion,do_kaon
      logical foundE,foundP,hasPhi
      integer je,jp
      integer keep_pdg, keep_pdgabs_flag, keep_final_flag
      integer excl_flag, excl_rad_flag, ist
      logical do_keep_pdg, keep_abs_pdg, keep_final_only, has_keep_pdg
      logical do_exclusive, excl_allow_rad
      logical is_exclusive_parent
      real*8 excl_mass_tol, excl_mass_flag
      COMMON /EXCLCUTS/ excl_mass_tol
      real*8 pxE,pyE,pzE,eE,pxP,pyP,pzP,eP,phiE,phiP,delPhi,bestPhi
      real*8 vxloc

c---------------------------------------------------------------------
! ... force block data modules to be read
C       external pydata
c ---------------------------------------------------------------------
*  ---- user switches --------------------------------------------------
      do_simc_out = .false.
      do_pythia_event_record = .true.
      do_proton    = .true.  ! default behaviour
      do_pion      = .false.
      do_kaon      = .false.

       iModel=0
       pbeam=100. 
       ebeam=5.0 
       ebeam_original=5.0
       ltype=11
       masse=PYMASS(11)
       massp=PYMASS(2212)
      ievent=0
      genevent=0
      lastgenevent=0
      tracknr=0
      RadState=0 ! Initial or Final State Radiation variable (use with RADGEN - 0 --> No Radiation)
      keep_pdg        = 0
      do_keep_pdg     = .false.
      keep_abs_pdg    = .true.
      keep_final_only = .false.
      do_exclusive    = .false.
      excl_allow_rad  = .false.
      excl_mass_tol   = 0.05d0

C...Read output file name
       READ(*,*) outname
C...Read lepton beam type 
       READ(*,*) ltype 
C...Read parameters for PYINIT call (beam and target particle energy).
       READ(*,*) pbeam, ebeam       
C...Read number of events to generate, and to print.
       READ(*,*) NEV,NPRT
C...Read min/max x of radgen lookup table
       READ(*,*) mcSet_XMin, mcSet_XMax
C...Read min/max y of generation range      
       READ(*,*) mcSet_YMin, mcSet_YMax
C...Read min/max Q2 of generation range      
       READ(*,*) mcSet_Q2Min, mcSet_Q2Max
C...Read information for cross section used in radgen
       READ(*,*) genSet_FStruct, genSet_R
C...Read parameters of radcorr: do radcorr (1), generate look-up table (2)
       READ(*,*) qedrad
C...Read parameters for PYTHIA-Model = which generation is done     
       READ(*,*) iModel
C...Read target type mass and charge
       READ(*,*) mcSet_TarA, mcSet_TarZ
C...Read nuclear pdf parameter mass number A, charge number Z
       READ(*,*) INUMOD, CHANUM
C...Read nuclear pdf correction order
       READ(*,*) ORDER
C...Read information for cross section used in radgen
  100  READ(*,'(A)',END=200) PARAM
       PARAM2 = ADJUSTL(PARAM)

C--- Custom cards (do NOT forward to PYGIVE) ---------------------------
C    Syntax (space-separated; optional '=' allowed):
C      KEEP_PDG        223        (0 disables)
C      KEEP_PDGABS     1          (1=abs match, 0=exact sign)
C      KEEP_PDG_FINAL  1          (1=final-state only, i.e. K(I,1)=1)
C      KEEP_EXCLUSIVE  1          (require exclusive ep -> ep + KEEP_PDG)
C      EXCL_ALLOW_RAD  1          (allow radiative events if RadState!=0)
C      EXCL_MASS_TOL   0.05       (GeV, hadronic system mass tolerance)
C      DO_SIMC         1          (1=write SIMC .dat output)
C      DO_PYTHIA_EVENT_RECORD 1   (1=write PYTHIA event record .txt)
C-----------------------------------------------------------------------

       if (PARAM2(1:8) .eq. 'KEEP_PDG') then
          ist = 9
          if (PARAM2(9:9) .eq. '=') ist = 10
          read(PARAM2(ist:),*) keep_pdg
          do_keep_pdg = (keep_pdg .ne. 0)
          goto 100
       else if (PARAM2(1:11) .eq. 'KEEP_PDGABS') then
          ist = 12
          if (PARAM2(12:12) .eq. '=') ist = 13
          read(PARAM2(ist:),*) keep_pdgabs_flag
          keep_abs_pdg = (keep_pdgabs_flag .ne. 0)
          goto 100
       else if (PARAM2(1:14) .eq. 'KEEP_PDG_FINAL') then
          ist = 15
          if (PARAM2(15:15) .eq. '=') ist = 16
          read(PARAM2(ist:),*) keep_final_flag
          keep_final_only = (keep_final_flag .ne. 0)
          goto 100
       else if (PARAM2(1:14) .eq. 'KEEP_EXCLUSIVE') then
          ist = 15
          if (PARAM2(15:15) .eq. '=') ist = 16
          read(PARAM2(ist:),*) excl_flag
          do_exclusive = (excl_flag .ne. 0)
          goto 100
       else if (PARAM2(1:13) .eq. 'EXCL_ALLOW_RAD') then
          ist = 14
          if (PARAM2(14:14) .eq. '=') ist = 15
          read(PARAM2(ist:),*) excl_rad_flag
          excl_allow_rad = (excl_rad_flag .ne. 0)
          goto 100
       else if (PARAM2(1:13) .eq. 'EXCL_MASS_TOL') then
          ist = 14
          if (PARAM2(14:14) .eq. '=') ist = 15
          read(PARAM2(ist:),*) excl_mass_flag
          excl_mass_tol = excl_mass_flag
          goto 100
       else if (PARAM2(1:7) .eq. 'DO_SIMC') then
          ist = 8
          if (PARAM2(8:8) .eq. '=') ist = 9
          read(PARAM2(ist:),*) keep_final_flag
          do_simc_out = (keep_final_flag .ne. 0)
          goto 100
       else if (PARAM2(1:22) .eq. 'DO_PYTHIA_EVENT_RECORD') then
          ist = 23
          if (PARAM2(23:23) .eq. '=') ist = 24
          read(PARAM2(ist:),*) keep_final_flag
          do_pythia_event_record = (keep_final_flag .ne. 0)
          goto 100
       end if

C--- Normal PYTHIA cards ------------------------------------------------
       CALL PYGIVE(PARAM)
       GOTO 100
c ---------------------------------------------------------------------
C...Initialize PYTHIA.      
c ---------------------------------------------------------------------
  200  write(*,*) '*********************************************'
       write(*,*) 'NOW all parameters are read by PYTHIA'
       write(*,*) '*********************************************'
C       call PYLIST(11)
C       call PYLIST(12)

C     Getting the date and time of the event generation
        
      call idate(today)   ! today(1)=day, (2)=month, (3)=year
      call itime(now)     ! now(1)=hour, (2)=minute, (3)=second
        
!     Take date as the SEED for the random number generation
       
      initseed = get_random_seed()
      write(6,*) 'SEED = ', initseed
      call rndmq (idum1,idum2,initseed,' ')

        
C     proton is defined in positive z and as target
      P(2,1)=0.0  
      P(2,2)=0.0  
      P(2,3)=pbeam
C     lepton is defined in negative z and as beam
      P(1,1)=0.0  
      P(1,2)=0.0  
      P(1,3)=-ebeam

      if (mcSet_TarZ.eq.0) then
        massp=PYMASS(2112)
      else
        massp=PYMASS(2212)
      endif
      masse=PYMASS(ltype)

      pbeamE=sqrt(pbeam**2+massp**2)
      pbeta=pbeam/pbeamE
      pgamma=pbeamE/massp
      ebeamE=sqrt(ebeam**2+masse**2)
      ebeamEnucl=pgamma*ebeamE-pgamma*pbeta*(-ebeam)
      epznucl=-pgamma*pbeta*(ebeamE)+pgamma*(-ebeam)
      write(*,*) ebeamEnucl, ebeamE, epznucl, -ebeam
      mcSet_EneBeam=sngl(ebeamEnucl)

      sqrts=sqrt(2.*pbeam*ebeam+2.*pbeamE*ebeamE+massp**2+masse**2)
      write(*,*) '*********************************************'
      write(*,*) 'proton beam momentum:', pbeam, 'GeV'
      write(*,*) 'lepton beam momentum:', ebeam, 'GeV'
      write(*,*) 'resulting sqrt(s):', sqrts, 'GeV'
      write(*,*) '*********************************************'

      ebeam_original = ebeam
      ! write(*,*)""
      ! write(*,*)"Initial 'ebeam_original' = ", ebeam_original
      ! write(*,*)"Initial 'mcSet_EneBeam'  = ", mcSet_EneBeam
      ! write(*,*)""

      ! MSTP(199) = 1 runs PYGAGA(5,WTGAGA) for RADGEN
      ! MSTP(199) = 0 runs PYGAGA(3,WTGAGA) for NO RADGEN
       if (iModel.eq.0) then
           UseLUT=.false.
           GenLUT=.false.
           qedrad=0
           MSTP(199)=0
           mcRadCor_EBrems=0.
       elseif (iModel.eq.1) then
         if (qedrad.eq.0) then
             mcRadCor_EBrems=0.
             UseLUT=.false.
             GenLUT=.false.
             MSTP(199)=1
         elseif (qedrad.eq.1) then
             mcRadCor_EBrems=0.
             UseLUT=.true.
             GenLUT=.false.
             MSTP(199)=1
             call radgen_init(UseLUT,GenLUT)
             write(*,*) 'I have initialized radgen'
         elseif (qedrad.eq.2) then
             write(*,*) 'radgen lookup table will be generated'
             mcRadCor_EBrems=0.
             UseLUT=.true.
             GenLUT=.true.
             MSTP(199)=1
             call radgen_init(UseLUT,GenLUT)
             goto 500
         endif
       endif

C     Initilization for collider mode
       IF ((mcSet_TarZ .GE. 1) .AND. (ltype .EQ. 11)) THEN
          CALL pyinit('3MOM', 'gamma/e-', 'p+', WIN)
    
       ELSE IF ((mcSet_TarZ .GE. 1) .AND. (ltype .EQ. 11) .AND. (pbeam .
     &   LT. 1E-3)) THEN
          CALL pyinit('FIXT', 'gamma/e-', 'p+', WIN)
          ! This is not going to get called because the line above will always be true at the same time this one is...
          
       ELSE IF ((mcSet_TarZ .GE. 1) .AND. (ltype .EQ. -11)) THEN
          CALL pyinit('3MOM', 'gamma/e+', 'p+', WIN)
          
       ELSE IF ((mcSet_TarZ .EQ. 0) .AND. (ltype .EQ. -11)) THEN
          CALL pyinit('3MOM', 'gamma/e+', 'n0', WIN)
          
       ELSE IF ((mcSet_TarZ .EQ. 0) .AND. (ltype .EQ. 11)) THEN
          CALL pyinit('3MOM', 'gamma/e-', 'n0', WIN)
          
       ENDIF

C      If we ever want to simulate fixed target we need to change this
C      win=ebeam
C      call pyinit('fixt','gamma/e-','p+', WIN)

c ---------------------------------------------------------------------
c     Open ascii output file(s)
c ---------------------------------------------------------------------
       baseout = outname
       if (len_trim(outname) .ge. 4) then
          if (outname(len_trim(outname)-3:len_trim(outname)) .eq.
     &        '.dat') then
             baseout = outname(1:len_trim(outname)-4)
          else if (outname(len_trim(outname)-3:len_trim(outname)) .eq.
     &        '.txt') then
             baseout = outname(1:len_trim(outname)-4)
          end if
       end if

       if (.not. do_simc_out .and. .not. do_pythia_event_record) then
          write(*,*) 'WARNING: do_simc=0 and do_pythia_event_record=0.'
          write(*,*) 'WARNING: no output file will be written.'
       end if

       if (do_simc_out) then
          simcfilename = trim(baseout)//'.dat'
          open(simcLun, file=simcfilename)
          write(*,*) 'SIMC output file: ', trim(simcfilename)
       end if

       if (do_pythia_event_record) then
          pythiafilename = trim(baseout)//'.txt'
          open(pythiaLun, file=pythiafilename)
          write(*,*) 'PYTHIA event record file: ', trim(pythiafilename)
       end if

c ---------------------------------------------------------------------
C...Event generation loop
c ---------------------------------------------------------------------

C   This is what we write in the ascii-file
       if (do_pythia_event_record) then
        write(pythiaLun,*)' PYTHIA EVENT FILE '
        write(pythiaLun,*)'============================================'
        write(pythiaLun,30)
30      format('I, ievent, genevent, subprocess, nucleon,
     &  targetparton, xtargparton, beamparton, xbeamparton,
     &  thetabeamprtn, truey, trueQ2, truex, trueW2, trueNu, leptonphi, 
     &  s_hat, t_hat, u_hat, pt2_hat, Q2_hat, F2, F1, R, sigma_rad, 
     &  SigRadCor, EBrems, photonflux, t-diff, nrTracks')
        write(pythiaLun,*)'============================================'

        write(pythiaLun,*)' I  K(I,1)  K(I,2)  K(I,3)  K(I,4)  K(I,5)
     &  P(I,1)  P(I,2)  P(I,3)  P(I,4)  P(I,5)  V(I,1)  V(I,2)  V(I,3)'
        write(pythiaLun,*)'============================================'
      end if
        ! write(*,*) ''
        ! write(*,*) 'New Event:'
        ! write(*,*) 'mcRadCor_EBrems                       = ', mcRadCor_EBrems

       DO 300 IEV=1,NEV
         ! write(*,*) ""
C         if(NEV.LE.500) then ! These lines are mainly meant for debugging
C             write(*,*) ""
C             write(*,*) "event generated:", IEV
C         end if
         IF(IEV.EQ.1) THEN
             goto 999
         ELSE
             ! write(*,*)"ebeam (before new event)=",ebeam
             ebeam = ebeam_original
             P(2,1)=0.0  
             P(2,2)=0.0  
             P(2,3)=pbeam
             P(1,1)=0.0  
             P(1,2)=0.0  
             P(1,3)=-ebeam
             pbeamE=sqrt(pbeam**2+massp**2)
             pbeta=pbeam/pbeamE
             pgamma=pbeamE/massp
             ebeamE=sqrt(ebeam**2+masse**2)
             ebeamEnucl=pgamma*ebeamE-pgamma*pbeta*(-ebeam)
             epznucl=-pgamma*pbeta*(ebeamE)+pgamma*(-ebeam)
             mcSet_EneBeam=sngl(ebeamEnucl)
             sqrts=sqrt(2.*pbeam*ebeam+2.*pbeamE*ebeamE+massp**2
     &          +masse**2)
             RadState=0
             ! WRITE(*,*)"RESETTING PYINKI(0)"
             CALL PYINKI(0)
             IF(MINT(141).NE.0.OR.MINT(142).NE.0) CALL PYGAGA(1,WTGAGA)
         END IF
999      CALL PYEVNT
         if (MSTI(61).eq.1) then
            ! WRITE(*,*)"MSTI(61).eq.1"
            ebeam = ebeam_original
            P(2,1)=0.0  
            P(2,2)=0.0  
            P(2,3)=pbeam
            P(1,1)=0.0  
            P(1,2)=0.0  
            P(1,3)=-ebeam
            pbeamE=sqrt(pbeam**2+massp**2)
            pbeta=pbeam/pbeamE
            pgamma=pbeamE/massp
            ebeamE=sqrt(ebeam**2+masse**2)
            ebeamEnucl=pgamma*ebeamE-pgamma*pbeta*(-ebeam)
            epznucl=-pgamma*pbeta*(ebeamE)+pgamma*(-ebeam)
            mcSet_EneBeam=sngl(ebeamEnucl)
            sqrts=sqrt(2.*pbeam*ebeam+2.*pbeamE*ebeamE+massp**2
     &          +masse**2)
            RadState=0
            ! WRITE(*,*)"RESETTING PYINKI(0)"
            CALL PYINKI(0)
            IF(MINT(141).NE.0.OR.MINT(142).NE.0) CALL PYGAGA(1,WTGAGA)
            write(*,*) 'Go back to PYEVNT call'
            write(*,*) 'Rerunning event:', IEV
            write(*,*) ''
            goto 999
         endif
         IF (IEV.LE.NPRT) THEN
             CALL PYLIST(2)
             ! write(*,*)""
             ! write(*,*)""
             ! write(*,*)""
             ! write(*,*) '*********************************************'
             ! write(*,*) 'Current lepton beam momentum:', ebeam, 'GeV'
             write(*,*) "event:     ", IEV
             ! ! write(*,*) "K(1,1) Sta=", K(1,1)
             ! ! write(*,*) "K(1,2) PID=", K(1,2)
             ! write(*,*) "P(1,1)  px=", P(1,1)
             ! write(*,*) "P(1,2)  py=", P(1,2)
             ! write(*,*) "P(1,3)  pz=", P(1,3)
             ! write(*,*) "P(1,4)   E=", P(1,4)
             ! ! write(*,*) "P(1,5) Mas=", P(1,5)
             write(*,*) "RadState  =", RadState
             ! write(*,*) '*********************************************'
             ! write(*,*)""
             write(*,*)""
             write(*,*)""
         END IF
         ievent=IEV
         genevent=NGEN(0,3)-lastgenevent

       trueX =  VINT(307)/VINT(309)/(sqrts**2-massp**2)
       trueW2 = massp**2 + VINT(307)*(1/trueX-1)
       trueNu = (trueW2 + VINT(307) - massp**2)/(2.*massp)
      if (mcRadCor_EBrems.gt.0.) then
         radgamEnucl=sqrt(dplabg(1)**2+dplabg(2)**2+dplabg(3)**2)
         radgamE=pgamma*radgamEnucl-pgamma*pbeta*dplabg(3)
         radgamp=-pgamma*pbeta*radgamEnucl+pgamma*dplabg(3)
C         write(*,*) radgamEnucl, radgamE, dplabg(3), radgamp
      else
        radgamEnucl=0D0
        radgamE=0D0
        radgamp=0D0 
      endif

         tracknr=N
         if (mcRadCor_EBrems.gt.0.) then
            nrtrack=tracknr+1
         else
            nrtrack=tracknr
         endif

C----- Optional event filter: require PDG to appear in event -----------
         if (do_keep_pdg) then
            has_keep_pdg = .false.
            do I = 1, tracknr
               if (keep_final_only) then
                  if (K(I,1) .ne. 1) cycle
               end if

               if (keep_abs_pdg) then
                  if (abs(K(I,2)) .eq. abs(keep_pdg)) then
                     has_keep_pdg = .true.
                     exit
                  end if
               else
                  if (K(I,2) .eq. keep_pdg) then
                     has_keep_pdg = .true.
                     exit
                  end if
               end if
            end do

            if (.not.has_keep_pdg) cycle
         end if

C----- Exclusive selection: ep -> ep + KEEP_PDG (no other particles) ---
         if (do_exclusive) then
            if (.not. do_keep_pdg) then
               write(*,*) 'ERROR: KEEP_EXCLUSIVE requires KEEP_PDG != 0'
               stop
            end if

            if (.not. excl_allow_rad) then
               if (RadState .ne. 0) cycle
            end if

            if (.not. is_exclusive_parent(keep_pdg, .false.)) cycle
         end if

         if ((msti(1).ge.91).and.(msti(1).le.94)) msti(16)=0
        if (do_pythia_event_record) then    
         write(pythiaLun,32) 0, ievent, genevent, msti(1), msti(12), 
     &        msti(16), pari(34), msti(15), pari(33), pari(53), 
     &        VINT(309), VINT(307), trueX, trueW2, trueNu,
     &        VINT(313), pari(14), pari(15), pari(16), 
     &        pari(18),  pari(22), sngl(py6f2), sngl(py6f1), 
     &        py6r, mcRadCor_Sigrad, mcRadCor_sigcor, radgamEnucl,
     &        VINT(319), VINT(45), nrtrack 
 32      format((I4,1x),(I10,1x),3(I4,1x),(I10,1x),f9.6,1x,
     &         I12,1x,
     &         2(f12.6,1x),7(f18.11,3x),12(f19.9,3x),I12)
         write(pythiaLun,*)'===================='//
     &        '========================'

         DO I=1,tracknr
         if (K(I,3).le.nrtrack) then
         write(pythiaLun,34) I,K(I,1),K(I,2),K(I,3),K(I,4),K(I,5),
     &        P(I,1),P(I,2),P(I,3),P(I,4),P(I,5),
     &        V(I,1),V(I,2),V(I,3)
         endif
         ENDDO
!          if (mcRadCor_EBrems.gt.0.) then
!             if (RadState.eq.1) then
!                 ! Initial State Radiated Photon ==> RadState = 1, Rad_Photon_Status = 55
!                 write(*,*)""
!                 write(*,*) CHAR(27)//'[1;34m'//
!      &      'RadState.eq.1'//CHAR(27)//'[0m'
!                 write(*,*)""
!                 write(29,34) nrtrack, 55, 22, 1, 0, 0,
!      &      sngl(dplabg(1)),sngl(dplabg(2)),sngl(-radgamp),
!      &      sngl(radgamE), 0., 0., 0., 0.
!             else
!                 ! Final State Radiated Photon ==> RadState = 2, Rad_Photon_Status = 56
!                 write(29,34) nrtrack, 56, 22, 1, 0, 0,
!      &      sngl(dplabg(1)),sngl(dplabg(2)),sngl(-radgamp),
!      &      sngl(radgamE), 0., 0., 0., 0.
!             endif
!          endif
         if (mcRadCor_EBrems.gt.0.) then
            if (RadState.eq.0) then
                ! Non-Radiated Radiated Photon (error in status code = 54 - see below)
                write(*,*)""
                write(*,*) CHAR(27)//'[1;34m'//
     &      'RadState.eq.0'//CHAR(27)//'[0m'
                write(*,*)""
            endif
            ! ERROR IN Radiated Photon      ==> RadState = 0, Rad_Photon_Status = 54 (Did not identify photon in ISR/FSR calculation --> need to fix bugs where these can occur)
            ! Initial State Radiated Photon ==> RadState = 1, Rad_Photon_Status = 55 
            ! Final State Radiated Photon   ==> RadState = 2, Rad_Photon_Status = 56
            ! Unknown State Radiated Photon ==> RadState = 3, Rad_Photon_Status = 57
                ! The beam or target PIDs changed before getting to radgen - cannot modify in case of ISR
                ! Specifically for if the beam and target are not an electron and proton (other configurations not available yet)
                ! If the beam's PID changed, but the event would be FSR, then the photon's status will still be 56 (57 is only used for if it should have been 55 if not for the different PID)
            write(pythiaLun,34) nrtrack, 54+RadState, 22, 1, 0, 0,
     &      sngl(dplabg(1)),sngl(dplabg(2)),sngl(-radgamp),
     &      sngl(radgamE), 0., 0., 0., 0.
         else if (RadState.NE.0) then
             write(*,*)""
             write(*,*)""
             write(*,*)""
             write(*,*) CHAR(27)//'[1;34m'//
     &      "ERROR: RadState=", RadState, " but no photon was written"//
     &      CHAR(27)//'[0m'
             write(*,*)""
             write(*,*)""
             write(*,*)""
         endif
 34      format(2(I6,1x),I10,1x,3(I8,1x),8(f15.6,1x))
         write(pythiaLun,*)'=============== Event finished =============='
      end if
      if (do_simc_out) then
!----- 1) veto φ and find the scattered e⁻ exactly as before ----------
         foundE = .false.
         hasPhi = .false.
         do I = 1, tracknr
            if (abs(K(I,2)) .eq. 333) then
               hasPhi = .true.
               exit
            end if
            if (.not.foundE .and. K(I,1) .eq. 1           
     &        .and. abs(K(I,2)) .eq. 11                
     &        .and. K(I,3) .eq. 3) then
               foundE = .true.
               pxE = P(I,1)
               pyE = P(I,2)
               pzE = P(I,3)
               eE  = P(I,4)
            end if
         ENDDO
         if (hasPhi .or. .not.foundE) cycle

!----- 2) rotate & flip the electron once -----------------------------
         phiE0 = atan2(pyE,pxE)
         if (phiE0 .lt. 0d0) phiE0 = phiE0 + twoPi
         sectorWidth = twoPi / 12.d0
         deltaphi    = 1.57079632679d0                          
     &              - (int(phiE0/sectorWidth) + 0.5d0)*sectorWidth
         cosd = cos(deltaphi)
         sind = sin(deltaphi)
         oldpxE = pxE
         oldpyE = pyE
         pxE = oldpxE*cosd - oldpyE*sind
         pyE = oldpxE*sind + oldpyE*cosd
         pzE = -pzE

!----- 3) pick which PID to output ------------------------------------
         if (do_pion) then
            targetPID = 211
         else if (do_kaon) then
            targetPID = 321
         else
            targetPID = 2212
         end if
      ! ———— compute angular window in radians ——————————————
      deg2rad = pi2/180.d0
      maxOff  = 20.d0*deg2rad

      do I = 1, tracknr

         ! skip if not final-state hadron of chosen PID
         if (K(I,1) /= 1 .or. K(I,2) /= targetPID) cycle

         ! compute Δφ relative to unrotated electron
         phiP0 = atan2(P(I,2),P(I,1))
         delPhi = phiE0 - phiP0
         if (delPhi > pi2) then
            delPhi = delPhi - twoPi
         end if
         if (delPhi < -pi2) then
            delPhi = delPhi + twoPi
         end if

         off = abs(abs(delPhi) - pi2)

         if (off > maxOff) cycle

         ! apply rotation & flip to hadron
         pxP = P(I,1);  pyP = P(I,2)
         oldpxP = pxP;  oldpyP = pyP
         pxP = oldpxP*cosd - oldpyP*sind
         pyP = oldpxP*sind + oldpyP*cosd
         pzP = -P(I,3)
         eP  = P(I,4)

         call random_number(vxloc)
         vxloc = -7.5d0 + 15.d0 * vxloc

         write(simcLun,'(11F12.6)') 
     &      pxP,pyP,pzP,eP,vxloc,
     &      pxE,pyE,pzE,eE,vxloc,1.0

      end do


      end if
      lastgenevent=NGEN(0,3)

  300  CONTINUE
      
C...Print cross sections.
       CALL PYSTAT(1)
       CALL PYSTAT(4)

       write(*,*)"The charm mass used is: ", PMAS(4,1)

C...Print the Pythia cross section which is needed to get an absolut 
C   normalisation the number is in microbarns
       write(*,*)'==================================================='
       write(*,*)'Pythia total cross section normalisation:',
     &            pari(1)*1000, ' microbarn'
       write(*,*)'Total Number of generated events', MSTI(5)
       write(*,*)'Total Number of trials', NGEN(0,3)
       write(*,*)'==================================================='
       if (do_simc_out) then
          close(simcLun)
       end if
       if (do_pythia_event_record) then
          close(pythiaLun)
       end if

  500  if (qedrad.eq.2) then
         write(*,*) 'lookup table is generated;'
         write(*,*) 'to run now pythia change parameter qedrad to 1'
       endif

C...Check pdf status       
       call PDFSTA
       END

      logical function is_exclusive_parent(parent_pdg, allow_photons)
      include 'pythia-radcorr/pythia.inc'
      integer parent_pdg
      logical allow_photons

      integer i, mom
      integer nE, nP, nParent, nParentFinal
      logical ok_from_parent
      real*8 sum_px, sum_py, sum_pz, sum_e, mass2, inv_mass
      real*8 excl_mass_tol
      COMMON /EXCLCUTS/ excl_mass_tol

      nE = 0
      nP = 0
      nParent = 0
      nParentFinal = 0
      sum_px = 0d0
      sum_py = 0d0
      sum_pz = 0d0
      sum_e = 0d0

C-- count parent occurrences anywhere in record (decayed parent may be status 2)
      do i = 1, N
         if (abs(K(i,2)) .eq. abs(parent_pdg)) then
            nParent = nParent + 1
            if (K(i,1) .eq. 1) nParentFinal = nParentFinal + 1
         end if
      end do

      if (nParent .lt. 1) then
         is_exclusive_parent = .false.
         return
      end if
      if (nParentFinal .gt. 1) then
         is_exclusive_parent = .false.
         return
      end if

C-- loop over final-state particles only
      do i = 1, N
         if (K(i,1) .ne. 1) cycle

C-- allow the scattered lepton and recoil proton
         if (abs(K(i,2)) .eq. 11) then
            nE = nE + 1
            cycle
         end if
         if (K(i,2) .eq. 2212) then
            nP = nP + 1
            cycle
         end if

C-- optionally allow final-state photons (from PYJETS)
         if (K(i,2) .eq. 22 .and. allow_photons) cycle

C-- allow the parent itself if you kept it stable (parent as final-state)
         if (abs(K(i,2)) .eq. abs(parent_pdg)) then
            sum_px = sum_px + dble(P(i,1))
            sum_py = sum_py + dble(P(i,2))
            sum_pz = sum_pz + dble(P(i,3))
            sum_e  = sum_e  + dble(P(i,4))
            cycle
         end if

C-- otherwise: require this final particle to descend from the parent
         ok_from_parent = .false.
         mom = K(i,3)
  10     continue
         if (mom .le. 0) goto 20
         if (abs(K(mom,2)) .eq. abs(parent_pdg)) then
            ok_from_parent = .true.
            goto 20
         end if
         mom = K(mom,3)
         goto 10

  20     continue
         if (.not. ok_from_parent) then
            is_exclusive_parent = .false.
            return
         end if
         sum_px = sum_px + dble(P(i,1))
         sum_py = sum_py + dble(P(i,2))
         sum_pz = sum_pz + dble(P(i,3))
         sum_e  = sum_e  + dble(P(i,4))
      end do

C-- require exactly one scattered lepton and one proton
      if (nE .eq. 1 .and. nP .eq. 1) then
         if (sum_e .le. 0d0) then
            is_exclusive_parent = .false.
            return
         end if
         mass2 = sum_e**2 - sum_px**2 - sum_py**2 - sum_pz**2
         if (mass2 .lt. 0d0) mass2 = 0d0
         inv_mass = sqrt(mass2)
         if (abs(inv_mass - dble(PYMASS(parent_pdg))) .le.
     &       excl_mass_tol) then
            is_exclusive_parent = .true.
         else
            is_exclusive_parent = .false.
         end if
      else
         is_exclusive_parent = .false.
      end if

      return
      end

      ! ================================================================
      ! Return a kind=8 integer from the system clock.
      ! Use the system_clock() intrinsic function if possible.
      ! In case of failure there, manually construct a time in
      ! milliseconds using date_and_time().
      ! Returns 0 in case of any problems.
      ! ================================================================
      function query_clock()

         implicit none

         integer*8 query_clock

         integer date(8)

         ! Get processor clock using system_clock() intrinsic.
         ! Resultant value can be in milli-, micro- or nano-seconds
         ! depending on the platform and the compiler.
         ! Don't use an integer(kind=4) to avoid wraparound.
         call system_clock(query_clock)
         ! In case of no clock, or failure to query the clock,
         ! the value is set to -huge(). In this case, manually
         ! construct a time using date information.
         if (query_clock .eq. -huge(query_clock)) then
            ! date_and_time was introduced in the Fortran95 standard.
            ! date(1) to date(8) are, in order:
            ! year, month, day, time difference with UTC in minutes,
            ! hour, minutes, seconds, milliseconds.
            call date_and_time(values=date)
            ! Milliseconds since 1970
            query_clock =
     +         (date(1) - 1970) * 365_8 * 24 * 60 * 60 * 1000 +
     +          date(2) * 31_8 * 24 * 60 * 60 * 1000 +
     +          date(3) * 24 * 60 * 60  * 1000 +
     +          date(5) * 60 * 60 * 1000 +
     +          date(6) * 60 * 1000 +
     +          date(7) * 1000 +
     +          date(8)
         end if
         ! Requre clock to be positive, otherwise set
         ! an error value.
         if (query_clock .lt. 0) then
            query_clock = 0
         end if
      end function query_clock
      ! ================================================================


      ! ================================================================
      ! Returns an integer seed from an integer(kind=8) clock value.
      ! ================================================================
      function make_clock_seed(clock)

         implicit none

         integer make_clock_seed

         integer*8 clock
         integer temp(2)

         ! Print a warning if we got a zero clock value.
         if (clock .eq. 0) then
            write(*,*) 'Warning: making clock seed from clock value 0'
         end if
         ! Transfer clock bitwise representation to pair of integers
         temp = transfer(clock, temp)
         ! XOR the two components
         make_clock_seed = ieor(temp(1), temp(2))
         ! Additional XOR with current process ID, plus a prime!
         make_clock_seed = ieor(make_clock_seed, getpid() + 1099279)
      end function make_clock_seed
      ! ================================================================


      ! ================================================================
      ! Get a random number to seed the random number generator.
      ! Impementation based on
      ! http://gcc.gnu.org/onlinedocs/gfortran/RANDOM_005fSEED.html
      ! ================================================================
      function get_random_seed()

         implicit none

         integer get_random_seed

         integer io_status
         integer unit_number
         parameter (unit_number=123)
         ! Functions
         integer*8 query_clock
         integer make_clock_seed

         ! First we try to use the OS-provided random number if available
         open(unit=unit_number, file='/dev/urandom', access='stream',
     +        form='unformatted', action='read', status='old',
     +        iostat=io_status)
         ! If /dev/urandom was opened then read a number from it
         if (io_status == 0) then
            read(unit_number) get_random_seed
            close(unit_number)
         ! That didn't work, so make a seed from the clock and
         ! the process ID (in case we have simultaneous jobs)
         else
            get_random_seed = make_clock_seed(query_clock())
         end if
         if (get_random_seed .lt. 0) then
            get_random_seed = -get_random_seed
         end if
      end function get_random_seed
      ! ================================================================
