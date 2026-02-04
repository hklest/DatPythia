# PYTHIA-RAD-CORR

## About

Documentation in the BNL Wiki
* https://wiki.bnl.gov/eic/index.php/PYTHIA
Note that repository and other technical information in the wiki are superseded by
this document.

Contacts
* Kolja Kauder <kkauder@bnl.gov>
* Elke Aschenauer <elke@bnl.gov>

Based on PYTHIA 6.4.28 with radiative corrections and output mdifications.

## Prerequisites
- gfortran
- LHAPDF5
- ROOT
- cmake
- CERNLIB

# Basic installation (on RCF, adapt as appropriate):

```shell
cd ${EICDIRECTORY}/PACKAGES
git clone https://gitlab.com/eic/mceg/PYTHIA-RAD-CORR.git
cd ${EICDIRECTORY}/PACKAGES/PYTHIA-RAD-CORR
mkdir build; cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=${EICDIRECTORY}
make -j 8 install
```
This will produce warnings of the form
```
Warning: $ should be the last specifier in format at (1)
```
which should be okay (it is a g77 extension allowed by gfortran).

# Example Run
To run, the pythia6-eic executable reads options from standard input
and is controlled via steering cards and input redirection, with
optional output redirection to a log file. The name of the output file
is specified in the steering card.
(Note: For backwards compatibility we create a symbolic link pythiaeRHIC
  to the renamed executable)
```
pythia6-eic < STEER_FILE > out.log
```

A practical example, assuming $EICDIRECTORY was set and the package
installed as above, is:
```shell
$EICDIRECTORY/bin/pythia6-eic < $EICDIRECTORY/PACKAGES/PYTHIA-RAD-CORR/STEER-FILES-Other/ep_noradcor.20x250.quicktest
```

# Steering files
A variety of steering files can be found in the STEER-*
directories. Which ones appropriately reflect the physics you wish to
implement is beyond the scope of this document. Please contact the EIC
UG for advice.

Note: The first few lines of the STEER file ultimately take
precedence, but if you want to for example change Q^2, it may be
helpful to also change CKIN(65) too at the bottom, otherwise it will be very slow.


#### Drivers
- pythia6-eic.f
- pythiaMain.cpp uses mainly class interface functions for setup
- UsingCardPythiaMain.cpp can be used in place of the fortran driver


## CHANGES wrt build_pythia6.sh script :

### LHAPDF-related:
- pdfset.f structm.f structp.f are removed (renamed)
- sugra.f is emptied out

### Fix for older pythia versions:
- pyalps.f doesn't need to be changed in 6.4.28
   (build script only affects commented out old code)

### Random number generation change:
Replace default routine by call to ranlux

- pyr.f
Added ranlux.f v 1.2 (Deleted an include line) to avoid issues linking
against cernlib

### PHYSICS:
- pydiff.f
- pydisg.f
- pygaga.f
- pyrand.f
- pyremn.f
- pysgqc.f
- pysigh.f
- pyxtot.f

### KK Note: On 2/26/2020, I pulled in
pydata.f
pydisg.f
pyremn.f
pyptdirc.f
from the original package (they seem to have been overlooked initially)


### Added in top directory:
#### Physics:
- radgen.f
- radgen_event.f
- radgen_init.f
- gmc_random.f
- pyth_xsec.f
- pythia_radgen_extras.f

#  Hard-coding against lhapdf5
Executables require LHAPDF5 installed and links executables against it.


####################### Modified Sep 27, 2024 by Henry Klest ################################

Updated CMakeLists.txt to find an existing version of LHAPDF5. This was very challenging because the majority of the versions of LHAPDF5 available on CVMFS were compiled with the wrong version of FORTRAN. It's crucial to source the setup_radcorr.sh or setup_radcorr.csh scripts before running so the proper version of LHAPDF5 can be found.

The actual edit I made to the code was updating the calculation of Bjorken x in src/radcorr/radgen_event.f. Previously it was calculated as x = Q2/(y*s), where the squared center of mass energy s was calculated as 4*p_p*p_e, which is a valid approximation for a collider but in fixed target mode where the proton mass is not negligible it breaks everything. 

Build instructions on the iFarm:

First modify setup_radcorr.sh to use your preferred directory locations instead of my work directory

Then:

source setup_radcorr.sh;
mkdir build;
cd build;
cmake ../ -DCMAKE_INSTALL_PREFIX=../;
make -j 8 install

# Radiated Beam Energy Corrections
- Edits by Richard Capobianco (as of 1-17-2025)
- Changes Made:
    - Events that use RADGEN to (potentially) produce a radiated photon will now differentiate between initial-state (ISR) and final-state radiation (FSR)
        - ISR vs FSR is determined by checking whether the radiated photon given by RADGEN is closer to the incoming beam or to the scattering electron (using the polar angles)
        - Status code of the radiated photons will be 55 for ISR and 56 for FSR
            - Status code 57 is used for events that would be ISR, but because the beam (K(1,2)) changed particles, the beam energy was not corrected (still investigating how to handle)
                - The change to K(1,2) happens occasionally if PYRAND is called more than once, so the investigation is centered around determining if the new PID has any impact on the event (if so, then the ISR corrections cannot be performed)
                - This caused a related issue in cases where PYRAND was called again after the beam energy was already changed, as the ISR correction code would not be able to reset the beam energy despite RADGEN being called again to generate new photons. Currently handled by killing these events unless the beam energy is still the original requested value
                    - This issue caused some events to be suck in infinite loops
        - No changes to the particle kinematics are (currently) made for FSR events, though these events will need to update the scattered electron's energy/momentum to account for emitting the radiated photon before the events are given to the detector simulator
        - For ISR events, momentum/energy conservation is used to modify the incoming beam's kinematics such that the requested beam emits the radiated photon
            - The beam at the vertex/Born level of the event will be lower than the initial requested amount and is free to move along the x and y axis (instead of just along z)
                - The rest of the event should conserve energy and momentum with the new beam
            - The beam is reset/modified using the PYINKI(0) routine followed by rerunning lines used by calling PYGAGA(1,WTGAGA)
                - Modifications were made to PYGAGA(5,WTGAGA)
                - The beam is reset before every call of 'radgen_event' and after writing the event to the output file (i.e., before calling PYEVNT in pythia6-eic.f)
- Known Bug/Issue(s):
    - Currently designed for just events with an electron beam and a proton target
        - Should not matter whether the proton is fixed or not, though all testing was done assuming 'pbeam' was 0
        - To update the incoming beam energy, the modifications to PYGAGA(5,WTGAGA) require the beam be an electron and the target be a proton
            - PID requirements are hard coded into the new lines of code
            - Some processes may change the electron beam into a photon before reaching PYGAGA(5,WTGAGA) which leads to events where ISR radiation cannot be handled with the current code
    - Also only tested for beam initialization with '3MOM' to define the momentums with P(I,1-3)
    - Interactions with other initial conditions/setups are not currently clear (still open to investigation)
        - Will make an attempt to generalize more for future use in later updates

Modified the implementation of RADGEN to account for the radiated photon's impact on the incoming beam if the radiation occures before the collision between the beam and target (determined based on whether the radiated photon would be closer to the incoming beam or to scattering electron). Previously, the event generator would use the requested beam energy always, however, the modifications made now allow the beam energy to vary based on the radiated photon given by RADGEN such that the requested beam is what radiates the photon rather than what is produced by the initial-state radiative process. Another status marker has also been added to the radiated photon given by RADGEN where the status code 56 now corresponds to when the photon (likely) radiated off of the scattered electron, while the status code 55 corresponds to the initial state radiations off of the incoming beam which impact the beam's energy.


