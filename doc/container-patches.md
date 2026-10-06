#################################################################################################
Patches for SKSNSim to work the SNBurstPipeline project                                                    
SKSNSim forked from https://github.com/SKSNSim/SKSNSim (Using the most current version as of 13 July 2026)  #
#################################################################################################

##################
Nikolas TC Boily #
13 July 2026     #
##################

Upstream bases for the 2026-07 upgrade:
  SKSNSim: https://github.com/SKSNSim/SKSNSim.git
           branch main @ ea3e3a3 (2026-06-01, "bug modification: the configuration is
           defaulting to DSNB configuration")

download_data.sh:
    Update data dir to fit within the container and extract the values to the top level.

main_snburst.cc:
    Update path of TRandom3 in the #include block to work in the container root path.

reweight.cc:
    Update path of TFile and TTree in the #include block for container compatability of root.

GNUmakefile:
    SKINTERNAL flag is needed for vector files to generate fully. Without it enabled, the 
    vector files do not fill the MCInfo and SNEvInfo, which leaves empty branches and
    incompatabilities with SKG4.

    Define SKINTERNAL constantly to enable these methods (which needs specific undef later
    to bypass fortran code that is not available in the container).

    Additional edits:

        Statically link the SKOFL libraries at compile so it can be built with the container's
        structure.

        Add SKROOT_LIBDIR

        Include path is $(SKOFL_ROOT)/include (env-driven); the no-SKOFL_ROOT fallback
        branch links and rpaths against /mnt/sim/deps/skofl/r32697/lib (updated from the
        old merged cern/2005 tree in the v1.2.0 deps layout split).


SKSNSimCrosssection.cc:
    Undef SKINTERNAL to bypass Fortran routines and use the C++ implementations that are available
    within the SKOFL libraries in the container.


SKSNSimVectorGenerator.cc:
    Update path of TRandom3 to the container root path.

    Undefine SKINTERNAL to bypass the fortran routine for sn_sun_dir_

    Bool m_sn_dir_set added for the sn direction functionality

    Line 452:
        Fixed #ifdef SKITERNAL typo --> SKINTERNAL
        Caused the sn_sundir_ branch to be dead code within the SKOFL build.

    Line 441-451: 
        Added block to set sdir if m_sn_dir_set is true.

    Line 642-646: 
        Added progress report to Tbins and Ebins.

    Lines 1140-1142: 
        Added progress report for fillevent.

SKSNSimVectorGenerator.hh:
    Line 329:
        Added bool m_sn_dir_set

    Line 391-398:
        Added void function SetSNDir to set sn dir.
        Added GetSNDir and GetSNDirSet helper functions.

SKSNSimUserConfiguration.cc:
    Added help lines for the sndir implementations
    Added SNDir implementation for user configuration to be stored for vector generation.

SKSNSimUserConfiguration.hh:
    Added double m_sn_dir[3] and bool m_sn_dir_set for sn direction functionality 

    Lines 115-120: 
        Set default values for sn dir parameters

    Added GetDefaultSNDirX(Y/Z) for default direction

    Lines 225-230: 
        Added sndir SKSNSimUserConfiguration processes for setting the sn direction. Also added
        helper functions GetSNDir and GetSNDirSet to return their states.

SKSNSimFileIO.cc:
    Updated path for TFile and TTree to match container paths for root.

SKSNSimTools.cc:
    Undef skinteranl to bypass fortran routine for elapseday.


Relevant for the SNBurstPipeline container:
setup/build_libsnevtinfo.sh (not SKSNSim source, but part of its container build chain):
    Rebuilds libsnevtinfo.so from ${SKOFL_ROOT}/src/skroot inside the container; paths
    updated to the deps/skofl/<release> layout. The rootcint call needs the cint stub
    include dirs (cint/cint/{include,stl,lib}) to run on the container's ROOT install --
    same fix as the SKG4 GNUmakefile rootcint recipe (see SKG4-Patches.txt).