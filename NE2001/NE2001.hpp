///Wrappers for the NE2001 subroutines
#pragma once

#include<fortran.h>

extern "C" {
    extern void FORTRAN_NAME(dmdsm)(
    ///All arguments must exists and be passed as a pointer to them.
    ///Input
    float * p_l,       /// galactic longitude in radians 
    float * p_b,       /// galactic latitude in radians
    int * p_ndir,      /// >= 0 calculates dist from dmpsr, < 0 for dmpsr from dist

    ///Input or output
    float * p_dmpsr,   /// (dispersion measure in pc/cm^3)
    float * p_dist,    /// (distance in kpc)

    ///Output
    char * p_limit,    /// (set to '>' if only a lower distance limit can be given; otherwise set to ' ')
    int * size,        /// size of p_limit char array (compability with fortran)
    float * p_sm,      /// (scattering measure, uniform weighting) (kpc/m^{20/3})
    float * p_smtau,   /// (scattering measure, weighting for pulse broadening)
    float * p_smtheta, /// (scattering measure, weighting for angular broadening of galactic sources)
    float * p_smiso);  /// (scattering measure appropriate for calculating the isoplanatic angle at the source's location
}

///Clike wrapper for the dmdsm function
///Can'nt use const arg practice - not known effect in fortran's subrouitne call
float fort_dmdsm(float l, float b, int ndir, float dmpsr, float dist);

