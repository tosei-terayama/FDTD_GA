#include "fdtd3d.h"
#include "pml.h"

void PML_idx_initialize(
    pml* idx_Dr, pml* idx_Dth, pml* idx_Dphi,
    pml* idx_Hr, pml* idx_Hth, pml* idx_Hphi)
{
    //define PML index of D field//
    idx_Dr[0].set_point(1, L, 1, Nphi - 1);
    idx_Dr[1].set_point(Ntheta - L, Ntheta - 1, 1, Nphi - 1);
    idx_Dr[2].set_point(L + 1, Ntheta - L - 1, 1, L);
    idx_Dr[3].set_point(L + 1, Ntheta - L - 1, Nphi - L, Nphi - 1);

    idx_Dth[0].set_point(0, L - 1, 1, Nphi - 1);
    idx_Dth[1].set_point(Ntheta - L, Ntheta - 1, 1, Nphi - 1);
    idx_Dth[2].set_point(L, Ntheta -L - 1, 1, L);
    idx_Dth[3].set_point(L, Ntheta - L - 1, Nphi - L, Nphi - 1);

    idx_Dphi[0].set_point(1, L, 0, Nphi - 1);
    idx_Dphi[1].set_point(Ntheta - L, Ntheta - 1, 0, Nphi - 1);
    idx_Dphi[2].set_point(L + 1, Ntheta - L - 1, 0, L - 1);
    idx_Dphi[3].set_point(L + 1, Ntheta - L - 1, Nphi - L, Nphi - 1);

    //define PML index of J field//
    idx_Hr[0].set_point(0, L - 1, 0, Nphi - 1);
    idx_Hr[1].set_point(Ntheta - L, Ntheta - 1, 0, Nphi - 1);
    idx_Hr[2].set_point(L, Ntheta - L - 1, 0, L - 1);
    idx_Hr[3].set_point(L, Ntheta - L - 1, Nphi - L, Nphi - 1);

    idx_Hth[0].set_point(1, L, 0, Nphi - 1);
    idx_Hth[1].set_point(Ntheta - L, Ntheta - 1, 0, Nphi - 1);
    idx_Hth[2].set_point(L + 1, Ntheta - L - 1, 0, L - 1);
    idx_Hth[3].set_point(L + 1, Ntheta - L - 1, Nphi - L, Nphi - 1);

    idx_Hphi[0].set_point(0, L - 1, 1, Nphi - 1);
    idx_Hphi[1].set_point(Ntheta - L, Ntheta - 1, 1, Nphi - 1);
    idx_Hphi[2].set_point(L, Ntheta - L - 1, 1, L);
    idx_Hphi[3].set_point(L, Ntheta - L - 1, Nphi - L, Nphi - 1);
}