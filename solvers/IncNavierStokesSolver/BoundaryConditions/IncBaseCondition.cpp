///////////////////////////////////////////////////////////////////////////////
//
// File: IncBaseCondition.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: Abstract base class for Extrapolate.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/IncBaseCondition.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{

NekDouble IncBaseCondition::StifflyStable_Betaq_Coeffs[3][3] = {
    {1.0, 0.0, 0.0}, {2.0, -1.0, 0.0}, {3.0, -3.0, 1.0}};
NekDouble IncBaseCondition::StifflyStable_Alpha_Coeffs[3][3] = {
    {1.0, 0.0, 0.0}, {2.0, -0.5, 0.0}, {3.0, -1.5, 1.0 / 3.0}};
NekDouble IncBaseCondition::StifflyStable_Gamma0_Coeffs[3] = {1.0, 1.5,
                                                              11.0 / 6.0};

IncBCFactory &GetIncBCFactory()
{
    static IncBCFactory instance;
    return instance;
}

IncBaseCondition::IncBaseCondition(
    [[maybe_unused]] const LibUtilities::SessionReaderSharedPtr pSession,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    [[maybe_unused]] Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> exp,
    [[maybe_unused]] int nbnd, [[maybe_unused]] int spacedim,
    [[maybe_unused]] int bnddim)
    : m_spacedim(spacedim), m_bnddim(bnddim), m_nbnd(nbnd), m_field(pFields[0])
{
}

void IncBaseCondition::v_Initialise(
    const LibUtilities::SessionReaderSharedPtr &pSession)
{
    m_npoints  = m_BndExp.begin()->second->GetNpoints();
    m_numCalls = 0;
    if (pSession->DefinesParameter("ExtrapolateOrder"))
    {
        m_intSteps = std::round(pSession->GetParameter("ExtrapolateOrder"));
        m_intSteps = std::min(3, std::max(0, m_intSteps));
    }
    else if (pSession->DefinesTimeIntScheme())
    {
        m_intSteps = pSession->GetTimeIntScheme().order;
    }
    else
    {
        m_intSteps = 0;
    }
}

void IncBaseCondition::SetNumPointsOnPlane0(int &npointsPlane0)
{
    if (m_BndExp.begin()->second->GetExpType() == MultiRegions::e2DH1D)
    {
        if (m_field->GetZIDs()[0] == 0)
        {
            npointsPlane0 = m_BndExp.begin()->second->GetPlane(0)->GetNpoints();
        }
        else
        {
            npointsPlane0 = 0;
        }
    }
    else
    {
        npointsPlane0 = m_BndExp.begin()->second->GetNpoints();
    }
}

void IncBaseCondition::InitialiseCoords(
    std::map<std::string, NekDouble> &params)
{
    MultiRegions::ExpListSharedPtr bndexp = m_BndExp.begin()->second;
    if (m_coords.size() == 0)
    {
        m_coords = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        for (size_t k = 0; k < m_spacedim; ++k)
        {
            m_coords[k] = Array<OneD, NekDouble>(m_npoints, 0.0);
        }
        if (m_spacedim == 2)
        {
            bndexp->GetCoords(m_coords[0], m_coords[1]);
        }
        else
        {
            bndexp->GetCoords(m_coords[0], m_coords[1], m_coords[2]);
        }
        // move the centre to the location of pivot
        std::vector<std::string> xyz = {"X0", "Y0", "Z0"};
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (params.find(xyz[i]) != params.end())
            {
                Vmath::Sadd(m_npoints, -params[xyz[i]], m_coords[i], 1,
                            m_coords[i], 1);
            }
        }
    }
}

void IncBaseCondition::ExtrapolateArray(
    const int numCalls, Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &array)
{
    if (m_intSteps == 1)
    {
        return;
    }
    int nint    = std::min(numCalls, m_intSteps);
    int nlevels = array.size();
    int dim     = array[0].size();
    int nPts    = array[0][0].size();
    // Check integer for time levels
    // Note that ExtrapolateArray assumes m_pressureCalls is >= 1
    // meaning v_EvaluatePressureBCs has been called previously
    ASSERTL0(nint > 0, "nint must be > 0 when calling ExtrapolateArray.");
    // Update array
    RollOver(array);
    // Extrapolate to outarray
    for (int i = 0; i < dim; ++i)
    {
        Vmath::Smul(nPts, StifflyStable_Betaq_Coeffs[nint - 1][nint - 1],
                    array[nint - 1][i], 1, array[nlevels - 1][i], 1);
    }
    for (int n = 0; n < nint - 1; ++n)
    {
        for (int i = 0; i < dim; ++i)
        {
            Vmath::Svtvp(nPts, StifflyStable_Betaq_Coeffs[nint - 1][n],
                         array[n][i], 1, array[nlevels - 1][i], 1,
                         array[nlevels - 1][i], 1);
        }
    }
}

void IncBaseCondition::AddRigidBodyAcc(Array<OneD, Array<OneD, NekDouble>> &N,
                                       std::map<std::string, NekDouble> &params,
                                       int npts0)
{
    if (npts0 == 0)
    {
        return;
    }
    InitialiseCoords(params);

    const auto getParam = [&](const std::string &name) {
        auto it = params.find(name);
        return it == params.end() ? 0.0 : it->second;
    };

    const NekDouble u0 = getParam("U");
    const NekDouble v0 = getParam("V");
    const NekDouble w0 = getParam("W");
    const NekDouble Ax = getParam("A_x");
    const NekDouble Ay = getParam("A_y");
    const NekDouble Az = getParam("A_z");

    const NekDouble Wx = getParam("Omega_x");
    const NekDouble Wy = getParam("Omega_y");
    const NekDouble Wz = getParam("Omega_z");
    const NekDouble DWx = getParam("DOmega_x");
    const NekDouble DWy = getParam("DOmega_y");
    const NekDouble DWz = getParam("DOmega_z");

    if (m_bnddim > 0)
    {
        Vmath::Svtvp(npts0, Wy * Wy + Wz * Wz, m_coords[0], 1, N[0], 1,
                     N[0], 1);
        if (m_spacedim > 1)
        {
            Vmath::Svtvp(npts0, -Wx * Wy + DWz, m_coords[1], 1, N[0], 1,
                         N[0], 1);
        }
        if (m_spacedim > 2)
        {
            Vmath::Svtvp(npts0, -Wx * Wz - DWy, m_coords[2], 1, N[0], 1,
                         N[0], 1);
        }
        Vmath::Sadd(npts0, -Ax - Wy * w0 + Wz * v0, N[0], 1, N[0], 1);
    }

    if (m_bnddim > 1)
    {
        Vmath::Svtvp(npts0, -Wx * Wy - DWz, m_coords[0], 1, N[1], 1, N[1],
                     1);
        if (m_spacedim > 1)
        {
            Vmath::Svtvp(npts0, Wx * Wx + Wz * Wz, m_coords[1], 1, N[1], 1,
                         N[1], 1);
        }
        if (m_spacedim > 2)
        {
            Vmath::Svtvp(npts0, -Wy * Wz + DWx, m_coords[2], 1, N[1], 1,
                         N[1], 1);
        }
        Vmath::Sadd(npts0, -Ay - Wz * u0 + Wx * w0, N[1], 1, N[1], 1);
    }

    if (m_bnddim > 2)
    {
        Vmath::Svtvp(npts0, -Wx * Wz + DWy, m_coords[0], 1, N[2], 1, N[2],
                     1);
        Vmath::Svtvp(npts0, -Wy * Wz - DWx, m_coords[1], 1, N[2], 1, N[2],
                     1);
        Vmath::Svtvp(npts0, Wx * Wx + Wy * Wy, m_coords[2], 1, N[2], 1,
                     N[2], 1);
        Vmath::Sadd(npts0, -Az - Wx * v0 + Wy * u0, N[2], 1, N[2], 1);
    }
}

void IncBaseCondition::EnsureVisPressureScratch(const int nq)
{
    if (m_visPressureVelocity.size() != m_spacedim)
    {
        m_visPressureVelocity =
            Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
        m_visPressureQ = Array<OneD, Array<OneD, NekDouble>>(m_spacedim);
    }
    for (int i = 0; i < m_spacedim; ++i)
    {
        if (m_visPressureVelocity[i].size() != nq)
        {
            m_visPressureVelocity[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        if (m_visPressureQ[i].size() != nq)
        {
            m_visPressureQ[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
        else
        {
            Vmath::Zero(nq, m_visPressureQ[i], 1);
        }
    }

    const int nTmp = m_spacedim == 2 ? 3 : 7;
    if (m_visPressureCurlTmp.size() != nTmp)
    {
        m_visPressureCurlTmp =
            Array<OneD, Array<OneD, NekDouble>>(nTmp);
    }
    for (int i = 0; i < nTmp; ++i)
    {
        if (m_visPressureCurlTmp[i].size() != nq)
        {
            m_visPressureCurlTmp[i] = Array<OneD, NekDouble>(nq, 0.0);
        }
    }

    if (m_visPressureTemp.size() != m_npoints)
    {
        m_visPressureTemp = Array<OneD, NekDouble>(m_npoints, 0.0);
    }
}

void IncBaseCondition::CurlCurlWithScratch(
    Array<OneD, Array<OneD, NekDouble>> &vel,
    Array<OneD, Array<OneD, NekDouble>> &q)
{
    int nq = m_bndElmtExps->GetTotPoints();

    Array<OneD, NekDouble> Vx    = m_visPressureCurlTmp[0];
    Array<OneD, NekDouble> Uy    = m_visPressureCurlTmp[1];
    Array<OneD, NekDouble> Dummy = m_visPressureCurlTmp[2];

    bool halfMode = false;
    if (m_bndElmtExps->GetExpType() == MultiRegions::e3DH1D)
    {
        m_field->GetSession()->MatchSolverInfo("ModeType", "HalfMode",
                                               halfMode, false);
    }

    switch (m_bndElmtExps->GetExpType())
    {
        case MultiRegions::e2D:
        {
            m_bndElmtExps->PhysDeriv(MultiRegions::eX, vel[1], Vx);
            m_bndElmtExps->PhysDeriv(MultiRegions::eY, vel[0], Uy);

            Vmath::Vsub(nq, Vx, 1, Uy, 1, Dummy, 1);

            m_bndElmtExps->PhysDeriv(Dummy, q[1], q[0]);

            Vmath::Smul(nq, -1.0, q[1], 1, q[1], 1);
        }
        break;

        case MultiRegions::e3D:
        case MultiRegions::e3DH1D:
        case MultiRegions::e3DH2D:
        {
            Array<OneD, NekDouble> Vz = m_visPressureCurlTmp[3];
            Array<OneD, NekDouble> Uz = m_visPressureCurlTmp[4];
            Array<OneD, NekDouble> Wx = m_visPressureCurlTmp[5];
            Array<OneD, NekDouble> Wy = m_visPressureCurlTmp[6];

            m_bndElmtExps->PhysDeriv(vel[0], Dummy, Uy, Uz);
            m_bndElmtExps->PhysDeriv(vel[1], Vx, Dummy, Vz);
            m_bndElmtExps->PhysDeriv(vel[2], Wx, Wy, Dummy);

            Vmath::Vsub(nq, Wy, 1, Vz, 1, q[0], 1);
            Vmath::Vsub(nq, Uz, 1, Wx, 1, q[1], 1);
            Vmath::Vsub(nq, Vx, 1, Uy, 1, q[2], 1);

            m_bndElmtExps->PhysDeriv(q[0], Dummy, Uy, Uz);
            m_bndElmtExps->PhysDeriv(q[1], Vx, Dummy, Vz);
            m_bndElmtExps->PhysDeriv(q[2], Wx, Wy, Dummy);

            if (halfMode)
            {
                Vmath::Neg(nq, Uz, 1);
                Vmath::Neg(nq, Vz, 1);
            }

            Vmath::Vsub(nq, Wy, 1, Vz, 1, q[0], 1);
            Vmath::Vsub(nq, Uz, 1, Wx, 1, q[1], 1);
            Vmath::Vsub(nq, Vx, 1, Uy, 1, q[2], 1);
        }
        break;
        default:
            ASSERTL0(false, "Dimension not supported");
            break;
    }
}

void IncBaseCondition::AddVisPressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    Array<OneD, Array<OneD, NekDouble>> &N,
    std::map<std::string, NekDouble> &params)
{
    if (m_intSteps == 0 || params.find("Kinvis") == params.end() ||
        params["Kinvis"] <= 0. || fields.size() == 0)
    {
        return;
    }
    NekDouble kinvis = params["Kinvis"];
    m_bndElmtExps->SetWaveSpace(m_field->GetWaveSpace());
    int nq = m_bndElmtExps->GetTotPoints();
    EnsureVisPressureScratch(nq);

    for (int i = 0; i < m_spacedim; i++)
    {
        m_field->ExtractPhysToBndElmt(m_nbnd, fields[i],
                                      m_visPressureVelocity[i]);
    }

    CurlCurlWithScratch(m_visPressureVelocity, m_visPressureQ);

    for (int i = 0; i < m_bnddim; i++)
    {
        m_field->ExtractElmtToBndPhys(m_nbnd, m_visPressureQ[i],
                                      m_visPressureTemp);
        Vmath::Svtvp(m_npoints, -kinvis, m_visPressureTemp, 1, N[i], 1, N[i],
                     1);
    }
}

void IncBaseCondition::AddExtrapVisPressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    Array<OneD, Array<OneD, NekDouble>> &N,
    std::map<std::string, NekDouble> &params)
{
    for (int i = 0; i < m_bnddim; ++i)
    {
        Vmath::Zero(m_npoints, m_extrapArray[m_intSteps - 1][i], 1);
    }
    AddVisPressureBCs(fields, m_extrapArray[m_intSteps - 1], params);
    ExtrapolateArray(m_numCalls, m_extrapArray);
    for (int i = 0; i < m_bnddim; i++)
    {
        Vmath::Vadd(m_npoints, m_extrapArray[m_intSteps - 1][i], 1, N[i], 1,
                    N[i], 1);
    }
}

void IncBaseCondition::RollOver(
    Array<OneD, Array<OneD, Array<OneD, NekDouble>>> &input)
{
    int nlevels = input.size();
    Array<OneD, Array<OneD, NekDouble>> tmp;
    tmp = input[nlevels - 1];
    for (int n = nlevels - 1; n > 0; --n)
    {
        input[n] = input[n - 1];
    }
    input[0] = tmp;
}

void IncBaseCondition::RigidBodyVelocity(
    Array<OneD, Array<OneD, NekDouble>> &velocities,
    std::map<std::string, NekDouble> &params, int npts0)
{
    if (npts0 == 0)
    {
        return;
    }
    // for the wall we need to calculate:
    // [V_wall]_xyz = [V_frame]_xyz + [Omega X r]_xyz
    // Note all vectors must be in moving frame coordinates xyz
    // not in inertial frame XYZ
    InitialiseCoords(params);

    // vx = OmegaY*z-OmegaZ*y
    // vy = OmegaZ*x-OmegaX*z
    // vz = OmegaX*y-OmegaY*x
    if (params.find("Omega_z") != params.end())
    {
        NekDouble Wz = params["Omega_z"];
        if (m_BndExp.find(0) != m_BndExp.end())
        {
            Vmath::Smul(npts0, -Wz, m_coords[1], 1, velocities[0], 1);
        }
        if (m_BndExp.find(1) != m_BndExp.end())
        {
            Vmath::Smul(npts0, Wz, m_coords[0], 1, velocities[1], 1);
        }
    }
    if (m_bnddim == 3)
    {
        if (params.find("Omega_x") != params.end())
        {
            NekDouble Wx = params["Omega_x"];
            if (m_BndExp.find(2) != m_BndExp.end())
            {
                Vmath::Smul(npts0, Wx, m_coords[1], 1, velocities[2], 1);
            }
            if (m_BndExp.find(1) != m_BndExp.end())
            {
                Vmath::Svtvp(npts0, -Wx, m_coords[2], 1, velocities[1], 1,
                             velocities[1], 1);
            }
        }
        if (params.find("Omega_y") != params.end())
        {
            NekDouble Wy = params["Omega_y"];
            if (m_BndExp.find(0) != m_BndExp.end())
            {
                Vmath::Svtvp(npts0, Wy, m_coords[2], 1, velocities[0], 1,
                             velocities[0], 1);
            }
            if (m_BndExp.find(2) != m_BndExp.end())
            {
                Vmath::Svtvp(npts0, -Wy, m_coords[0], 1, velocities[2], 1,
                             velocities[2], 1);
            }
        }
    }

    // add the translation velocity
    std::vector<std::string> vars = {"U", "V", "W"};
    for (int k = 0; k < m_bnddim; ++k)
    {
        if (params.find(vars[k]) != params.end() &&
            m_BndExp.find(k) != m_BndExp.end())
        {
            Vmath::Sadd(npts0, params[vars[k]], velocities[k], 1, velocities[k],
                        1);
        }
    }
}

} // namespace Nektar
