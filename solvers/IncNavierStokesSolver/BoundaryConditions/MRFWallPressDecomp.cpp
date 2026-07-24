///////////////////////////////////////////////////////////////////////////////
//
// File: MRFWallPressDecomp.cpp
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
// Description: Wall boundary condition of moving reference frame with pressure
// decomposition.
//
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/BoundaryConditions/MRFWallPressDecomp.h>
#include <LibUtilities/Communication/Comm.h>

namespace Nektar
{
std::string MRFWallPressDecomp::className =
    GetIncBCFactory().RegisterCreatorFunction(
        "MRFWallPressDecomp", MRFWallPressDecomp::create,
        "Wall boundary condition of moving reference frame with pressure "
        "decomposition");

MRFWallPressDecomp::MRFWallPressDecomp(
    [[maybe_unused]] const LibUtilities::SessionReaderSharedPtr pSession,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> pFields,
    [[maybe_unused]] Array<OneD, SpatialDomains::BoundaryConditionShPtr> cond,
    [[maybe_unused]] Array<OneD, MultiRegions::ExpListSharedPtr> exp,
    [[maybe_unused]] int nbnd, [[maybe_unused]] int spacedim,
    [[maybe_unused]] int bnddim)
    : StaticWall(pSession, pFields, cond, exp, nbnd, spacedim, bnddim)
{
    classname = "MRFWallPressDecomp";
    m_hasVels = false;
    m_velocityBasisReady = false;
    m_velocityBasisNpts0 = 0;
    for (size_t i = 0; i < m_bnddim; ++i)
    {
        if (cond[i]->GetUserDefined() == classname)
        {
            m_hasVels   = true;
            m_BndExp[i] = exp[i];
        }
    }
    m_pressure    = cond.size() - 1;
    m_hasPressure = cond[m_pressure]->GetUserDefined() == classname;
    if (m_hasPressure)
    {
        m_BndExp[m_pressure] = exp[m_pressure];
    }
}

MRFWallPressDecomp::~MRFWallPressDecomp()
{
}

void MRFWallPressDecomp::v_Initialise(
    const LibUtilities::SessionReaderSharedPtr &pSession)
{
    IncBaseCondition::v_Initialise(pSession);
    if (!(pSession->GetSolverInfo("SolverType") == "PressDecompVCSFSI"))
    {
        ASSERTL0(false,
                 "The boundary condition MRFWallPressDecomp is only "
                 "supported for use in the SolverType 'PressDecompVCSFSI'");
    }
    m_field->GetBndElmtExpansion(m_nbnd, m_bndElmtExps, false);
    if (m_hasPressure)
    {
        m_extrapArray =
            Array<OneD, Array<OneD, Array<OneD, NekDouble>>>(m_intSteps);
        for (int n = 0; n < m_intSteps; ++n)
        {
            m_extrapArray[n] = Array<OneD, Array<OneD, NekDouble>>(m_bnddim);
            for (int i = 0; i < m_bnddim; ++i)
            {
                m_extrapArray[n][i] = Array<OneD, NekDouble>(m_npoints, 0.0);
            }
        }
    }
}

void MRFWallPressDecomp::EnsurePressureRhs()
{
    if (m_pressureRhs.size() != m_bnddim)
    {
        m_pressureRhs = Array<OneD, Array<OneD, NekDouble>>(m_bnddim);
    }
    for (int i = 0; i < m_bnddim; ++i)
    {
        if (m_pressureRhs[i].size() != m_npoints)
        {
            m_pressureRhs[i] = Array<OneD, NekDouble>(m_npoints, 0.0);
        }
        else
        {
            Vmath::Zero(m_npoints, m_pressureRhs[i], 1);
        }
    }
}

/// @brief v_Update set correct BCs (in wavespace)
/// @param fields
/// @param Adv is in wavespace for 3DH1D
/// @param params
void MRFWallPressDecomp::v_Update(
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &fields,
    [[maybe_unused]] const Array<OneD, const Array<OneD, NekDouble>> &Adv,
    std::map<std::string, NekDouble> &params)
{
    int nptsPlane0 = 0;
    SetNumPointsOnPlane0(nptsPlane0);
    if (params.find("Omega_x") != params.end() ||
        params.find("Omega_y") != params.end() ||
        params.find("Omega_z") != params.end())
    {
        InitialiseCoords(params);
    }
    // pressure, do if not define velocity
    if (params.find("pressure") != params.end() && m_hasPressure &&
        fields.size() > 0)
    {
        ++m_numCalls;

        EnsurePressureRhs();
        // add viscous term and centripetal acceleration
        AddExtrapCentVisPressureBCs(fields, m_pressureRhs, params, nptsPlane0);
        m_BndExp[m_pressure]->NormVectorIProductWRTBase(
            m_pressureRhs, m_BndExp[m_pressure]->UpdateCoeffs());
    }
    // velocity, do if define velocity
    if (params.find("velocity") != params.end() && m_hasVels && nptsPlane0)
    {
        EnsureVelocityCoeffBasis(params, nptsPlane0);
        ApplyVelocityCoeffBasis(params);
    }
}

int MRFWallPressDecomp::GetVelocityBasisDofCount() const
{
    return m_bnddim == 3 ? 6 : 3;
}

int MRFWallPressDecomp::GetVelocityBasisIndex(int component, int dof) const
{
    return component * GetVelocityBasisDofCount() + dof;
}

std::string MRFWallPressDecomp::GetVelocityBasisParamName(int dof) const
{
    if (m_bnddim == 3)
    {
        switch (dof)
        {
            case 0:
                return "U";
            case 1:
                return "V";
            case 2:
                return "W";
            case 3:
                return "Omega_x";
            case 4:
                return "Omega_y";
            case 5:
                return "Omega_z";
            default:
                ASSERTL0(false, "Invalid 3D velocity basis dof");
        }
    }
    else
    {
        switch (dof)
        {
            case 0:
                return "U";
            case 1:
                return "V";
            case 2:
                return "Omega_z";
            default:
                ASSERTL0(false, "Invalid 2D velocity basis dof");
        }
    }
    return "";
}

NekDouble MRFWallPressDecomp::GetVelocityBasisParamValue(
    int dof, const std::map<std::string, NekDouble> &params) const
{
    const std::string name = GetVelocityBasisParamName(dof);
    auto it                = params.find(name);
    return it == params.end() ? 0.0 : it->second;
}

MultiRegions::ExpListSharedPtr MRFWallPressDecomp::GetVelocityBoundaryExpansion(
    int component) const
{
    auto it = m_BndExp.find(component);
    if (it == m_BndExp.end())
    {
        return MultiRegions::ExpListSharedPtr();
    }

    if (it->second->GetExpType() == MultiRegions::e2DH1D)
    {
        return it->second->GetPlane(0);
    }
    return it->second;
}

void MRFWallPressDecomp::EnsureVelocityCoeffBasis(
    std::map<std::string, NekDouble> &params, int npts0)
{
    ASSERTL0(npts0 > 0, "Velocity boundary basis requires non-zero points");

    if (m_velocityBasisReady && m_velocityBasisNpts0 == npts0)
    {
        bool valid = true;
        for (int k = 0; k < m_bnddim; ++k)
        {
            MultiRegions::ExpListSharedPtr exp =
                GetVelocityBoundaryExpansion(k);
            if (!exp)
            {
                continue;
            }

            for (int dof = 0; dof < GetVelocityBasisDofCount(); ++dof)
            {
                const int idx = GetVelocityBasisIndex(k, dof);
                if (idx >= m_velocityCoeffBasis.size() ||
                    m_velocityCoeffBasis[idx].size() != exp->GetNcoeffs())
                {
                    valid = false;
                    break;
                }
            }
        }
        if (valid)
        {
            return;
        }
    }

    InitialiseCoords(params);

    const int nDofs = GetVelocityBasisDofCount();
    m_velocityCoeffBasis =
        Array<OneD, Array<OneD, NekDouble>>(m_bnddim * nDofs);
    m_velocityBasisNpts0 = npts0;

    for (int k = 0; k < m_bnddim; ++k)
    {
        MultiRegions::ExpListSharedPtr exp = GetVelocityBoundaryExpansion(k);
        if (!exp)
        {
            continue;
        }

        const int nCoeffs = exp->GetNcoeffs();
        for (int dof = 0; dof < nDofs; ++dof)
        {
            m_velocityCoeffBasis[GetVelocityBasisIndex(k, dof)] =
                Array<OneD, NekDouble>(nCoeffs, 0.0);
        }
    }

    Array<OneD, Array<OneD, NekDouble>> velocities(m_bnddim);
    for (int k = 0; k < m_bnddim; ++k)
    {
        if (m_BndExp.find(k) != m_BndExp.end())
        {
            velocities[k] = Array<OneD, NekDouble>(npts0, 0.0);
        }
    }

    for (int dof = 0; dof < nDofs; ++dof)
    {
        for (int k = 0; k < m_bnddim; ++k)
        {
            if (velocities[k].size() > 0)
            {
                Vmath::Zero(npts0, velocities[k], 1);
            }
        }

        std::map<std::string, NekDouble> unitParams;
        unitParams[GetVelocityBasisParamName(dof)] = 1.0;
        RigidBodyVelocity(velocities, unitParams, npts0);

        for (int k = 0; k < m_bnddim; ++k)
        {
            MultiRegions::ExpListSharedPtr exp =
                GetVelocityBoundaryExpansion(k);
            if (!exp)
            {
                continue;
            }

            exp->FwdTransBndConstrained(
                velocities[k],
                m_velocityCoeffBasis[GetVelocityBasisIndex(k, dof)]);
        }

    }

    m_velocityBasisReady = true;
}

void MRFWallPressDecomp::ApplyVelocityCoeffBasis(
    std::map<std::string, NekDouble> &params)
{
    ASSERTL0(m_velocityBasisReady,
             "Velocity boundary basis must be built before applying it");

    const int nDofs = GetVelocityBasisDofCount();
    for (int k = 0; k < m_bnddim; ++k)
    {
        MultiRegions::ExpListSharedPtr exp = GetVelocityBoundaryExpansion(k);
        if (!exp)
        {
            continue;
        }

        Array<OneD, NekDouble> coeffs = exp->UpdateCoeffs();
        const int nCoeffs             = coeffs.size();
        Vmath::Zero(nCoeffs, coeffs, 1);

        for (int dof = 0; dof < nDofs; ++dof)
        {
            const NekDouble value = GetVelocityBasisParamValue(dof, params);
            if (value != 0.0)
            {
                Vmath::Svtvp(
                    nCoeffs, value,
                    m_velocityCoeffBasis[GetVelocityBasisIndex(k, dof)], 1,
                    coeffs, 1, coeffs, 1);
            }
        }
    }
}

void MRFWallPressDecomp::AddCentripetalAcc(
    Array<OneD, Array<OneD, NekDouble>> &N,
    std::map<std::string, NekDouble> &params, int npts0)
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
    const NekDouble Wx = getParam("Omega_x");
    const NekDouble Wy = getParam("Omega_y");
    const NekDouble Wz = getParam("Omega_z");

    if (Wx == 0.0 && Wy == 0.0 && Wz == 0.0)
    {
        return;
    }

    if (m_bnddim > 0)
    {
        Vmath::Svtvp(npts0, Wy * Wy + Wz * Wz, m_coords[0], 1, N[0], 1,
                     N[0], 1);
        if (m_spacedim > 1)
        {
            Vmath::Svtvp(npts0, -Wx * Wy, m_coords[1], 1, N[0], 1, N[0],
                         1);
        }
        if (m_spacedim > 2)
        {
            Vmath::Svtvp(npts0, -Wx * Wz, m_coords[2], 1, N[0], 1, N[0],
                         1);
        }
        Vmath::Sadd(npts0, -Wy * w0 + Wz * v0, N[0], 1, N[0], 1);
    }

    if (m_bnddim > 1)
    {
        Vmath::Svtvp(npts0, -Wx * Wy, m_coords[0], 1, N[1], 1, N[1], 1);
        if (m_spacedim > 1)
        {
            Vmath::Svtvp(npts0, Wx * Wx + Wz * Wz, m_coords[1], 1, N[1],
                         1, N[1], 1);
        }
        if (m_spacedim > 2)
        {
            Vmath::Svtvp(npts0, -Wy * Wz, m_coords[2], 1, N[1], 1, N[1],
                         1);
        }
        Vmath::Sadd(npts0, -Wz * u0 + Wx * w0, N[1], 1, N[1], 1);
    }

    if (m_bnddim > 2)
    {
        Vmath::Svtvp(npts0, -Wx * Wz, m_coords[0], 1, N[2], 1, N[2], 1);
        Vmath::Svtvp(npts0, -Wy * Wz, m_coords[1], 1, N[2], 1, N[2], 1);
        Vmath::Svtvp(npts0, Wx * Wx + Wy * Wy, m_coords[2], 1, N[2], 1,
                     N[2], 1);
        Vmath::Sadd(npts0, -Wx * v0 + Wy * u0, N[2], 1, N[2], 1);
    }
}

void MRFWallPressDecomp::AddExtrapCentVisPressureBCs(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    Array<OneD, Array<OneD, NekDouble>> &N,
    std::map<std::string, NekDouble> &params, int npts0)
{
    for (int i = 0; i < m_bnddim; ++i)
    {
        Vmath::Zero(m_npoints, m_extrapArray[m_intSteps - 1][i], 1);
    }
    AddVisPressureBCs(fields, m_extrapArray[m_intSteps - 1], params);
    AddCentripetalAcc(m_extrapArray[m_intSteps - 1], params, npts0);
    ExtrapolateArray(m_numCalls, m_extrapArray);
    for (int i = 0; i < m_bnddim; i++)
    {
        Vmath::Vadd(m_npoints, m_extrapArray[m_intSteps - 1][i], 1, N[i], 1,
                    N[i], 1);
    }
}

} // namespace Nektar
