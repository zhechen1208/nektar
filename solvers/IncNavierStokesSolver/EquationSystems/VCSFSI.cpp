///////////////////////////////////////////////////////////////////////////////
//
// File: VCSFSI.cpp
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
// Description: Velocity Correction Scheme for fluid-structure interaction
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/VCSFSI.h>
#include <LibUtilities/BasicUtils/Timer.h>
#include <MultiRegions/ContField.h>
#include <SolverUtils/Core/Misc.h>

#include <boost/algorithm/string.hpp>

namespace Nektar
{
std::string VCSFSI::className =
    SolverUtils::GetEquationSystemFactory().RegisterCreatorFunction(
        "VCSFSI", VCSFSI::create);

std::string VCSFSI::solverTypeLookupId =
    LibUtilities::SessionReader::RegisterEnumValue("SolverType", "VCSFSI",
                                                   eVCSFSI);

/**
 * Constructor. Creates ...
 *
 * \param
 * \param
 */
VCSFSI::VCSFSI(const LibUtilities::SessionReaderSharedPtr &pSession,
               const SpatialDomains::MeshGraphSharedPtr &pGraph)
    : UnsteadySystem(pSession, pGraph),
      VelocityCorrectionScheme(pSession, pGraph),
      m_enablePressureDecomposition(false),
      m_pressureDecompWriteFld(true),
      m_pressureDecompOutputFrequency(0),
      m_pressureDecompOutputIndex(0)
{
}

void VCSFSI::v_InitObject(bool DeclareField)
{
    VelocityCorrectionScheme::v_InitObject(DeclareField);
    Array<OneD, NekDouble> tmp = m_movingFrameData + 18;
    m_rigidSolver.InitObject(m_session, m_fields[0], tmp);
    m_rigidSolver.SetMovableDoFs(m_movableDoFs);
    InitialisePressureDecomposition();
}

/**
 * Destructor
 */
VCSFSI::~VCSFSI(void)
{
}

void VCSFSI::v_DoInitialise(bool dumpInitialConditions)
{
    m_rigidSolver.SetInitialConditions(m_session, m_movingFrameData);
    VelocityCorrectionScheme::v_DoInitialise(dumpInitialConditions);
    Array<OneD, NekDouble> AddedMass;
    m_rigidSolver.SetNewmarkBetaSolver(AddedMass);
    Array<OneD, NekDouble> aeroforce(12, 0.);
    InitialiseFilter(aeroforce);
    m_rigidSolver.SetOldFvis(aeroforce);
}

void VCSFSI::InitialisePressureDecomposition()
{
    m_session->MatchSolverInfo("PressureDecomposition", "True",
                               m_enablePressureDecomposition, false);
    m_session->MatchSolverInfo("PressureDecompWriteFld", "True",
                               m_pressureDecompWriteFld, true);
    m_session->LoadParameter("PressureDecompOutputFrequency",
                             m_pressureDecompOutputFrequency, 0);

    if (m_enablePressureDecomposition)
    {
        auto pressure = std::dynamic_pointer_cast<MultiRegions::ContField>(
            m_pressure);
        ASSERTL0(pressure,
                 "Pressure decomposition requires a continuous pressure "
                 "expansion.");

        // This constructor recreates boundary expansions for p while sharing
        // only the immutable mesh and assembly-map data with m_pressure.
        m_pressureDecomp =
            MemoryManager<MultiRegions::ContField>::AllocateSharedPtr(
                *pressure, m_graph,
                m_session->GetVariable(m_nConvectiveFields));

        m_paCoeff = Array<OneD, NekDouble>(m_pressure->GetNcoeffs(), 0.0);
        m_pqCoeff = Array<OneD, NekDouble>(m_pressure->GetNcoeffs(), 0.0);
        m_pvisCoeff =
            Array<OneD, NekDouble>(m_pressure->GetNcoeffs(), 0.0);
        m_paPhys  = Array<OneD, NekDouble>(m_pressure->GetTotPoints(), 0.0);
        m_pqPhys  = Array<OneD, NekDouble>(m_pressure->GetTotPoints(), 0.0);
        m_pvisPhys =
            Array<OneD, NekDouble>(m_pressure->GetTotPoints(), 0.0);
        m_paForce    = Array<OneD, NekDouble>(3, 0.0);
        m_pqForce    = Array<OneD, NekDouble>(3, 0.0);
        m_pvisForce  = Array<OneD, NekDouble>(3, 0.0);
        m_pForce     = Array<OneD, NekDouble>(3, 0.0);
        m_presForce  = Array<OneD, NekDouble>(3, 0.0);
        m_pressurePoissonRhs =
            Array<OneD, NekDouble>(m_pressure->GetTotPoints(), 0.0);
    }
}

void VCSFSI::v_SetUpPressureForcing(
    const Array<OneD, const Array<OneD, NekDouble>> &fields,
    Array<OneD, Array<OneD, NekDouble>> &Forcing, NekDouble aiiDt)
{
    VelocityCorrectionScheme::v_SetUpPressureForcing(fields, Forcing, aiiDt);

    if (m_enablePressureDecomposition)
    {
        if (m_pressurePoissonRhs.size() != Forcing[0].size())
        {
            m_pressurePoissonRhs =
                Array<OneD, NekDouble>(Forcing[0].size(), 0.0);
        }
        Vmath::Vcopy(Forcing[0].size(), Forcing[0], 1,
                     m_pressurePoissonRhs, 1);
    }
}

void VCSFSI::UpdatePressureDecomposition(NekDouble time)
{
    if (!m_enablePressureDecomposition)
    {
        return;
    }

    ++m_pressureDecompOutputIndex;
    ComputePaFull(time);
    ComputePq(time);
    ComputePvis(time);
    EvaluatePressureComponentForces(time);
    OutputPressureComponents(time);
}

void VCSFSI::CorrectPressureAfterSolid()
{
}

void VCSFSI::ZeroPressureBoundaryConditions()
{
    Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
        m_pressureDecomp->GetBndCondExpansions();

    for (int i = 0; i < bndExp.size(); ++i)
    {
        Vmath::Zero(bndExp[i]->GetNcoeffs(), bndExp[i]->UpdateCoeffs(), 1);
    }
}

void VCSFSI::ComputePaFull([[maybe_unused]] NekDouble time)
{
    auto isPaWall = [](const std::string &bcType) {
        return boost::iequals(bcType, "MRFWallPressDecomp") ||
               boost::iequals(bcType, "MRFWall");
    };

    Vmath::Zero(m_paCoeff.size(), m_paCoeff, 1);
    Vmath::Zero(m_paPhys.size(), m_paPhys, 1);

    // HelmSolve reads these coefficients directly. Each component starts from
    // homogeneous data and then supplies only its own boundary contribution.
    ZeroPressureBoundaryConditions();

    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds =
        m_pressureDecomp->GetBndConditions();
    Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
        m_pressureDecomp->GetBndCondExpansions();

    int numPts = 0;
    for (int i = 0; i < bndConds.size(); ++i)
    {
        if (isPaWall(bndConds[i]->GetUserDefined()))
        {
            numPts += bndExp[i]->GetTotPoints();
        }
    }

    int globalNumPts = numPts;
    m_session->GetComm()->AllReduce(globalNumPts, LibUtilities::ReduceSum);
    if (globalNumPts == 0)
    {
        return;
    }

    const NekDouble u0     = m_movingFrameData[6];
    const NekDouble v0     = m_movingFrameData[7];
    const NekDouble ax     = m_movingFrameData[12];
    const NekDouble ay     = m_movingFrameData[13];
    const NekDouble az     = m_movingFrameData[14];
    const NekDouble omega  = m_movingFrameData[11];
    const NekDouble dOmega = m_movingFrameData[17];

    if (std::abs(ax) < NekConstants::kNekZeroTol &&
        std::abs(ay) < NekConstants::kNekZeroTol &&
        std::abs(az) < NekConstants::kNekZeroTol &&
        std::abs(omega) < NekConstants::kNekZeroTol &&
        std::abs(u0) < NekConstants::kNekZeroTol &&
        std::abs(v0) < NekConstants::kNekZeroTol &&
        std::abs(dOmega) < NekConstants::kNekZeroTol)
    {
        return;
    }

    Array<OneD, NekDouble> bc(numPts, 0.0);
    const int ndim = m_fields[0]->GetCoordim(0);
    int offset = 0;
    for (int i = 0; i < bndConds.size(); ++i)
    {
        if (!isPaWall(bndConds[i]->GetUserDefined()))
        {
            continue;
        }

        const int npts = bndExp[i]->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble>> n(3);
        for (int j = 0; j < 3; ++j)
        {
            n[j] = Array<OneD, NekDouble>(npts, 0.0);
        }
        bndExp[i]->GetNormals(n);

        Array<OneD, NekDouble> bcSeg = bc + offset;
        if (ndim >= 2)
        {
            Array<OneD, NekDouble> x0(npts, 0.0), x1(npts, 0.0);
            Array<OneD, NekDouble> x2;
            if (ndim == 2)
            {
                bndExp[i]->GetCoords(x0, x1);
            }
            else
            {
                x2 = Array<OneD, NekDouble>(npts, 0.0);
                bndExp[i]->GetCoords(x0, x1, x2);
            }
            if (m_movingFrameData.size() >= 21)
            {
                Vmath::Sadd(npts, -m_movingFrameData[18], x0, 1, x0, 1);
                Vmath::Sadd(npts, -m_movingFrameData[19], x1, 1, x1, 1);
                if (ndim == 3)
                {
                    Vmath::Sadd(npts, -m_movingFrameData[20], x2, 1, x2, 1);
                }
            }

            NekDouble omega2 = omega * omega;
            Array<OneD, NekDouble> accX(npts, ax - omega * v0);
            Array<OneD, NekDouble> accY(npts, ay + omega * u0);

            if (std::abs(dOmega) >= NekConstants::kNekZeroTol)
            {
                Vmath::Svtvp(npts, -dOmega, x1, 1, accX, 1, accX, 1);
                Vmath::Svtvp(npts, dOmega, x0, 1, accY, 1, accY, 1);
            }
            if (std::abs(omega2) >= NekConstants::kNekZeroTol)
            {
                Vmath::Svtvp(npts, -omega2, x0, 1, accX, 1, accX, 1);
                Vmath::Svtvp(npts, -omega2, x1, 1, accY, 1, accY, 1);
            }

            Vmath::Vmul(npts, accX, 1, n[0], 1, bcSeg, 1);
            Array<OneD, NekDouble> tmp(npts, 0.0);
            Vmath::Vmul(npts, accY, 1, n[1], 1, tmp, 1);
            Vmath::Vadd(npts, tmp, 1, bcSeg, 1, bcSeg, 1);
        }
        else if (ndim == 1)
        {
            Vmath::Svtvp(npts, ax, n[0], 1, bcSeg, 1, bcSeg, 1);
        }

        if (ndim > 2)
        {
            Vmath::Svtvp(npts, az, n[2], 1, bcSeg, 1, bcSeg, 1);
        }

        Vmath::Smul(npts, -1.0, bcSeg, 1, bcSeg, 1);
        bndExp[i]->IProductWRTBase(bcSeg, bndExp[i]->UpdateCoeffs());
        offset += npts;
    }

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;
    Array<OneD, NekDouble> forcing(m_pressureDecomp->GetTotPoints(), 0.0);
    m_pressureDecomp->HelmSolve(forcing, m_paCoeff, factors);
    m_pressureDecomp->BwdTrans(m_paCoeff, m_paPhys);
}

void VCSFSI::ComputePq([[maybe_unused]] NekDouble time)
{
    Vmath::Zero(m_pqCoeff.size(), m_pqCoeff, 1);
    Vmath::Zero(m_pqPhys.size(), m_pqPhys, 1);

    ZeroPressureBoundaryConditions();

    const int npts = m_fields[0]->GetTotPoints();
    Array<OneD, NekDouble> forcing(npts, 0.0);
    Vmath::Vcopy(npts, m_pressurePoissonRhs, 1, forcing, 1);

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;
    m_pressureDecomp->HelmSolve(forcing, m_pqCoeff, factors);
    m_pressureDecomp->BwdTrans(m_pqCoeff, m_pqPhys);
}

void VCSFSI::ComputePvis([[maybe_unused]] NekDouble time)
{
    auto isPvisWall = [](const std::string &bcType) {
        return boost::iequals(bcType, "MRFWallPressDecomp") ||
               boost::iequals(bcType, "MRFWall");
    };

    Vmath::Zero(m_pvisCoeff.size(), m_pvisCoeff, 1);
    Vmath::Zero(m_pvisPhys.size(), m_pvisPhys, 1);

    if (m_kinvis <= 0.0)
    {
        return;
    }

    ZeroPressureBoundaryConditions();

    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds =
        m_pressureDecomp->GetBndConditions();
    Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
        m_pressureDecomp->GetBndCondExpansions();

    int numPts = 0;
    for (int i = 0; i < bndConds.size(); ++i)
    {
        if (isPvisWall(bndConds[i]->GetUserDefined()))
        {
            numPts += bndExp[i]->GetTotPoints();
        }
    }

    int globalNumPts = numPts;
    m_session->GetComm()->AllReduce(globalNumPts, LibUtilities::ReduceSum);
    if (globalNumPts == 0)
    {
        return;
    }

    const int ndim = m_velocity.size();
    Array<OneD, NekDouble> bc(numPts, 0.0);
    int offset = 0;
    for (int i = 0; i < bndConds.size(); ++i)
    {
        if (!isPvisWall(bndConds[i]->GetUserDefined()))
        {
            continue;
        }

        const int npts = bndExp[i]->GetTotPoints();
        MultiRegions::ExpListSharedPtr bndElmtExp;
        m_pressureDecomp->GetBndElmtExpansion(i, bndElmtExp, false);

        const int nq = bndElmtExp->GetTotPoints();
        Array<OneD, Array<OneD, NekDouble>> velocity(ndim), curlCurl(ndim);
        for (int j = 0; j < ndim; ++j)
        {
            velocity[j] = Array<OneD, NekDouble>(nq, 0.0);
            curlCurl[j] = Array<OneD, NekDouble>(nq, 0.0);
            m_pressureDecomp->ExtractPhysToBndElmt(
                i, m_fields[m_velocity[j]]->GetPhys(), velocity[j]);
        }

        bndElmtExp->SetWaveSpace(m_pressureDecomp->GetWaveSpace());
        bndElmtExp->CurlCurl(velocity, curlCurl);

        Array<OneD, Array<OneD, NekDouble>> n(3);
        for (int j = 0; j < 3; ++j)
        {
            n[j] = Array<OneD, NekDouble>(npts, 0.0);
        }
        bndExp[i]->GetNormals(n);

        Array<OneD, NekDouble> bcSeg = bc + offset;
        for (int j = 0; j < ndim; ++j)
        {
            Array<OneD, NekDouble> tmp(npts, 0.0);
            Array<OneD, NekDouble> viscComp(npts, 0.0);
            m_pressureDecomp->ExtractElmtToBndPhys(i, curlCurl[j], tmp);
            Vmath::Vmul(npts, tmp, 1, n[j], 1, viscComp, 1);
            Vmath::Svtvp(npts, -m_kinvis, viscComp, 1, bcSeg, 1, bcSeg, 1);
        }

        bndExp[i]->IProductWRTBase(bcSeg, bndExp[i]->UpdateCoeffs());
        offset += npts;
    }

    StdRegions::ConstFactorMap factors;
    factors[StdRegions::eFactorLambda] = 0.0;
    Array<OneD, NekDouble> forcing(m_pressureDecomp->GetTotPoints(), 0.0);
    m_pressureDecomp->HelmSolve(forcing, m_pvisCoeff, factors);
    m_pressureDecomp->BwdTrans(m_pvisCoeff, m_pvisPhys);
}

void VCSFSI::EvaluatePressureComponentForces([[maybe_unused]] NekDouble time)
{
    if (!m_pressureForceOutputInitialised)
    {
        InitialisePressureComponentForceOutput();
    }

    if (!m_pressureForceHasGlobalBoundary)
    {
        return;
    }

    IntegratePressureForce(m_paPhys, m_paForce);
    IntegratePressureForce(m_pqPhys, m_pqForce);
    IntegratePressureForce(m_pvisPhys, m_pvisForce);
    IntegratePressureForce(m_pressure->GetPhys(), m_pForce);

    // Surface integration is linear, so compute the residual force without
    // a fifth boundary traversal for the residual pressure field.
    Vmath::Vcopy(m_pForce.size(), m_pForce, 1, m_presForce, 1);
    Vmath::Vsub(m_presForce.size(), m_presForce, 1, m_paForce, 1,
                m_presForce, 1);
    Vmath::Vsub(m_presForce.size(), m_presForce, 1, m_pqForce, 1,
                m_presForce, 1);
    Vmath::Vsub(m_presForce.size(), m_presForce, 1, m_pvisForce, 1,
                m_presForce, 1);

    // FilterAeroForces reports forces in directions that rotate with the
    // moving frame. Apply the same z-axis rotation to every pressure force.
    if (m_movingFrameData.size() >= 6 && m_movingFrameData[5] != 0.0)
    {
        const NekDouble c = std::cos(m_movingFrameData[5]);
        const NekDouble s = std::sin(m_movingFrameData[5]);
        auto projectForce = [c, s](Array<OneD, NekDouble> &force) {
            const NekDouble fx = force[0];
            const NekDouble fy = force[1];
            force[0]          = c * fx - s * fy;
            force[1]          = s * fx + c * fy;
        };

        projectForce(m_pForce);
        projectForce(m_paForce);
        projectForce(m_pqForce);
        projectForce(m_pvisForce);
        projectForce(m_presForce);
    }

    if (m_pressureForceStream.is_open() &&
        m_session->GetComm()->TreatAsRankZero())
    {
        m_pressureForceStream << std::scientific << std::setprecision(10)
                              << time;
        for (int i = 0; i < 3; ++i)
        {
            m_pressureForceStream << " " << m_pForce[i];
        }
        for (int i = 0; i < 3; ++i)
        {
            m_pressureForceStream << " " << m_paForce[i];
        }
        for (int i = 0; i < 3; ++i)
        {
            m_pressureForceStream << " " << m_pqForce[i];
        }
        for (int i = 0; i < 3; ++i)
        {
            m_pressureForceStream << " " << m_pvisForce[i];
        }
        for (int i = 0; i < 3; ++i)
        {
            m_pressureForceStream << " " << m_presForce[i];
        }
        m_pressureForceStream << std::endl;
    }
}

void VCSFSI::OutputPressureComponents([[maybe_unused]] NekDouble time)
{
    if (!m_pressureDecompWriteFld)
    {
        return;
    }

    std::vector<std::string> vars = {"p", "p_a", "p_q", "p_vis"};
    std::vector<Array<OneD, NekDouble>> fields = {
        m_pressure->GetCoeffs(), m_paCoeff, m_pqCoeff, m_pvisCoeff};
    std::string outname = m_sessionName + "_pdec_" +
                          std::to_string(m_pressureDecompOutputIndex) + ".fld";
    WriteFld(outname, m_pressure, fields, vars);
}

void VCSFSI::InitialisePressureComponentForceOutput()
{
    m_pressureForceOutputInitialised = true;

    unsigned int numBoundaryRegions = m_pressure->GetBndConditions().size();
    m_pressureForceBoundaryIsInList.assign(numBoundaryRegions, false);
    Array<OneD, const SpatialDomains::BoundaryConditionShPtr> bndConds =
        m_pressure->GetBndConditions();

    int cnt = 0;
    for (int i = 0; i < bndConds.size(); ++i)
    {
        std::string bcType = bndConds[i]->GetUserDefined();
        if (boost::iequals(bcType, "MRFWallPressDecomp") ||
            boost::iequals(bcType, "MRFWall"))
        {
            m_pressureForceBoundaryIsInList[i] = true;
            ++cnt;
        }
    }

    int globalCnt = cnt;
    m_session->GetComm()->AllReduce(globalCnt, LibUtilities::ReduceSum);
    m_pressureForceHasGlobalBoundary = globalCnt > 0;

    if (!m_pressureForceHasGlobalBoundary)
    {
        return;
    }

    std::string outputBase = m_sessionName;

    std::string suffix = ".mrf";
    if (outputBase.size() >= suffix.size() &&
        outputBase.substr(outputBase.size() - suffix.size()) == suffix)
    {
        outputBase =
            outputBase.substr(0, outputBase.size() - suffix.size());
    }
    outputBase += "_pdec.fce";

    if (m_session->GetComm()->TreatAsRankZero())
    {
        m_pressureForceStream.open(outputBase.c_str());
        m_pressureForceStream
            << "Variables = t, Fpx, Fpy, Fpz, Fax, Fay, Faz, Fqx, Fqy, "
               "Fqz, Fvisx, Fvisy, Fvisz, Fresx, Fresy, Fresz"
            << std::endl;
    }
}

void VCSFSI::IntegratePressureForce(
    const Array<OneD, NekDouble> &pressurePhys,
    Array<OneD, NekDouble> &force) const
{
    Vmath::Zero(force.size(), force, 1);

    const int expdim = m_fields[0]->GetGraph()->GetMeshDimension();
    Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
        m_pressure->GetBndCondExpansions();

    for (int n = 0; n < bndExp.size(); ++n)
    {
        if (n >= static_cast<int>(m_pressureForceBoundaryIsInList.size()) ||
            !m_pressureForceBoundaryIsInList[n])
        {
            continue;
        }

        const int nbc = bndExp[n]->GetTotPoints();
        if (nbc == 0)
        {
            continue;
        }

        Array<OneD, Array<OneD, NekDouble>> normals(3);
        for (int i = 0; i < 3; ++i)
        {
            normals[i] = Array<OneD, NekDouble>(nbc, 0.0);
        }
        bndExp[n]->GetNormals(normals);

        MultiRegions::ExpListSharedPtr bndElmtExp;
        m_pressure->GetBndElmtExpansion(n, bndElmtExp, false);

        Array<OneD, NekDouble> pElm(bndElmtExp->GetTotPoints(), 0.0);
        Array<OneD, NekDouble> pBnd(nbc, 0.0);
        Array<OneD, NekDouble> fp(nbc, 0.0);
        m_pressure->ExtractPhysToBndElmt(n, pressurePhys, pElm);
        m_pressure->ExtractElmtToBndPhys(n, pElm, pBnd);

        int offset = 0;
        for (int e = 0; e < bndExp[n]->GetExpSize(); ++e)
        {
            const int npts = bndExp[n]->GetExp(e)->GetTotPoints();
            for (int i = 0; i < expdim; ++i)
            {
                Array<OneD, NekDouble> pBndSeg = pBnd + offset;
                Array<OneD, NekDouble> normalSeg = normals[i] + offset;
                Array<OneD, NekDouble> fpSeg = fp + offset;
                Vmath::Vmul(npts, pBndSeg, 1, normalSeg, 1, fpSeg, 1);
                force[i] += bndExp[n]->GetExp(e)->Integral(fpSeg);
            }
            offset += npts;
        }

        ASSERTL0(offset == nbc, "Boundary quadrature-point count mismatch.");
    }

    LibUtilities::CommSharedPtr vComm   = m_fields[0]->GetComm();
    LibUtilities::CommSharedPtr rowComm = vComm->GetRowComm();
    LibUtilities::CommSharedPtr colComm =
        m_session->DefinesSolverInfo("HomoStrip")
            ? vComm->GetColumnComm()->GetColumnComm()
            : vComm->GetColumnComm();

    rowComm->AllReduce(force, LibUtilities::ReduceSum);
    colComm->AllReduce(force, LibUtilities::ReduceSum);
}

void VCSFSI::InitialiseFilter(Array<OneD, NekDouble> aeroforce)
{
    if (!m_rigidSolver.HasFreeDoFs())
    {
        return;
    }
    std::map<std::string, std::string> vParams;
    m_rigidSolver.GetFilterInfo(m_session, vParams);
    vParams["OutputFile"] = ".dummyMRFForceFile";
    m_aeroforceFilter =
        MemoryManager<SolverUtils::FilterAeroForces>::AllocateSharedPtr(
            m_session, shared_from_this(), vParams);
    m_aeroforceFilter->Initialise(m_fields, 0.0);
    m_aeroforceFilter->GetForces(m_fields, NullNekDouble1DArray, 0.);
    GetAeroForce(aeroforce);
}

void VCSFSI::v_SolveSolid(NekDouble time)
{
    // call rigid solver
    Array<OneD, NekDouble> aeroforce(12, 0.);
    if (m_rigidSolver.HasFreeDoFs())
    {
        m_aeroforceFilter->GetForces(m_fields, NullNekDouble1DArray, time);
        GetAeroForce(aeroforce);
    }

    // 0-5 pressure force at n+1; 6-11 viscous force at n
    m_rigidSolver.UpdateFrameVelocity(aeroforce, time, m_movingFrameData);
    // update velocity boundary condition
    UpdateVelocityBCs(time);

    CorrectPressureAfterSolid();
}

} // namespace Nektar
