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
#include <LibUtilities/BasicUtils/Filesystem.hpp>
#include <LibUtilities/BasicUtils/Timer.h>
#include <LibUtilities/Foundations/ManagerAccess.h>
#include <MultiRegions/ContField.h>
#include <SolverUtils/Core/Misc.h>

#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <iomanip>
#include <sstream>
#include <vector>

namespace Nektar
{
namespace
{
// Evaluate dF/ds on internal spanwise Gauss sections after integrating the
// quadrature-weighted surface traction only in the chordwise direction.
void AddSpanwiseSectionForce(
    const StdRegions::StdExpansionSharedPtr &face,
    const Array<OneD, const NekDouble> &spanCoords,
    const Array<OneD, Array<OneD, NekDouble>> &weighted,
    const std::vector<NekDouble> &sectionCoords, const int nSections,
    const NekDouble spanTol, Array<OneD, NekDouble> &force)
{
    const int nq0 = face->GetNumPoints(0), nq1 = face->GetNumPoints(1);
    NekDouble dv0 = 0.0, dv1 = 0.0;
    for (int j = 0; j < nq1; ++j)
    {
        for (int i = 0; i < nq0; ++i)
        {
            const int q = i + nq0 * j;
            if (i + 1 < nq0) dv0 = std::max(dv0, std::abs(spanCoords[q]-spanCoords[q+1]));
            if (j + 1 < nq1) dv1 = std::max(dv1, std::abs(spanCoords[q]-spanCoords[q+nq0]));
        }
    }
    const int spanDim = dv0 > dv1 ? 0 : 1;
    const NekDouble spanVariation = std::max(dv0, dv1);
    const NekDouble transverseVariation = std::min(dv0, dv1);
    if (spanVariation <= spanTol)
    {
        // Root/tip end caps have no finite measure in the span direction.
        // Their resultant belongs to the global force, not to dF/ds.
        return;
    }
    ASSERTL0(transverseVariation <= spanTol,
             "Spanwise force output requires faces aligned with SpanForceDir.");
    const int ns = spanDim == 0 ? nq0 : nq1;
    const int nc = spanDim == 0 ? nq1 : nq0;
    const auto &z = face->GetBasis(spanDim)->GetZ();
    const auto &w = face->GetBasis(spanDim)->GetW();
    auto Interpolate = [&z, ns](const Array<OneD, NekDouble> &values,
                                 NekDouble xi)
    {
        NekDouble value = 0.0;
        for (int k = 0; k < ns; ++k)
        {
            NekDouble l = 1.0;
            for (int m = 0; m < ns; ++m)
            {
                if (m != k)
                {
                    l *= (xi - z[m]) / (z[k] - z[m]);
                }
            }
            value += l * values[k];
        }
        return value;
    };
    auto InterpolateDerivative =
        [&z, ns](const Array<OneD, NekDouble> &values, NekDouble xi)
    {
        NekDouble derivative = 0.0;
        for (int k = 0; k < ns; ++k)
        {
            NekDouble basisDerivative = 0.0;
            for (int m = 0; m < ns; ++m)
            {
                if (m == k)
                {
                    continue;
                }
                NekDouble term = 1.0 / (z[k] - z[m]);
                for (int l = 0; l < ns; ++l)
                {
                    if (l != k && l != m)
                    {
                        term *= (xi - z[l]) / (z[k] - z[l]);
                    }
                }
                basisDerivative += term;
            }
            derivative += basisDerivative * values[k];
        }
        return derivative;
    };

    Array<OneD, NekDouble> spanValues(ns, 0.0);
    for (int is = 0; is < ns; ++is)
    {
        spanValues[is] = spanCoords[spanDim == 0 ? is : nq0 * is];
    }

    LibUtilities::PointsKey sectionKey(
        nSections, LibUtilities::eGaussGaussLegendre);
    auto sectionPoints = LibUtilities::PointsManager()[sectionKey];
    const auto &sectionZ = sectionPoints->GetZ();
    for (int d = 0; d < 3; ++d)
    {
        Array<OneD, NekDouble> g(ns, 0.0);
        for (int is = 0; is < ns; ++is)
        {
            for (int ic = 0; ic < nc; ++ic)
            {
                const int q = spanDim == 0 ? is + nq0 * ic : ic + nq0 * is;
                g[is] += weighted[d][q] / w[is];
            }
        }
        for (int section = 0; section < nSections; ++section)
        {
            const NekDouble xi = sectionZ[section];
            const NekDouble span = Interpolate(spanValues, xi);
            const NekDouble dSpanDxi =
                std::abs(InterpolateDerivative(spanValues, xi));
            ASSERTL0(dSpanDxi > spanTol,
                     "Invalid spanwise coordinate metric on section.");
            auto point = std::lower_bound(
                sectionCoords.begin(), sectionCoords.end(), span);
            if (point == sectionCoords.end() ||
                std::abs(*point - span) > spanTol)
            {
                ASSERTL0(point != sectionCoords.begin() &&
                             std::abs(*(point - 1) - span) <= spanTol,
                         "Unable to locate spanwise force section.");
                --point;
            }
            const int pointId =
                static_cast<int>(point - sectionCoords.begin());
            force[3 * pointId + d] +=
                Interpolate(g, xi) / dSpanDxi;
        }
    }
}
}
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
    m_session->MatchSolverInfo("SpanForceOutput", "True",
                               m_spanForceOutput, false);
    m_session->LoadParameter("SpanForceDir", m_spanForceDir, 2);
    m_session->LoadParameter("SpanForceSubstrips", m_spanForceSubstrips, 6);
    ASSERTL0(m_spanForceSubstrips > 0,
             "SpanForceSubstrips must be a positive integer.");

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
        if (m_spanForceOutput)
        {
            ASSERTL0(m_spanForceDir >= 0 && m_spanForceDir < 3,
                     "SpanForceDir must be 0 (x), 1 (y), or 2 (z).");
            m_spanForceOutputDir = m_sessionName + "_spanforce";
        }
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

    if (m_spanForceOutput)
    {
        if (!m_spanForceStripsInitialised)
        {
            InitialiseSpanwiseForceStrips();
        }
        // Integrate once at sub-strip resolution. The coarse strip forces are
        // formed from these same integrals, guaranteeing conservation.
        IntegratePressureForceSpanwisePoints(m_paPhys, m_spanFaPointForce);
        IntegratePressureForceSpanwisePoints(m_pqPhys, m_spanFqPointForce);
        IntegratePressureForceSpanwisePoints(m_pvisPhys,
                                             m_spanFvisPrePointForce);
        IntegratePressureForceSpanwisePoints(m_pressure->GetPhys(),
                                             m_spanFtotalPointForce);
        IntegrateFrictionForceSpanwisePoints(m_spanFfrcPointForce);
        Vmath::Vadd(m_spanFtotalPointForce.size(), m_spanFtotalPointForce, 1,
                    m_spanFfrcPointForce, 1, m_spanFtotalPointForce, 1);

        // Validate the sub-strip integrals against the original spectral-face
        // quadrature. Coarse strips are retained only for this internal check.
        IntegratePressureForceSpanwise(m_paPhys, m_spanFaForce);
        IntegratePressureForceSpanwise(m_pqPhys, m_spanFqForce);
        IntegratePressureForceSpanwise(m_pvisPhys, m_spanFvisPreForce);
        IntegratePressureForceSpanwise(m_pressure->GetPhys(),
                                       m_spanFtotalForce);
        IntegrateFrictionForceSpanwise(m_spanFfrcForce);
        Vmath::Vadd(m_spanFtotalForce.size(), m_spanFtotalForce, 1,
                    m_spanFfrcForce, 1, m_spanFtotalForce, 1);
        CheckSpanwiseSubstripConservation(
            m_spanFaPointForce, m_spanFaForce, "Fa");
        CheckSpanwiseSubstripConservation(
            m_spanFqPointForce, m_spanFqForce, "Fq");
        CheckSpanwiseSubstripConservation(
            m_spanFvisPrePointForce, m_spanFvisPreForce, "Fvis-pre");
        CheckSpanwiseSubstripConservation(
            m_spanFfrcPointForce, m_spanFfrcForce, "Ffrc");
        CheckSpanwiseSubstripConservation(
            m_spanFtotalPointForce, m_spanFtotalForce, "Ftotal");
        if (m_movingFrameData.size() >= 6 && m_movingFrameData[5] != 0.0)
        {
            const NekDouble c = std::cos(m_movingFrameData[5]);
            const NekDouble s = std::sin(m_movingFrameData[5]);
            Array<OneD, NekDouble> *pointForces[] = {
                &m_spanFaPointForce, &m_spanFqPointForce,
                &m_spanFvisPrePointForce, &m_spanFfrcPointForce,
                &m_spanFtotalPointForce};
            for (auto pointForce : pointForces)
            {
                for (int point = 0;
                     point < static_cast<int>(m_spanForcePoints.size());
                     ++point)
                {
                    const NekDouble fx = (*pointForce)[3 * point];
                    const NekDouble fy = (*pointForce)[3 * point + 1];
                    (*pointForce)[3 * point]     = c * fx - s * fy;
                    (*pointForce)[3 * point + 1] = s * fx + c * fy;
                }
            }
        }
        WriteSpanwiseForcePoints(time);
    }

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

void VCSFSI::InitialiseSpanwiseForceStrips()
{
    Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
        m_pressure->GetBndCondExpansions();
    std::vector<NekDouble> localEdges;
    std::vector<NekDouble> localPoints;
    std::vector<NekDouble> localPointWeights;

    for (int n = 0; n < bndExp.size(); ++n)
    {
        if (n >= static_cast<int>(m_pressureForceBoundaryIsInList.size()) ||
            !m_pressureForceBoundaryIsInList[n])
        {
            continue;
        }

        const int nbc = bndExp[n]->GetTotPoints();
        Array<OneD, NekDouble> x(nbc), y(nbc), z(nbc);
        bndExp[n]->GetCoords(x, y, z);
        Array<OneD, NekDouble> coords[] = {x, y, z};
        int offset = 0;
        for (int e = 0; e < bndExp[n]->GetExpSize(); ++e)
        {
            const int npts = bndExp[n]->GetExp(e)->GetTotPoints();
            NekDouble lower = coords[m_spanForceDir][offset];
            NekDouble upper = lower;
            for (int q = 1; q < npts; ++q)
            {
                const NekDouble span = coords[m_spanForceDir][offset + q];
                lower = std::min(lower, span);
                upper = std::max(upper, span);
            }
            const NekDouble faceTol =
                1.0e-10 * std::max(
                    1.0, std::max(std::abs(lower), std::abs(upper)));
            if (upper - lower > faceTol)
            {
                localEdges.push_back(lower);
                localEdges.push_back(upper);
                LibUtilities::PointsKey sectionKey(
                    m_spanForceSubstrips,
                    LibUtilities::eGaussGaussLegendre);
                auto sectionPoints =
                    LibUtilities::PointsManager()[sectionKey];
                const auto &sectionZ = sectionPoints->GetZ();
                const auto &sectionW = sectionPoints->GetW();
                for (int section = 0;
                     section < m_spanForceSubstrips; ++section)
                {
                    const NekDouble span =
                        0.5 * ((1.0 - sectionZ[section]) * lower +
                               (1.0 + sectionZ[section]) * upper);
                    localPoints.push_back(span);
                    localPointWeights.push_back(
                        0.5 * (upper - lower) * sectionW[section]);
                }
            }
            offset += npts;
        }
    }

    auto comm = m_fields[0]->GetComm();
    Array<OneD, int> localSize(1, localEdges.size());
    Array<OneD, int> sizes;
    comm->AllGather(localSize, sizes);
    Array<OneD, int> offsets(sizes.size(), 0);
    int totalSize = 0;
    for (int i = 0; i < sizes.size(); ++i)
    {
        offsets[i] = totalSize;
        totalSize += sizes[i];
    }
    std::vector<NekDouble> globalEdges(totalSize);
    comm->AllGatherv(localEdges, globalEdges, sizes, offsets);

    ASSERTL0(!globalEdges.empty(),
             "No finite-width spanwise boundary faces were found.");
    std::sort(globalEdges.begin(), globalEdges.end());
    const NekDouble scale =
        std::max(1.0, globalEdges.back() - globalEdges.front());
    const NekDouble tol = 1.0e-10 * scale;
    for (const auto edge : globalEdges)
    {
        if (m_spanForceEdges.empty() ||
            std::abs(edge - m_spanForceEdges.back()) > tol)
        {
            m_spanForceEdges.push_back(edge);
        }
    }
    ASSERTL0(m_spanForceEdges.size() >= 2,
             "At least two spanwise mesh lines are required.");

    Array<OneD, int> localPointSize(1, localPoints.size());
    Array<OneD, int> pointSizes;
    comm->AllGather(localPointSize, pointSizes);
    Array<OneD, int> pointOffsets(pointSizes.size(), 0);
    int totalPointSize = 0;
    for (int i = 0; i < pointSizes.size(); ++i)
    {
        pointOffsets[i] = totalPointSize;
        totalPointSize += pointSizes[i];
    }
    std::vector<NekDouble> globalPoints(totalPointSize);
    std::vector<NekDouble> globalPointWeights(totalPointSize);
    comm->AllGatherv(localPoints, globalPoints, pointSizes, pointOffsets);
    comm->AllGatherv(localPointWeights, globalPointWeights,
                     pointSizes, pointOffsets);
    std::vector<std::pair<NekDouble, NekDouble>> pointWeights;
    pointWeights.reserve(totalPointSize);
    for (int i = 0; i < totalPointSize; ++i)
    {
        pointWeights.emplace_back(globalPoints[i], globalPointWeights[i]);
    }
    std::sort(pointWeights.begin(), pointWeights.end());
    for (const auto &pointWeight : pointWeights)
    {
        const NekDouble point = pointWeight.first;
        if (m_spanForcePoints.empty() ||
            std::abs(point - m_spanForcePoints.back()) > tol)
        {
            m_spanForcePoints.push_back(point);
            m_spanForcePointWeights.push_back(pointWeight.second);
        }
    }
    ASSERTL0(m_spanForcePoints.size() >= 2,
             "No spanwise sub-strip boundaries were found.");
    m_spanForceSubstripToStrip.resize(m_spanForcePoints.size());
    for (int sub = 0;
         sub < static_cast<int>(m_spanForcePoints.size()); ++sub)
    {
        const NekDouble centre = m_spanForcePoints[sub];
        auto edge = std::upper_bound(m_spanForceEdges.begin(),
                                     m_spanForceEdges.end(), centre);
        int strip = static_cast<int>(edge - m_spanForceEdges.begin()) - 1;
        strip = std::max(
            0, std::min(strip, static_cast<int>(m_spanForceEdges.size()) - 2));
        ASSERTL0(centre >= m_spanForceEdges[strip] - tol &&
                     centre <= m_spanForceEdges[strip + 1] + tol,
                 "A force section lies outside its spanwise mesh strip.");
        m_spanForceSubstripToStrip[sub] = strip;
    }

    const int n = 3 * (m_spanForceEdges.size() - 1);
    m_spanFaForce      = Array<OneD, NekDouble>(n, 0.0);
    m_spanFqForce      = Array<OneD, NekDouble>(n, 0.0);
    m_spanFvisPreForce = Array<OneD, NekDouble>(n, 0.0);
    m_spanFfrcForce    = Array<OneD, NekDouble>(n, 0.0);
    m_spanFtotalForce  = Array<OneD, NekDouble>(n, 0.0);
    const int np = 3 * m_spanForcePoints.size();
    m_spanFaPointForce      = Array<OneD, NekDouble>(np, 0.0);
    m_spanFqPointForce      = Array<OneD, NekDouble>(np, 0.0);
    m_spanFvisPrePointForce = Array<OneD, NekDouble>(np, 0.0);
    m_spanFfrcPointForce    = Array<OneD, NekDouble>(np, 0.0);
    m_spanFtotalPointForce  = Array<OneD, NekDouble>(np, 0.0);
    if (m_session->GetComm()->TreatAsRankZero())
    {
        fs::create_directories(m_spanForceOutputDir);
    }
    m_spanForceStripsInitialised = true;
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

void VCSFSI::IntegratePressureForceSpanwise(
    const Array<OneD, NekDouble> &pressurePhys,
    Array<OneD, NekDouble> &force) const
{
    Vmath::Zero(force.size(), force, 1);
    const int expdim = m_fields[0]->GetGraph()->GetMeshDimension();
    const NekDouble spanTol = 1.0e-10 * std::max(
        1.0, m_spanForceEdges.back() - m_spanForceEdges.front());
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
        Array<OneD, NekDouble> pBnd(nbc, 0.0), x(nbc), y(nbc), z(nbc);
        m_pressure->ExtractPhysToBndElmt(n, pressurePhys, pElm);
        m_pressure->ExtractElmtToBndPhys(n, pElm, pBnd);
        bndExp[n]->GetCoords(x, y, z);
        Array<OneD, NekDouble> coords[] = {x, y, z};

        int offset = 0;
        for (int e = 0; e < bndExp[n]->GetExpSize(); ++e)
        {
            auto bc = bndExp[n]->GetExp(e);
            const int npts = bc->GetTotPoints();
            NekDouble spanLower = coords[m_spanForceDir][offset];
            NekDouble spanUpper = spanLower;
            for (int q = 0; q < npts; ++q)
            {
                const NekDouble span = coords[m_spanForceDir][offset + q];
                spanLower = std::min(spanLower, span);
                spanUpper = std::max(spanUpper, span);
            }
            if (spanUpper - spanLower <= spanTol)
            {
                offset += npts;
                continue;
            }
            const NekDouble spanCentre = 0.5 * (spanLower + spanUpper);
            auto upper = std::upper_bound(m_spanForceEdges.begin(),
                                          m_spanForceEdges.end(),
                                          spanCentre + spanTol);
            const int bin = std::max(
                0, std::min(static_cast<int>(upper - m_spanForceEdges.begin()) - 1,
                            static_cast<int>(m_spanForceEdges.size()) - 2));
            ASSERTL0(spanLower >= m_spanForceEdges[bin] - spanTol &&
                         spanUpper <= m_spanForceEdges[bin + 1] + spanTol,
                     "A spanwise boundary face lies outside the configured "
                     "spanwise mesh strip.");

            for (int d = 0; d < expdim; ++d)
            {
                Array<OneD, NekDouble> value(npts, 0.0);
                for (int q = 0; q < npts; ++q)
                {
                    value[q] = pBnd[offset + q] * normals[d][offset + q];
                }
                force[3 * bin + d] += bc->Integral(value);
            }
            offset += npts;
        }
        ASSERTL0(offset == nbc, "Boundary quadrature-point count mismatch.");
    }

    auto comm = m_fields[0]->GetComm();
    comm->GetRowComm()->AllReduce(force, LibUtilities::ReduceSum);
    auto colComm = m_session->DefinesSolverInfo("HomoStrip")
                       ? comm->GetColumnComm()->GetColumnComm()
                       : comm->GetColumnComm();
    colComm->AllReduce(force, LibUtilities::ReduceSum);
}

void VCSFSI::CheckSpanwiseSubstripConservation(
    const Array<OneD, NekDouble> &substripForce,
    const Array<OneD, NekDouble> &stripForce,
    const std::string &forceName) const
{
    ASSERTL0(substripForce.size() ==
                 3 * static_cast<int>(m_spanForceSubstripToStrip.size()),
             "Unexpected spanwise section-force array size.");
    ASSERTL0(m_spanForcePointWeights.size() ==
                 m_spanForceSubstripToStrip.size(),
             "Unexpected spanwise section-weight array size.");
    ASSERTL0(stripForce.size() ==
                 3 * (static_cast<int>(m_spanForceEdges.size()) - 1),
             "Unexpected spanwise mesh-strip force array size.");
    Array<OneD, NekDouble> accumulated(stripForce.size(), 0.0);
    for (int sub = 0;
         sub < static_cast<int>(m_spanForceSubstripToStrip.size()); ++sub)
    {
        const int strip = m_spanForceSubstripToStrip[sub];
        for (int d = 0; d < 3; ++d)
        {
            accumulated[3 * strip + d] +=
                substripForce[3 * sub + d] *
                m_spanForcePointWeights[sub];
        }
    }
    NekDouble maxForce = 0.0, maxError = 0.0;
    for (int i = 0; i < stripForce.size(); ++i)
    {
        maxForce = std::max(maxForce, std::abs(stripForce[i]));
        maxError = std::max(
            maxError, std::abs(accumulated[i] - stripForce[i]));
    }
    const NekDouble tolerance = 1.0e-10 * std::max(1.0, maxForce);
    ASSERTL0(maxError <= tolerance,
             "Spanwise section-force conservation check failed for " +
                 forceName + ": error = " + std::to_string(maxError) +
                 ", tolerance = " + std::to_string(tolerance));
}

void VCSFSI::IntegratePressureForceSpanwisePoints(
    const Array<OneD, NekDouble> &pressurePhys,
    Array<OneD, NekDouble> &force) const
{
    Vmath::Zero(force.size(), force, 1);
    const int expdim = m_fields[0]->GetGraph()->GetMeshDimension();
    const NekDouble spanTol = 1.0e-10 * std::max(
        1.0, m_spanForceEdges.back() - m_spanForceEdges.front());
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
        for (int d = 0; d < 3; ++d)
        {
            normals[d] = Array<OneD, NekDouble>(nbc, 0.0);
        }
        bndExp[n]->GetNormals(normals);

        MultiRegions::ExpListSharedPtr bndElmtExp;
        m_pressure->GetBndElmtExpansion(n, bndElmtExp, false);
        Array<OneD, NekDouble> pElm(bndElmtExp->GetTotPoints(), 0.0);
        Array<OneD, NekDouble> pBnd(nbc, 0.0), x(nbc), y(nbc), z(nbc);
        m_pressure->ExtractPhysToBndElmt(n, pressurePhys, pElm);
        m_pressure->ExtractElmtToBndPhys(n, pElm, pBnd);
        bndExp[n]->GetCoords(x, y, z);
        Array<OneD, NekDouble> coords[] = {x, y, z};

        int offset = 0;
        for (int e = 0; e < bndExp[n]->GetExpSize(); ++e)
        {
            auto bc = bndExp[n]->GetExp(e);
            const int nq0 = bc->GetNumPoints(0);
            const int nq1 = bc->GetNumPoints(1);
            const int npts = bc->GetTotPoints();
            ASSERTL0(npts == nq0 * nq1,
                     "Spanwise point output requires quadrilateral boundary "
                     "faces.");

            NekDouble variation0 = 0.0, variation1 = 0.0;
            NekDouble spanLower = coords[m_spanForceDir][offset];
            NekDouble spanUpper = spanLower;
            for (int j = 0; j < nq1; ++j)
            {
                for (int i = 0; i < nq0; ++i)
                {
                    const int q = i + nq0 * j;
                    const NekDouble span = coords[m_spanForceDir][offset + q];
                    spanLower = std::min(spanLower, span);
                    spanUpper = std::max(spanUpper, span);
                    if (i + 1 < nq0)
                    {
                        variation0 = std::max(
                            variation0, std::abs(span - coords[m_spanForceDir]
                                                       [offset + q + 1]));
                    }
                    if (j + 1 < nq1)
                    {
                        variation1 = std::max(
                            variation1, std::abs(span - coords[m_spanForceDir]
                                                       [offset + q + nq0]));
                    }
                }
            }
            if (spanUpper - spanLower <= spanTol)
            {
                // End-cap force is retained in the global surface force, but
                // cannot be represented as a finite spanwise line load.
                offset += npts;
                continue;
            }
            const int spanDim = variation0 > variation1 ? 0 : 1;
            const NekDouble transverseVariation =
                spanDim == 0 ? variation1 : variation0;
            ASSERTL0(transverseVariation < spanTol,
                     "Spanwise point output requires boundary faces aligned "
                     "with SpanForceDir.");
            Array<OneD, Array<OneD, NekDouble>> substripWeighted(3);
            for (int d = 0; d < 3; ++d)
            {
                substripWeighted[d] = Array<OneD, NekDouble>(npts, 0.0);
                if (d < expdim)
                {
                    Array<OneD, NekDouble> value(npts, 0.0);
                    for (int q = 0; q < npts; ++q)
                    {
                        value[q] = pBnd[offset + q] * normals[d][offset + q];
                    }
                    bc->MultiplyByQuadratureMetric(value, substripWeighted[d]);
                }
            }
            Array<OneD, const NekDouble> span =
                coords[m_spanForceDir] + offset;
            AddSpanwiseSectionForce(bc, span, substripWeighted,
                                    m_spanForcePoints, m_spanForceSubstrips,
                                    spanTol, force);
            offset += npts;
        }
        ASSERTL0(offset == nbc, "Boundary quadrature-point count mismatch.");
    }

    auto comm = m_fields[0]->GetComm();
    comm->GetRowComm()->AllReduce(force, LibUtilities::ReduceSum);
    auto colComm = m_session->DefinesSolverInfo("HomoStrip")
                       ? comm->GetColumnComm()->GetColumnComm()
                       : comm->GetColumnComm();
    colComm->AllReduce(force, LibUtilities::ReduceSum);
}

void VCSFSI::WriteSpanwiseForcePoints(NekDouble time) const
{
    if (!m_session->GetComm()->TreatAsRankZero())
    {
        return;
    }
    std::ostringstream name;
    name << m_spanForceOutputDir << "/spanforce_" << std::setw(6)
         << std::setfill('0') << m_spanForceOutputIndex++ << ".dat";
    std::ofstream out(name.str());
    ASSERTL0(out.good(), "Unable to open spanwise sub-strip output file.");
    out << std::scientific << std::setprecision(12);
    out << "# time " << time << "\n";
    out << "# s Fa_x Fa_y Fa_z Fq_x Fq_y Fq_z "
           "Fvispre_x Fvispre_y Fvispre_z Ffrc_x Ffrc_y Ffrc_z "
           "Ftotal_x Ftotal_y Ftotal_z\n";
    const Array<OneD, NekDouble> *forces[] = {
        &m_spanFaPointForce, &m_spanFqPointForce, &m_spanFvisPrePointForce,
        &m_spanFfrcPointForce, &m_spanFtotalPointForce};
    for (int point = 0; point < static_cast<int>(m_spanForcePoints.size());
         ++point)
    {
        out << m_spanForcePoints[point];
        for (const auto force : forces)
        {
            for (int d = 0; d < 3; ++d)
            {
                out << " " << (*force)[3 * point + d];
            }
        }
        out << "\n";
    }
}

void VCSFSI::IntegrateFrictionForceSpanwise(
    Array<OneD, NekDouble> &force) const
{
    Vmath::Zero(force.size(), force, 1);
    const int expdim = m_fields[0]->GetGraph()->GetMeshDimension();
    ASSERTL0(expdim == 3,
             "Spanwise force output currently requires a full 3D mesh.");
    const NekDouble spanTol = 1.0e-10 * std::max(
        1.0, m_spanForceEdges.back() - m_spanForceEdges.front());
    const NekDouble rho = m_session->DefinesParameter("rho")
                              ? m_session->GetParameter("rho")
                              : 1.0;
    const NekDouble mu = m_session->DefinesParameter("Kinvis")
                             ? rho * m_session->GetParameter("Kinvis")
                             : m_session->GetParameter("mu");

    Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
        m_fields[0]->GetBndCondExpansions();
    Array<OneD, int> bcToElmt, bcToTrace;
    m_fields[0]->GetBoundaryToElmtMap(bcToElmt, bcToTrace);
    int cnt = 0;

    for (int n = 0; n < bndExp.size(); ++n)
    {
        if (n >= static_cast<int>(m_pressureForceBoundaryIsInList.size()) ||
            !m_pressureForceBoundaryIsInList[n])
        {
            cnt += bndExp[n]->GetExpSize();
            continue;
        }
        for (int i = 0; i < bndExp[n]->GetExpSize(); ++i, ++cnt)
        {
            auto elmt = m_fields[0]->GetExp(bcToElmt[cnt]);
            const int trace = bcToTrace[cnt];
            const int nq = elmt->GetTotPoints();
            const int offset = m_fields[0]->GetPhys_Offset(bcToElmt[cnt]);
            Array<OneD, Array<OneD, NekDouble>> grad(9);
            for (int j = 0; j < 3; ++j)
            {
                Array<OneD, const NekDouble> vel =
                    m_fields[m_velocity[j]]->GetPhys() + offset;
                for (int k = 0; k < 3; ++k)
                {
                    grad[3 * j + k] = Array<OneD, NekDouble>(nq, 0.0);
                    elmt->PhysDeriv(k, vel, grad[3 * j + k]);
                }
            }

            auto bc = bndExp[n]->GetExp(i);
            const int nbc = bc->GetTotPoints();
            auto normals = elmt->GetTraceNormal(trace);
            Array<OneD, Array<OneD, NekDouble>> gradb(9);
            Array<OneD, Array<OneD, NekDouble>> coords(3), coordsb(3);
            for (int d = 0; d < 3; ++d)
            {
                coords[d] = Array<OneD, NekDouble>(nq, 0.0);
                coordsb[d] = Array<OneD, NekDouble>(nbc, 0.0);
            }
            elmt->GetCoords(coords[0], coords[1], coords[2]);
            for (int j = 0; j < 9; ++j)
            {
                gradb[j] = Array<OneD, NekDouble>(nbc, 0.0);
                elmt->GetTracePhysVals(trace, bc, grad[j], gradb[j]);
            }
            for (int d = 0; d < 3; ++d)
            {
                elmt->GetTracePhysVals(trace, bc, coords[d], coordsb[d]);
            }

            NekDouble spanLower = coordsb[m_spanForceDir][0];
            NekDouble spanUpper = spanLower;
            for (int q = 0; q < nbc; ++q)
            {
                const NekDouble span = coordsb[m_spanForceDir][q];
                spanLower = std::min(spanLower, span);
                spanUpper = std::max(spanUpper, span);
            }
            if (spanUpper - spanLower <= spanTol)
            {
                continue;
            }
            const NekDouble spanCentre = 0.5 * (spanLower + spanUpper);
            auto upper = std::upper_bound(m_spanForceEdges.begin(),
                                          m_spanForceEdges.end(),
                                          spanCentre + spanTol);
            const int bin = std::max(
                0, std::min(static_cast<int>(upper - m_spanForceEdges.begin()) - 1,
                            static_cast<int>(m_spanForceEdges.size()) - 2));
            ASSERTL0(spanLower >= m_spanForceEdges[bin] - spanTol &&
                         spanUpper <= m_spanForceEdges[bin + 1] + spanTol,
                     "A spanwise boundary face lies outside the configured "
                     "spanwise mesh strip.");

            for (int d = 0; d < 3; ++d)
            {
                Array<OneD, NekDouble> traction(nbc, 0.0);
                for (int k = 0; k < 3; ++k)
                {
                    Vmath::Vvtvp(nbc, gradb[3 * k + d], 1, normals[k], 1,
                                 traction, 1, traction, 1);
                    Vmath::Vvtvp(nbc, gradb[3 * d + k], 1, normals[k], 1,
                                 traction, 1, traction, 1);
                }
                Vmath::Smul(nbc, -mu, traction, 1, traction, 1);
                force[3 * bin + d] += bc->Integral(traction);
            }
        }
    }

    auto comm = m_fields[0]->GetComm();
    comm->GetRowComm()->AllReduce(force, LibUtilities::ReduceSum);
    auto colComm = m_session->DefinesSolverInfo("HomoStrip")
                       ? comm->GetColumnComm()->GetColumnComm()
                       : comm->GetColumnComm();
    colComm->AllReduce(force, LibUtilities::ReduceSum);
}

void VCSFSI::IntegrateFrictionForceSpanwisePoints(
    Array<OneD, NekDouble> &force) const
{
    Vmath::Zero(force.size(), force, 1);
    const int expdim = m_fields[0]->GetGraph()->GetMeshDimension();
    ASSERTL0(expdim == 3,
             "Spanwise force output currently requires a full 3D mesh.");
    const NekDouble spanTol = 1.0e-10 * std::max(
        1.0, m_spanForceEdges.back() - m_spanForceEdges.front());
    const NekDouble rho = m_session->DefinesParameter("rho")
                              ? m_session->GetParameter("rho")
                              : 1.0;
    const NekDouble mu = m_session->DefinesParameter("Kinvis")
                             ? rho * m_session->GetParameter("Kinvis")
                             : m_session->GetParameter("mu");

    Array<OneD, MultiRegions::ExpListSharedPtr> bndExp =
        m_fields[0]->GetBndCondExpansions();
    Array<OneD, int> bcToElmt, bcToTrace;
    m_fields[0]->GetBoundaryToElmtMap(bcToElmt, bcToTrace);
    int cnt = 0;
    for (int n = 0; n < bndExp.size(); ++n)
    {
        if (n >= static_cast<int>(m_pressureForceBoundaryIsInList.size()) ||
            !m_pressureForceBoundaryIsInList[n])
        {
            cnt += bndExp[n]->GetExpSize();
            continue;
        }
        for (int i = 0; i < bndExp[n]->GetExpSize(); ++i, ++cnt)
        {
            auto elmt = m_fields[0]->GetExp(bcToElmt[cnt]);
            const int trace = bcToTrace[cnt];
            const int nq = elmt->GetTotPoints();
            const int offset = m_fields[0]->GetPhys_Offset(bcToElmt[cnt]);
            Array<OneD, Array<OneD, NekDouble>> grad(9);
            for (int j = 0; j < 3; ++j)
            {
                Array<OneD, const NekDouble> vel =
                    m_fields[m_velocity[j]]->GetPhys() + offset;
                for (int k = 0; k < 3; ++k)
                {
                    grad[3 * j + k] = Array<OneD, NekDouble>(nq, 0.0);
                    elmt->PhysDeriv(k, vel, grad[3 * j + k]);
                }
            }

            auto bc = bndExp[n]->GetExp(i);
            const int nq0 = bc->GetNumPoints(0);
            const int nq1 = bc->GetNumPoints(1);
            const int nbc = bc->GetTotPoints();
            ASSERTL0(nbc == nq0 * nq1,
                     "Spanwise point output requires quadrilateral boundary "
                     "faces.");
            auto normals = elmt->GetTraceNormal(trace);
            Array<OneD, Array<OneD, NekDouble>> gradb(9);
            Array<OneD, Array<OneD, NekDouble>> coords(3), coordsb(3);
            for (int d = 0; d < 3; ++d)
            {
                coords[d] = Array<OneD, NekDouble>(nq, 0.0);
                coordsb[d] = Array<OneD, NekDouble>(nbc, 0.0);
            }
            elmt->GetCoords(coords[0], coords[1], coords[2]);
            for (int j = 0; j < 9; ++j)
            {
                gradb[j] = Array<OneD, NekDouble>(nbc, 0.0);
                elmt->GetTracePhysVals(trace, bc, grad[j], gradb[j]);
            }
            for (int d = 0; d < 3; ++d)
            {
                elmt->GetTracePhysVals(trace, bc, coords[d], coordsb[d]);
            }

            NekDouble variation0 = 0.0, variation1 = 0.0;
            NekDouble spanLower = coordsb[m_spanForceDir][0];
            NekDouble spanUpper = spanLower;
            for (int j = 0; j < nq1; ++j)
            {
                for (int k = 0; k < nq0; ++k)
                {
                    const int q = k + nq0 * j;
                    const NekDouble span = coordsb[m_spanForceDir][q];
                    spanLower = std::min(spanLower, span);
                    spanUpper = std::max(spanUpper, span);
                    if (k + 1 < nq0)
                    {
                        variation0 = std::max(
                            variation0, std::abs(span - coordsb[m_spanForceDir]
                                                       [q + 1]));
                    }
                    if (j + 1 < nq1)
                    {
                        variation1 = std::max(
                            variation1, std::abs(span - coordsb[m_spanForceDir]
                                                       [q + nq0]));
                    }
                }
            }
            if (spanUpper - spanLower <= spanTol)
            {
                continue;
            }
            const int spanDim = variation0 > variation1 ? 0 : 1;
            const NekDouble transverseVariation =
                spanDim == 0 ? variation1 : variation0;
            ASSERTL0(transverseVariation < spanTol,
                     "Spanwise point output requires boundary faces aligned "
                     "with SpanForceDir.");
            Array<OneD, Array<OneD, NekDouble>> substripWeighted(3);
            for (int d = 0; d < 3; ++d)
            {
                substripWeighted[d] = Array<OneD, NekDouble>(nbc, 0.0);
                Array<OneD, NekDouble> traction(nbc, 0.0);
                for (int k = 0; k < 3; ++k)
                {
                    Vmath::Vvtvp(nbc, gradb[3 * k + d], 1, normals[k], 1,
                                 traction, 1, traction, 1);
                    Vmath::Vvtvp(nbc, gradb[3 * d + k], 1, normals[k], 1,
                                 traction, 1, traction, 1);
                }
                Vmath::Smul(nbc, -mu, traction, 1, traction, 1);
                bc->MultiplyByQuadratureMetric(traction, substripWeighted[d]);
            }
            Array<OneD, const NekDouble> span =
                coordsb[m_spanForceDir];
            AddSpanwiseSectionForce(bc, span, substripWeighted,
                                    m_spanForcePoints, m_spanForceSubstrips,
                                    spanTol, force);
        }
    }

    auto comm = m_fields[0]->GetComm();
    comm->GetRowComm()->AllReduce(force, LibUtilities::ReduceSum);
    auto colComm = m_session->DefinesSolverInfo("HomoStrip")
                       ? comm->GetColumnComm()->GetColumnComm()
                       : comm->GetColumnComm();
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
