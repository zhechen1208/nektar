///////////////////////////////////////////////////////////////////////////////
//
// File: VCSFSI.h
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
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VCSFSI_H
#define NEKTAR_SOLVERS_VCSFSI_H

#include <IncNavierStokesSolver/EquationSystems/RigidSolver.h>
#include <IncNavierStokesSolver/EquationSystems/VelocityCorrectionScheme.h>
#include <SolverUtils/Filters/FilterAeroForces.h>
#include <fstream>
#include <vector>
namespace Nektar
{
class VCSFSI : public VelocityCorrectionScheme
{
public:
    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create(
        const LibUtilities::SessionReaderSharedPtr &pSession,
        const SpatialDomains::MeshGraphSharedPtr &pGraph)
    {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<VCSFSI>::AllocateSharedPtr(pSession, pGraph);
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;

    /// Constructor.
    VCSFSI(const LibUtilities::SessionReaderSharedPtr &pSession,
           const SpatialDomains::MeshGraphSharedPtr &pGraph);

    ~VCSFSI() override;

    void v_InitObject(bool DeclareField = true) override;

protected:
    static std::string solverTypeLookupId;

    bool m_verbose;
    bool m_enablePressureDecomposition;
    bool m_pressureDecompWriteFld;
    int m_pressureDecompOutputFrequency;
    int m_pressureDecompOutputIndex;

    // Virtual functions
    void v_DoInitialise(bool dumpInitialConditions = true) override;
    void v_SetUpPressureForcing(
        const Array<OneD, const Array<OneD, NekDouble>> &fields,
        Array<OneD, Array<OneD, NekDouble>> &Forcing,
        NekDouble aiiDt) override;

    void v_SolveSolid(NekDouble time) override;
    void InitialiseFilter(Array<OneD, NekDouble> aeroforce);
    virtual void CorrectPressureAfterSolid();
    virtual void InitialisePressureDecomposition();
    virtual void UpdatePressureDecomposition(NekDouble time);
    virtual void ComputePaFull(NekDouble time);
    virtual void ComputePq(NekDouble time);
    virtual void ComputePvis(NekDouble time);
    virtual void EvaluatePressureComponentForces(NekDouble time);
    virtual void OutputPressureComponents(NekDouble time);
    virtual void InitialisePressureComponentForceOutput();
    virtual void IntegratePressureForce(
        const Array<OneD, NekDouble> &pressurePhys,
        Array<OneD, NekDouble> &force) const;
    virtual void IntegrateFrictionForce(
        Array<OneD, NekDouble> &force) const;
    virtual void IntegratePressureForceSpanwise(
        const Array<OneD, NekDouble> &pressurePhys,
        Array<OneD, NekDouble> &force) const;
    virtual void IntegratePressureForceTipCap(
        const Array<OneD, NekDouble> &pressurePhys,
        Array<OneD, NekDouble> &force) const;
    virtual void IntegrateFrictionForceSpanwise(
        Array<OneD, NekDouble> &force) const;
    virtual void IntegrateFrictionForceTipCap(
        Array<OneD, NekDouble> &force) const;
    virtual void InitialiseSpanwiseForceStrips();
    virtual void IntegratePressureForceSpanwisePoints(
        const Array<OneD, NekDouble> &pressurePhys,
        Array<OneD, NekDouble> &force) const;
    virtual void IntegrateFrictionForceSpanwisePoints(
        Array<OneD, NekDouble> &force) const;
    virtual void WriteSpanwiseForcePoints(NekDouble time) const;
    virtual void CheckSpanwiseSubstripConservation(
        const Array<OneD, NekDouble> &substripForce,
        const Array<OneD, NekDouble> &stripForce,
        const std::string &forceName) const;
    virtual void CheckSpanwiseTipCapConservation(
        const Array<OneD, NekDouble> &tipCapForce,
        const Array<OneD, NekDouble> &sideForce,
        const Array<OneD, NekDouble> &globalForce,
        const std::string &forceName) const;
    void ZeroPressureBoundaryConditions();

    RigidSolver m_rigidSolver;
    SolverUtils::FilterAeroForcesSharedPtr m_aeroforceFilter;
    MultiRegions::ExpListSharedPtr m_pressureDecomp;
    Array<OneD, NekDouble> m_paCoeff;
    Array<OneD, NekDouble> m_pqCoeff;
    Array<OneD, NekDouble> m_pvisCoeff;
    Array<OneD, NekDouble> m_paPhys;
    Array<OneD, NekDouble> m_pqPhys;
    Array<OneD, NekDouble> m_pvisPhys;
    Array<OneD, NekDouble> m_pressurePoissonRhs;
    Array<OneD, NekDouble> m_paForce;
    Array<OneD, NekDouble> m_pqForce;
    Array<OneD, NekDouble> m_pvisForce;
    Array<OneD, NekDouble> m_pForce;
    Array<OneD, NekDouble> m_frictionForce;
    Array<OneD, NekDouble> m_totalForce;
    std::vector<bool> m_pressureForceBoundaryIsInList;
    bool m_pressureForceOutputInitialised = false;
    bool m_pressureForceHasGlobalBoundary = false;
    std::ofstream m_pressureForceStream;
    bool m_spanForceOutput = false;
    bool m_spanForceStripsInitialised = false;
    int m_spanForceDir = 2;
    int m_spanForceSubstrips = 6;
    std::vector<NekDouble> m_spanForceEdges;
    std::vector<NekDouble> m_spanForcePoints;
    std::vector<NekDouble> m_spanForcePointWeights;
    std::vector<int> m_spanForceSubstripToStrip;
    std::string m_spanForceOutputDir;
    mutable int m_spanForceOutputIndex = 0;
    Array<OneD, NekDouble> m_spanFaForce;
    Array<OneD, NekDouble> m_spanFqForce;
    Array<OneD, NekDouble> m_spanFvisPreForce;
    Array<OneD, NekDouble> m_spanFfrcForce;
    Array<OneD, NekDouble> m_spanFtotalForce;
    Array<OneD, NekDouble> m_spanFaPointForce;
    Array<OneD, NekDouble> m_spanFqPointForce;
    Array<OneD, NekDouble> m_spanFvisPrePointForce;
    Array<OneD, NekDouble> m_spanFfrcPointForce;
    Array<OneD, NekDouble> m_spanFtotalPointForce;
    Array<OneD, NekDouble> m_tipCapFaForce;
    Array<OneD, NekDouble> m_tipCapFqForce;
    Array<OneD, NekDouble> m_tipCapFvisPreForce;
    Array<OneD, NekDouble> m_tipCapFfrcForce;
    Array<OneD, NekDouble> m_tipCapFtotalForce;
};

typedef std::shared_ptr<VCSFSI> VCSFSISharedPtr;

} // namespace Nektar

#endif // VCSFSI_H
