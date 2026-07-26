///////////////////////////////////////////////////////////////////////////////
//
// File: RigidSolver.h
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
// Description: Allows for a moving frame of reference, through adding c * du/dx
// to the body force, where c is the frame velocity vector
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_RIGIDSOLVER
#define NEKTAR_SOLVERS_RIGIDSOLVER

#include <string>

#include <LibUtilities/BasicUtils/Equation.h>
#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/EquationSystem.h>
#include <SolverUtils/Forcing/Forcing.h>
#include <SolverUtils/SolverUtilsDeclspec.h>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <cmath>

namespace Nektar
{

enum RigidSolveType
{
    eInertialTranslation = 0,
    ePlanarRigidBody     = 2,
    ePrescribedMRF       = 4,
    eFree3D6DoF          = 5
};

/***
 * Solve the body's motion using Newmark-Beta method
 * M ddx + C dx + K x = F
 * In discrete form
 * CoeffMatrix dx = rhs
 ***/
class Newmark_BetaSolver
{
public:
    Newmark_BetaSolver() {};
    ~Newmark_BetaSolver() {};
    void SetNewmarkBeta(NekDouble beta, NekDouble gamma, NekDouble dt,
                        Array<OneD, NekDouble> M, Array<OneD, NekDouble> C,
                        Array<OneD, NekDouble> K, std::set<int> DirDoFs,
                        int solveType);
    void SolvePrescribed(Array<OneD, Array<OneD, NekDouble>> u,
                         std::map<int, NekDouble> motionPrescribed);
    void SolveFreeFixMat(Array<OneD, Array<OneD, NekDouble>> u,
                         Array<OneD, NekDouble> force);
    void SolveFreeVarMat6DoF(
        Array<OneD, Array<OneD, NekDouble>> u,
        const Array<OneD, NekDouble> &force,
        const Array<OneD, NekDouble> &nonlinearTerm,
        const Array<OneD, NekDouble> &nonlinearJacobian,
        const Array<OneD, NekDouble> &linearisationVelocity);
    void SolveFreeVarMatNDof(
        Array<OneD, Array<OneD, NekDouble>> u,
        const Array<OneD, NekDouble> &force,
        const Array<OneD, NekDouble> &nonlinearTerm,
        const Array<OneD, NekDouble> &nonlinearJacobian,
        const Array<OneD, NekDouble> &linearisationVelocity, int nDofs);
    void SolveFreeVarMatNDofConstrained(
        Array<OneD, Array<OneD, NekDouble>> u,
        const Array<OneD, NekDouble> &force,
        const Array<OneD, NekDouble> &nonlinearTerm,
        const Array<OneD, NekDouble> &nonlinearJacobian,
        const Array<OneD, NekDouble> &linearisationVelocity,
        const Array<OneD, NekDouble> &velocityConstraints,
        const Array<OneD, NekDouble> &constraintVelocity, int nDofs,
        int nConstraints);
    int m_rows;
    int m_motionDofs;
    std::vector<int> m_index;
    Array<OneD, NekDouble> m_coeffs;
    Array<OneD, Array<OneD, NekDouble>> m_Matrix;
    Array<OneD, Array<OneD, NekDouble>> m_M;
    Array<OneD, Array<OneD, NekDouble>> m_C;
    Array<OneD, Array<OneD, NekDouble>> m_K;
};

class FrameTransform
{
public:
    FrameTransform();
    ~FrameTransform() {};
    void SetAngle(const Array<OneD, NekDouble> theta);
    void BodyToInerital(const int dim, const Array<OneD, NekDouble> &body,
                        Array<OneD, NekDouble> &inertial);
    void IneritalToBody(const int dim, const Array<OneD, NekDouble> &inertial,
                        Array<OneD, NekDouble> &body);

private:
    Array<OneD, NekDouble> m_matrix;
};

class RigidSolver
{
public:
    void InitObject(const LibUtilities::SessionReaderSharedPtr session,
                    const MultiRegions::ExpListSharedPtr &pField,
                    Array<OneD, NekDouble> pivot);
    void SetMovableDoFs(std::vector<bool> &moveDoFs);
    void SetInitialConditions(
        const LibUtilities::SessionReaderSharedPtr session,
        Array<OneD, NekDouble> MRFData);
    void UpdateFrameVelocity(Array<OneD, NekDouble> &aeroforce,
                             const NekDouble &time,
                             Array<OneD, NekDouble> &MRFData);
    void SetNewmarkBetaSolver(Array<OneD, NekDouble> &AddedMass);
    void GetFilterInfo(const LibUtilities::SessionReaderSharedPtr session,
                       std::map<std::string, std::string> &vParams);
    inline bool HasFreeDoFs()
    {
        return m_hasFreeMotion;
    };
    inline bool IsFullFree3D6DoF() const
    {
        return m_free3D6DoF;
    }
    void SetOldFvis(Array<OneD, NekDouble> force);

protected:
    void LoadParameters(const LibUtilities::SessionReaderSharedPtr session,
                        const TiXmlElement *pSolver);
    void InitBodySolver(const LibUtilities::SessionReaderSharedPtr session,
                        const TiXmlElement *pSolver,
                        Array<OneD, NekDouble> pivot);
    NekDouble EvaluateExpression(
        const LibUtilities::SessionReaderSharedPtr session,
        std::string expression);
    void UpdatePrescribed(const NekDouble &time,
                          std::map<int, NekDouble> &Dirs);
    void UpdateMRFData(Array<OneD, NekDouble> &MRFData);
    void SetInitialConditions(std::map<int, NekDouble> &Dirs);
    void SolveBodyMotion(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                         const Array<OneD, NekDouble> &forcebody,
                         std::map<int, NekDouble> &Dirs);
    void ParserFunctionToMap(
        const bool allowzero,
        const LibUtilities::SessionReaderSharedPtr session,
        const std::string &FuncName, const std::string &var,
        std::map<int, LibUtilities::EquationSharedPtr> &result, const int i);
    void CheckParameters();
    bool HasFull3DPrescribedOrientation() const;
    void LoadThetaConvention(const TiXmlElement *pSolver);
    Array<OneD, NekDouble> QuaternionFromConfiguredTheta(
        const Array<OneD, NekDouble> &theta) const;
    void UpdatePrescribedMRFData(const NekDouble &time,
                                 Array<OneD, NekDouble> &MRFData);
    void UpdateFree3DMRFData(Array<OneD, NekDouble> &MRFData);
    void WritePrescribedMRFOutput(const NekDouble &time,
                                  const Array<OneD, NekDouble> &MRFData);
    void WriteFree3DMRFOutput(const NekDouble &time,
                              const Array<OneD, NekDouble> &MRFData);

private:
    RigidSolveType m_solveType;
    // eInertialTranslation: translation only;
    // ePlanarRigidBody: unified planar body-frame solve.
    void SolveInertialFrame(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                            const Array<OneD, NekDouble> &forcebody,
                            std::map<int, NekDouble> &Dirs);
    void SolveFree3D6DoF(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                         const Array<OneD, NekDouble> &forcebody,
                         std::map<int, NekDouble> &Dirs);
    void SolveFreeRigidBody2D(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                              const Array<OneD, NekDouble> &forcebody,
                              std::map<int, NekDouble> &Dirs);
    int m_index;
    NekDouble m_currentTime;
    NekDouble m_timestep;
    NekDouble m_beta;
    NekDouble m_gamma;
    NekDouble m_nonlinearTolerance;
    int m_nonlinearMaxIterations;
    int m_spacedim;
    bool m_isRoot;
    bool m_hasFreeMotion;
    bool m_hasRotation;
    // Fully prescribed motions in both 2D and 3D share the quaternion/MRF
    // state update.  Planar input occupies the z-rotation slots.
    bool m_prescribedMRF;
    bool m_free3D6DoF;
    bool m_hasCustomThetaConvention;
    int m_outputFrequency;
    std::ofstream m_outputStream;
    std::set<int> m_dirDoFs;
    std::set<int> m_inertialTransConstraints;
    std::set<int> m_bodyAngularConstraints;
    Array<OneD, NekDouble> m_inertialConstraintPosition;
    Array<OneD, NekDouble> m_inertialConstraintVelocity;
    Array<OneD, NekDouble> m_inertialConstraintAcceleration;
    Array<OneD, bool> m_hasInertialConstraintPosition;
    Array<OneD, int> m_thetaOrder;
    Array<OneD, bool> m_thetaBodyFrame;
    // position and velocity
    // m_vel[0][0,dim], displacement, dim translational + 1 rotation
    // m_vel[1][0,dim], velocity, dim translational + 1 rotation
    // m_vel[2][0,dim], acceleration, dim translational + 1 rotation
    Array<OneD, NekDouble> m_inertialPosition;
    Array<OneD, NekDouble> m_inertialVelocity;
    Array<OneD, NekDouble> m_inertialAcceleration;
    Array<OneD, NekDouble> m_quaternion;
    Array<OneD, Array<OneD, NekDouble>> m_vel;
    // externel force or moving velocity
    std::map<int, LibUtilities::EquationSharedPtr> m_extForceFunction;
    std::map<int, LibUtilities::EquationSharedPtr> m_frameVelFunction;
    Array<OneD, NekDouble> m_extForceXYZ;
    Array<OneD, NekDouble> m_pivot;
    // Body-frame vector from PIVOTPOINT to the centre of mass.
    Array<OneD, NekDouble> m_comOffset;
    NekDouble m_pivotdistance;
    // fluid force
    Array<OneD, NekDouble> m_oldFvis;
    // linear system
    NekDouble m_mass;
    NekDouble m_rotaionInertia;
    Array<OneD, NekDouble> m_rotationInertia;
    Array<OneD, NekDouble> m_M;
    Array<OneD, NekDouble> m_C;
    Array<OneD, NekDouble> m_K;
    // utility classes
    Newmark_BetaSolver m_bodySolver;
    FrameTransform m_frame;
};

} // namespace Nektar

#endif
