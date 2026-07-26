///////////////////////////////////////////////////////////////////////////////
//
// File: RigidSolver.cpp
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
// Description: Solving the absolute flow in a moving body frame,
// by adding (U0 + Omega X (x - x0)) . grad u - Omega X u
// as the body force.
// U0 is the translational velocity of the body frame.
// Omega is the angular velocity.
// x0 is the rotation pivot in the body frame.
// All vectors use the basis of the body frame.
// Translational motion is allowed for all dimensions.
// Rotation is not allowed for 1D, 2DH1D, 3DH2D.
// Rotation in z direction is allowed for 2D and 3DH1D.
// Rotation in 3 directions are allowed for 3D.
// Prescribed 3D rotation is represented internally using quaternions.
///////////////////////////////////////////////////////////////////////////////

#include <IncNavierStokesSolver/EquationSystems/RigidSolver.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/VmathArray.hpp>
#include <LibUtilities/LinearAlgebra/Lapack.hpp>
#include <MultiRegions/ExpList.h>
#include <SolverUtils/Core/MovingFrameTransforms.h>
#include <SolverUtils/Filters/FilterInterfaces.hpp>
#include <algorithm>
#include <array>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/format.hpp>

namespace Nektar
{

namespace
{

int AxisIndex(const char axis)
{
    switch (axis)
    {
        case 'x':
        case 'X':
            return 0;
        case 'y':
        case 'Y':
            return 1;
        case 'z':
        case 'Z':
            return 2;
        default:
            return -1;
    }
}

std::string AxisName(const int axis)
{
    static const std::string names[3] = {"x", "y", "z"};
    return names[axis];
}

bool IsThetaFrameBody(const std::string &frame)
{
    if (boost::iequals(frame, "Body") ||
        boost::iequals(frame, "BodyFrame"))
    {
        return true;
    }
    if (boost::iequals(frame, "Lab") ||
        boost::iequals(frame, "Inertial") ||
        boost::iequals(frame, "InertialFrame"))
    {
        return false;
    }

    ASSERTL0(false, "Theta frame must be Body, Lab or Inertial.");
    return false;
}

std::string GetRequiredElementText(const TiXmlElement *element,
                                   const std::string &tagName)
{
    const char *text = element->GetText();
    ASSERTL0(text, tagName + " must contain a value.");
    return text;
}

} // namespace

void RigidSolver::InitObject(const LibUtilities::SessionReaderSharedPtr session,
                             const MultiRegions::ExpListSharedPtr &pField,
                             Array<OneD, NekDouble> pivot)
{
    bool isH1d, isH2d;
    session->MatchSolverInfo("Homogeneous", "1D", isH1d, false);
    session->MatchSolverInfo("Homogeneous", "2D", isH2d, false);
    bool singleMode, halfMode;
    session->MatchSolverInfo("ModeType", "SingleMode", singleMode, false);
    session->MatchSolverInfo("ModeType", "HalfMode", halfMode, false);
    if (singleMode || halfMode)
    {
        isH1d = false;
    }
    int expdim            = isH2d ? 1 : pField->GetGraph()->GetMeshDimension();
    m_spacedim            = expdim + (isH1d ? 1 : 0) + (isH2d ? 2 : 0);
    m_isRoot              = pField->GetComm()->TreatAsRankZero();
    m_index               = 1;
    m_currentTime              = -1.;
    m_prescribed3DMRF          = false;
    m_free3D6DoF               = false;
    // The Newton-consistent path is the default for an unconstrained 2D body.
    // Set UseUnifiedFreeRigidBody = 0 only to reproduce legacy results.
    m_useUnifiedFreeRigidBody  = true;
    m_inertialTransConstraints.clear();
    m_bodyAngularConstraints.clear();
    m_hasCustomThetaConvention = false;
    m_thetaOrder               = Array<OneD, int>(3, 0);
    m_thetaBodyFrame           = Array<OneD, bool>(3, false);
    m_thetaOrder[0]            = 2;
    m_thetaOrder[1]            = 1;
    m_thetaOrder[2]            = 0;
    m_extForceXYZ              = Array<OneD, NekDouble>(6, 0.0);
    m_oldFvis                  = Array<OneD, NekDouble>(6, 0.0);
    m_quaternion = SolverUtils::MovingFrame::IdentityQuaternion();
    TiXmlElement *pSolver      = session->GetElement("Nektar/RIGIDSOLVER");
    LoadParameters(session, pSolver);
    InitBodySolver(session, pSolver, pivot);
    CheckParameters();
}

void RigidSolver::CheckParameters()
{
    if (m_free3D6DoF)
    {
        ASSERTL0(m_spacedim == 3,
                 "Full free rigid-body motion is available only in 3D.");
        ASSERTL0(fabs(m_pivotdistance) < NekConstants::kNekZeroTol,
                 "Full free 3D rigid-body motion currently requires "
                 "PIVOTDISTANCE = 0 (the pivot is the centre of mass).");
        NekDouble restoringTerms = 0.0;
        for (size_t i = 0; i < m_C.size(); ++i)
        {
            restoringTerms += fabs(m_C[i]) + fabs(m_K[i]);
        }
        ASSERTL0(restoringTerms < NekConstants::kNekZeroTol,
                 "Full free 3D 6DoF motion currently requires zero DAMPING "
                 "and RIGIDITY.");
        m_inertialPosition = Array<OneD, NekDouble>(3, 0.0);
        m_inertialVelocity = Array<OneD, NekDouble>(3, 0.0);
        m_inertialAcceleration = Array<OneD, NekDouble>(3, 0.0);
        m_inertialConstraintPosition = Array<OneD, NekDouble>(3, 0.0);
        m_inertialConstraintVelocity = Array<OneD, NekDouble>(3, 0.0);
        m_inertialConstraintAcceleration = Array<OneD, NekDouble>(3, 0.0);
        m_hasInertialConstraintPosition = Array<OneD, bool>(3, false);
        for (const int direction : m_inertialTransConstraints)
        {
            m_hasInertialConstraintPosition[direction] =
                m_frameVelFunction.find(direction + 6) !=
                m_frameVelFunction.end();
        }
        m_solveType = eFree3D6DoF;
        return;
    }
    m_prescribed3DMRF =
        m_spacedim == 3 && !m_hasFreeMotion && HasFull3DPrescribedOrientation();
    if (m_prescribed3DMRF)
    {
        m_solveType = ePrescribed3DMRF;
        return;
    }

    // Count free translational DoFs
    int nFreeTrans = 0;
    bool freeX     = false;

    for (int i = 0; i < m_spacedim; ++i)
    {
        if (m_dirDoFs.find(i) == m_dirDoFs.end())
        {
            ++nFreeTrans;
            if (i == 0)
            {
                freeX = true;
            }
        }
    }
    const bool rotationFree =
        m_hasRotation && m_dirDoFs.find(m_spacedim) == m_dirDoFs.end();
    if (m_spacedim == 2 && rotationFree)
    {
        // With a free rotation, translational constraints are imposed in the
        // inertial frame. Keep the body-frame components as unknowns and add
        // their inertial-frame velocity constraints to the Newton system.
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_dirDoFs.find(i) != m_dirDoFs.end())
            {
                m_inertialTransConstraints.insert(i);
                m_dirDoFs.erase(i);
            }
        }
    }
    if (!m_hasRotation)
    {
        m_solveType = 0; // inertial frame
    }
    else if (rotationFree)
    {
        // The nonlinear body-frame solver supports any subset of the two
        // translational DoFs when Omega_z is free.
        m_solveType = 2;
    }
    else if (nFreeTrans == 1 && freeX)
    {
        m_solveType = 1; // with given rotation and only x free
        if (m_spacedim >= 2)
        {
            m_dirDoFs.erase(1);
        }
    }
    else if ((nFreeTrans == 2 || nFreeTrans == 0) && m_hasRotation)
    {
        m_solveType = 2; // full planar rigid solver (rotation involved)
    }
    else
    {
        ASSERTL0(false, "Unsupported rigid motion: this solver supports only "
                        "(i) translation-only motion, "
                        "(ii) x-only translation with optional rotation, or "
                        "(iii) full planar rigid-body motion.");
    }
    if (m_hasRotation)
    {
        m_inertialPosition = Array<OneD, NekDouble>(2, 0.);
        m_inertialVelocity = Array<OneD, NekDouble>(2, 0.);
        m_inertialAcceleration = Array<OneD, NekDouble>(2, 0.);
        m_inertialConstraintPosition = Array<OneD, NekDouble>(2, 0.);
        m_inertialConstraintVelocity = Array<OneD, NekDouble>(2, 0.);
        m_inertialConstraintAcceleration = Array<OneD, NekDouble>(2, 0.);
        m_hasInertialConstraintPosition = Array<OneD, bool>(2, false);
        for (const int direction : m_inertialTransConstraints)
        {
            m_hasInertialConstraintPosition[direction] =
                m_frameVelFunction.find(direction + 6) !=
                m_frameVelFunction.end();
        }
        // m_K should be zero
        NekDouble sum = 0;
        for (size_t i = 0; i < m_K.size(); ++i)
        {
            sum += fabs(m_K[i]);
        }
        for (size_t i = 0; i < m_C.size(); ++i)
        {
            sum += fabs(m_C[i]);
        }
        ASSERTL0(sum < NekConstants::kNekZeroTol,
                 "C, K should be zero [in body frame solver].");
    }
}

bool RigidSolver::HasFull3DPrescribedOrientation() const
{
    if (m_spacedim != 3)
    {
        return false;
    }

    return m_frameVelFunction.find(3) != m_frameVelFunction.end() ||
           m_frameVelFunction.find(4) != m_frameVelFunction.end() ||
           m_frameVelFunction.find(9) != m_frameVelFunction.end() ||
           m_frameVelFunction.find(10) != m_frameVelFunction.end() ||
           m_frameVelFunction.find(15) != m_frameVelFunction.end() ||
           m_frameVelFunction.find(16) != m_frameVelFunction.end() ||
           m_frameVelFunction.find(
               SolverUtils::MovingFrame::kQuaternionOffset) !=
               m_frameVelFunction.end() ||
           m_frameVelFunction.find(
               SolverUtils::MovingFrame::kQuaternionOffset + 1) !=
               m_frameVelFunction.end() ||
           m_frameVelFunction.find(
               SolverUtils::MovingFrame::kQuaternionOffset + 2) !=
               m_frameVelFunction.end() ||
           m_frameVelFunction.find(
               SolverUtils::MovingFrame::kQuaternionOffset + 3) !=
               m_frameVelFunction.end();
}

void RigidSolver::SetMovableDoFs(std::vector<bool> &moveDoFs)
{
    if (m_free3D6DoF)
    {
        ASSERTL0(moveDoFs.size() >= 6,
                 "Full free 3D motion requires six movable-DoF entries.");
        for (int i = 0; i < 6; ++i)
        {
            moveDoFs[i] = true;
        }
        return;
    }

    if (m_prescribed3DMRF)
    {
        const int nDoFs = std::min(6, static_cast<int>(moveDoFs.size()));
        for (int i = 0; i < nDoFs; ++i)
        {
            moveDoFs[i] = true;
        }
        return;
    }

    for (const auto &it : m_frameVelFunction)
    {
        if (it.first < 6)
        {
            moveDoFs[it.first] = true;
        }
    }
    for (int i = 0; i < m_spacedim; ++i)
    {
        if (m_dirDoFs.find(i) == m_dirDoFs.end() &&
            m_inertialTransConstraints.find(i) ==
                m_inertialTransConstraints.end())
        {
            moveDoFs[i] = true;
        }
    }
    if (m_dirDoFs.find(m_spacedim) == m_dirDoFs.end())
    {
        moveDoFs[5] = true;
    }
    if (moveDoFs[5])
    {
        // since m_frameVelFunction is defined in inerital frame
        // motions in 2 directions occur if there is rotation
        if (moveDoFs[0] || moveDoFs[1])
        {
            moveDoFs[0] = true;
            moveDoFs[1] = true;
        }
    }
}

NekDouble RigidSolver::EvaluateExpression(
    const LibUtilities::SessionReaderSharedPtr session, std::string expression)
{
    NekDouble value = 0.;
    try
    {
        LibUtilities::Equation expession(session->GetInterpreter(), expression);
        value = expession.Evaluate();
    }
    catch (const std::runtime_error &)
    {
        NEKERROR(ErrorUtil::efatal, "Error evaluating expression" + expression);
    }
    return value;
}

void RigidSolver::ParserFunctionToMap(
    const bool allowzero, const LibUtilities::SessionReaderSharedPtr session,
    const std::string &FuncName, const std::string &var,
    std::map<int, LibUtilities::EquationSharedPtr> &result, const int i)
{
    if (session->DefinesFunction(FuncName, var))
    {
        LibUtilities::EquationSharedPtr equ =
            session->GetFunction(FuncName, var);
        if (allowzero || "0" != equ->GetExpression())
        {
            result[i] = equ;
        }
    }
}

void RigidSolver::LoadParameters(
    const LibUtilities::SessionReaderSharedPtr session,
    const TiXmlElement *pSolver)
{
    const TiXmlElement *funcNameElmt;
    // load frame velocity
    funcNameElmt = pSolver->FirstChildElement("FRAMEVELOCITY");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        ASSERTL0(session->DefinesFunction(FuncName),
                 "Function '" + FuncName + "' is not defined in the session.");
        // linear velocity
        for (int i = 0; i < m_spacedim; ++i)
        {
            std::string var = session->GetVariable(i);
            ParserFunctionToMap(false, session, FuncName, var,
                                m_frameVelFunction, i);
        }
        // linear displacement and acceleration
        std::vector<std::string> linearDispVar = {"X", "Y", "Z"};
        std::vector<std::string> linearAcceVar = {"A_x", "A_y", "A_z"};
        for (int i = 0; i < m_spacedim; ++i)
        {
            ParserFunctionToMap(true, session, FuncName, linearDispVar[i],
                                m_frameVelFunction, i + 6);
            ParserFunctionToMap(true, session, FuncName, linearAcceVar[i],
                                m_frameVelFunction, i + 12);
        }
        // angular velocities
        std::vector<std::string> angularVar = {"Omega_x", "Omega_y", "Omega_z"};
        m_hasRotation                       = false;
        for (int i = 0; i < 3; ++i)
        {
            std::string var = angularVar[i];
            ParserFunctionToMap(false, session, FuncName, var,
                                m_frameVelFunction, i + 3);
            m_hasRotation = m_hasRotation || m_frameVelFunction.find(i + 3) !=
                                                 m_frameVelFunction.end();
        }
        // angular displacement and acceleration
        std::vector<std::string> angularDispVar  = {"Theta_x", "Theta_y",
                                                    "Theta_z"};
        std::vector<std::string> angularAccepVar = {"DOmega_x", "DOmega_y",
                                                    "DOmega_z"};
        for (int i = 0; i < 3; ++i)
        {
            ParserFunctionToMap(true, session, FuncName, angularDispVar[i],
                                m_frameVelFunction, i + 3 + 6);
            ParserFunctionToMap(true, session, FuncName, angularAccepVar[i],
                                m_frameVelFunction, i + 3 + 12);
        }
        std::vector<std::string> quatVar = {"Q0", "Q1", "Q2", "Q3"};
        for (int i = 0; i < 4; ++i)
        {
            ParserFunctionToMap(true, session, FuncName, quatVar[i],
                                m_frameVelFunction,
                                SolverUtils::MovingFrame::kQuaternionOffset +
                                    i);
        }
    }

    LoadThetaConvention(pSolver);

    // load external force
    funcNameElmt = pSolver->FirstChildElement("EXTERNALFORCE");
    if (funcNameElmt)
    {
        std::string FuncName = funcNameElmt->GetText();
        if (session->DefinesFunction(FuncName))
        {
            std::vector<std::string> forceVar = {"Fx", "Fy", "Fz"};
            for (int i = 0; i < m_spacedim; ++i)
            {
                std::string var = forceVar[i];
                ParserFunctionToMap(false, session, FuncName, var,
                                    m_extForceFunction, i);
            }
            std::vector<std::string> momentVar = {"Mx", "My", "Mz"};
            for (int i = 0; i < 3; ++i)
            {
                std::string var = momentVar[i];
                ParserFunctionToMap(false, session, FuncName, var,
                                    m_extForceFunction, i + 3);
            }
        }
    }

    // OutputFile
    const TiXmlElement *mssgTag = pSolver->FirstChildElement("OutputFile");
    std::string filename;
    if (mssgTag)
    {
        filename = mssgTag->GetText();
    }
    else
    {
        filename = session->GetSessionName();
    }
    if (!(filename.length() >= 4 &&
          filename.substr(filename.length() - 4) == ".mrf"))
    {
        filename += ".mrf";
    }
    if (m_isRoot)
    {
        m_outputStream.open(filename.c_str());
        if (m_spacedim == 2)
        {
            m_outputStream
                << "Variables = t, x, ux, ax, y, uy, ay, theta, omega, domega"
                << std::endl;
        }
        else if (m_spacedim == 3)
        {
            const TiXmlElement *motion =
                pSolver->FirstChildElement("MOTIONPRESCRIBED");
            std::vector<std::string> motionValues;
            if (motion)
            {
                ParseUtils::GenerateVector(motion->GetText(), motionValues);
            }
            if (HasFull3DPrescribedOrientation() || motionValues.size() == 6)
            {
                m_outputStream
                    << "Variables = t, x, ux, ax, y, uy, ay, z, uz, az, "
                       "theta_x, omega_x, domega_x, theta_y, omega_y, "
                       "domega_y, theta_z, omega_z, domega_z"
                    << std::endl;
            }
            else
            {
                m_outputStream << "Variables = t, x, ux, ax, y, uy, ay, z, "
                                  "uz, az, theta, omega, domega"
                               << std::endl;
            }
        }
    }

    // output frequency
    m_outputFrequency = 1;
    mssgTag           = pSolver->FirstChildElement("OutputFrequency");
    if (mssgTag)
    {
        std::vector<std::string> values;
        std::string mssgStr = mssgTag->GetText();
        m_outputFrequency   = round(EvaluateExpression(session, mssgStr));
    }
    ASSERTL0(m_outputFrequency > 0,
             "OutputFrequency should be greater than zero.");
}

void RigidSolver::LoadThetaConvention(const TiXmlElement *pSolver)
{
    const TiXmlElement *thetaOrder = pSolver->FirstChildElement("ThetaOrder");
    if (thetaOrder)
    {
        std::string order = GetRequiredElementText(thetaOrder, "ThetaOrder");
        bool seen[3]     = {false, false, false};
        int cnt          = 0;

        for (size_t i = 0; i < order.size(); ++i)
        {
            const char c = order[i];
            if (c == ',' || c == ' ' || c == '\t')
            {
                continue;
            }

            const int axis = AxisIndex(c);
            ASSERTL0(axis >= 0,
                     "ThetaOrder must contain only x, y and z axes.");
            ASSERTL0(!seen[axis],
                     "ThetaOrder must not repeat the same axis.");
            ASSERTL0(cnt < 3, "ThetaOrder must contain exactly three axes.");

            m_thetaOrder[cnt++] = axis;
            seen[axis]          = true;
        }

        ASSERTL0(cnt == 3, "ThetaOrder must contain exactly three axes.");
        m_hasCustomThetaConvention = true;
    }

    const TiXmlElement *thetaFrame = pSolver->FirstChildElement("ThetaFrame");
    if (thetaFrame)
    {
        const bool isBody = IsThetaFrameBody(
            GetRequiredElementText(thetaFrame, "ThetaFrame"));
        for (int i = 0; i < 3; ++i)
        {
            m_thetaBodyFrame[i] = isBody;
        }
        m_hasCustomThetaConvention = true;
    }

    for (int i = 0; i < 3; ++i)
    {
        const std::string tagName = "Theta_" + AxisName(i) + "Frame";
        const TiXmlElement *axisFrame = pSolver->FirstChildElement(tagName);
        if (axisFrame)
        {
            m_thetaBodyFrame[i] = IsThetaFrameBody(
                GetRequiredElementText(axisFrame, tagName));
            m_hasCustomThetaConvention = true;
        }
    }
}

Array<OneD, NekDouble> RigidSolver::QuaternionFromConfiguredTheta(
    const Array<OneD, NekDouble> &theta) const
{
    using namespace SolverUtils::MovingFrame;

    if (!m_hasCustomThetaConvention)
    {
        return QuaternionFromEulerZYX(theta[0], theta[1], theta[2]);
    }

    Array<OneD, NekDouble> q = IdentityQuaternion();
    for (int i = 2; i >= 0; --i)
    {
        const int axis = m_thetaOrder[i];
        const Array<OneD, NekDouble> qAxis =
            QuaternionFromAxisAngle(axis, theta[axis]);

        if (m_thetaBodyFrame[axis])
        {
            q = MultiplyQuaternions(q, qAxis);
        }
        else
        {
            q = MultiplyQuaternions(qAxis, q);
        }
    }
    return q;
}

void RigidSolver::InitBodySolver(
    const LibUtilities::SessionReaderSharedPtr session,
    const TiXmlElement *pSolver, Array<OneD, NekDouble> pivot)
{
    int NumDof = m_spacedim + 1;
    const TiXmlElement *mssgTag;
    std::string mssgStr;
    std::vector<std::string> prescribedValues;
    mssgTag = pSolver->FirstChildElement("MOTIONPRESCRIBED");
    if (mssgTag)
    {
        ParseUtils::GenerateVector(mssgTag->GetText(), prescribedValues);
        if (m_spacedim == 3 && prescribedValues.size() == 6)
        {
            m_free3D6DoF = true;
            NumDof = 6;
        }
    }
    // allocate memory and initialise
    m_vel = Array<OneD, Array<OneD, NekDouble>>(3);
    for (size_t i = 0; i < 3; ++i)
    {
        m_vel[i] = Array<OneD, NekDouble>(NumDof, 0.);
    }
    // read prescribed motion DoFs
    m_dirDoFs.clear();
    for (int i = 0; i < NumDof; ++i)
    {
        m_dirDoFs.insert(i);
    }
    mssgTag = pSolver->FirstChildElement("MOTIONPRESCRIBED");
    if (mssgTag)
    {
        const std::vector<std::string> &values = prescribedValues;
        const bool full3DPrescribedInput =
            m_spacedim == 3 && HasFull3DPrescribedOrientation() &&
            values.size() == 6;
        ASSERTL0(values.size() == NumDof || full3DPrescribedInput,
                 "MOTIONPRESCRIBED vector should be of size " +
                     std::to_string(NumDof) +
                     (m_spacedim == 3 ? " or 6 for prescribed 3D MRF" : ""));
        if (m_free3D6DoF)
        {
            for (int i = 0; i < NumDof; ++i)
            {
                m_dirDoFs.erase(i);
                if (i < 3 && EvaluateExpression(session, values[i]) != 0.0)
                {
                    m_inertialTransConstraints.insert(i);
                }
                if (i >= 3 && EvaluateExpression(session, values[i]) != 0.0)
                {
                    m_bodyAngularConstraints.insert(i - 3);
                }
            }
        }
        else if (full3DPrescribedInput)
        {
            for (int i = 0; i < 6; ++i)
            {
                ASSERTL0(EvaluateExpression(session, values[i]) != 0,
                         "Full prescribed 3D MRF requires all six "
                         "MOTIONPRESCRIBED entries to be prescribed.");
            }
        }
        else
        {
            for (int i = 0; i < NumDof; ++i)
            {
                if (EvaluateExpression(session, values[i]) == 0)
                {
                    m_dirDoFs.erase(i);
                }
            }
        }
    }
    m_hasRotation = m_hasRotation || m_free3D6DoF ||
                    m_dirDoFs.find(m_spacedim) == m_dirDoFs.end();
    m_hasFreeMotion = m_dirDoFs.size() < NumDof;
    // read mass matrix
    m_M     = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag = pSolver->FirstChildElement("MASS");
    ASSERTL0(m_dirDoFs.size() == NumDof || mssgTag, "Mass matrix is required.");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == 1, "Mass matrix should be a scalar.");
        m_mass = EvaluateExpression(session, values[0]);
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_M[i + i * NumDof] = m_mass;
        }
    }
    mssgTag = pSolver->FirstChildElement("ROTATIONINERTIA");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == 1 || (m_free3D6DoF && values.size() == 3),
                 "Inertia should be a scalar, or three principal inertias for "
                 "full free 3D 6DoF motion.");
        m_rotaionInertia = EvaluateExpression(session, values[0]);
        if (m_free3D6DoF)
        {
            m_rotationInertia = Array<OneD, NekDouble>(3, 0.0);
            for (int i = 0; i < 3; ++i)
            {
                m_rotationInertia[i] = EvaluateExpression(
                    session, values[values.size() == 1 ? 0 : i]);
                ASSERTL0(m_rotationInertia[i] > 0.0,
                         "All principal moments of inertia must be positive.");
                m_M[(i + 3) + (i + 3) * NumDof] = m_rotationInertia[i];
            }
        }
        else
        {
            m_M[NumDof * NumDof - 1] = m_rotaionInertia;
        }
    }
    // read damping matrix
    m_C     = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag = pSolver->FirstChildElement("DAMPING");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof * NumDof,
                 "Damping matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        int count = 0;
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                m_C[count] = EvaluateExpression(session, values[count]);
                ++count;
            }
        }
    }
    // read rigidity matrix
    m_K     = Array<OneD, NekDouble>(NumDof * NumDof, 0.);
    mssgTag = pSolver->FirstChildElement("RIGIDITY");
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == NumDof * NumDof,
                 "Rigidity matrix should be of size " + std::to_string(NumDof) +
                     "X" + std::to_string(NumDof));
        int count = 0;
        for (int i = 0; i < NumDof; ++i)
        {
            for (int j = 0; j < NumDof; ++j)
            {
                m_K[count] = EvaluateExpression(session, values[count]);
                ++count;
            }
        }
    }
    // read pivot point
    mssgTag = pSolver->FirstChildElement("PIVOTPOINT");
    m_pivot = Array<OneD, NekDouble>(m_spacedim, 0.);
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == m_spacedim,
                 "PIVOTPOINT vector should be of size " +
                     std::to_string(m_spacedim));
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_pivot[i] = EvaluateExpression(session, values[i]);
        }
    }
    Vmath::Vcopy(m_spacedim, m_pivot, 1, pivot, 1);
    // read the distance between pivotpoint and masscenter
    mssgTag         = pSolver->FirstChildElement("PIVOTDISTANCE");
    m_pivotdistance = 0.;
    if (mssgTag)
    {
        std::vector<std::string> values;
        mssgStr = mssgTag->GetText();
        ParseUtils::GenerateVector(mssgStr, values);
        ASSERTL0(values.size() == 1, "PivotDistance should be a scalar.");
        m_pivotdistance = EvaluateExpression(session, values[0]);
    }
    // read Newmark Beta paramters
    m_timestep = session->GetParameter("TimeStep");
    m_beta     = 0.25;
    m_gamma    = 0.51;
    if (session->DefinesParameter("NewmarkBeta"))
    {
        m_beta = session->GetParameter("NewmarkBeta");
    }
    if (session->DefinesParameter("NewmarkGamma"))
    {
        m_gamma = session->GetParameter("NewmarkGamma");
    }
    m_nonlinearTolerance    = 1.0e-10;
    m_nonlinearMaxIterations = 8;
    if (session->DefinesParameter("RigidBodyNonlinearTolerance"))
    {
        m_nonlinearTolerance =
            session->GetParameter("RigidBodyNonlinearTolerance");
    }
    if (session->DefinesParameter("RigidBodyNonlinearMaxIterations"))
    {
        m_nonlinearMaxIterations = static_cast<int>(
            session->GetParameter("RigidBodyNonlinearMaxIterations"));
    }
    ASSERTL0(m_nonlinearTolerance > 0.0,
             "RigidBodyNonlinearTolerance must be positive.");
    ASSERTL0(m_nonlinearMaxIterations > 0,
             "RigidBodyNonlinearMaxIterations must be positive.");
    if (session->DefinesParameter("UseUnifiedFreeRigidBody"))
    {
        m_useUnifiedFreeRigidBody =
            session->GetParameter("UseUnifiedFreeRigidBody") != 0.0;
    }
}

void RigidSolver::UpdatePrescribed(const NekDouble &time,
                                   std::map<int, NekDouble> &Dirs)
{
    int NumDof = m_spacedim + 1;
    for (auto it : m_frameVelFunction)
    {
        if (it.first < 3)
        {
            const NekDouble value = it.second->Evaluate(0., 0., 0., time);
            if (it.first < m_spacedim &&
                m_inertialTransConstraints.find(it.first) !=
                    m_inertialTransConstraints.end())
            {
                m_inertialConstraintVelocity[it.first] = value;
            }
            else
            {
                Dirs[it.first] = value;
            }
        }
        else if (it.first == 5)
        {
            Dirs[m_spacedim] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first < 9)
        {
            const int direction = it.first - 6;
            const NekDouble value = it.second->Evaluate(0., 0., 0., time);
            if (direction < m_spacedim &&
                m_inertialTransConstraints.find(direction) !=
                    m_inertialTransConstraints.end())
            {
                m_inertialConstraintPosition[direction] = value;
            }
            else
            {
                Dirs[NumDof + direction] = value;
            }
        }
        else if (it.first == 11)
        {
            Dirs[NumDof + m_spacedim] = it.second->Evaluate(0., 0., 0., time);
        }
        else if (it.first < 15)
        {
            const int direction = it.first - 12;
            const NekDouble value = it.second->Evaluate(0., 0., 0., time);
            if (direction < m_spacedim &&
                m_inertialTransConstraints.find(direction) !=
                    m_inertialTransConstraints.end())
            {
                m_inertialConstraintAcceleration[direction] = value;
            }
            else
            {
                Dirs[(NumDof << 1) + direction] = value;
            }
        }
        else if (it.first == 17)
        {
            Dirs[(NumDof << 1) + m_spacedim] =
                it.second->Evaluate(0., 0., 0., time);
        }
    }
    for (auto i : m_dirDoFs)
    {
        if (Dirs.find(i) == Dirs.end())
        {
            Dirs[i]                 = 0.;
            Dirs[i + NumDof]        = 0.;
            Dirs[i + (NumDof << 1)] = 0.;
        }
    }
}

void RigidSolver::SetOldFvis(Array<OneD, NekDouble> force)
{
    for (int i = 0; i < 6; ++i)
    {
        m_oldFvis[i] = force[6 + i];
    }
}

/**
 * @brief Updates the forcing array with the current required forcing.
 * @param pFields
 * @param time
 */
void RigidSolver::UpdateFrameVelocity(Array<OneD, NekDouble> &aeroforce,
                                      const NekDouble &time,
                                      Array<OneD, NekDouble> &MRFData)
{
    if (m_currentTime >= time)
    {
        return;
    }

    if (m_prescribed3DMRF)
    {
        UpdatePrescribedMRFData(time, MRFData);
        if (m_isRoot && m_index % m_outputFrequency == 0)
        {
            WritePrescribedMRFOutput(time, MRFData);
        }
        m_currentTime = time;
        ++m_index;
        return;
    }

    for (int i = 0; i < 6; ++i)
    {
        aeroforce[i] +=
            2. * aeroforce[6 + i] - m_oldFvis[i]; // extrapolate viscous
        m_oldFvis[i] = aeroforce[6 + i];
    }
    m_currentTime = time;
    std::map<int, NekDouble> Dirs;
    UpdatePrescribed(time, Dirs);
    // compute the velocites whoes functions are provided in inertial frame
    if (m_hasFreeMotion)
    {
        for (auto it : m_extForceFunction)
        {
            m_extForceXYZ[it.first] = it.second->Evaluate(0., 0., 0., time);
        }
    }
    SolveBodyMotion(m_vel, aeroforce, Dirs);
    if (m_free3D6DoF)
    {
        UpdateFree3DMRFData(MRFData);
        if (m_isRoot && m_index % m_outputFrequency == 0)
        {
            WriteFree3DMRFOutput(time, MRFData);
        }
        ++m_index;
        return;
    }
    Array<OneD, Array<OneD, NekDouble>> tmpVel;
    if (m_hasRotation)
    {
        // update inertial position of the body
        tmpVel = Array<OneD, Array<OneD, NekDouble>>(m_vel.size());
        for (size_t i = 0; i < m_vel.size(); ++i)
        {
            tmpVel[i] = Array<OneD, NekDouble>(m_vel[i].size());
            Vmath::Vcopy(m_vel[i].size(), m_vel[i], 1, tmpVel[i], 1);
        }
        Array<OneD, NekDouble> angle(3, 0.);
        angle[2] = tmpVel[0][m_spacedim];
        m_frame.SetAngle(angle);
        m_frame.BodyToInerital(m_spacedim, tmpVel[1], tmpVel[1]);
        m_frame.BodyToInerital(m_spacedim, tmpVel[2], tmpVel[2]);
        tmpVel[2][0] -= tmpVel[1][m_spacedim] * tmpVel[1][1];
        tmpVel[2][1] += tmpVel[1][m_spacedim] * tmpVel[1][0];
        for (int i = 0; i < m_spacedim; ++i)
        {
            if (m_inertialTransConstraints.find(i) !=
                m_inertialTransConstraints.end())
            {
                // The constraint is defined in the inertial frame. A supplied
                // displacement takes precedence; otherwise advance its
                // position with the prescribed velocity/acceleration.
                if (m_hasInertialConstraintPosition[i])
                {
                    m_inertialPosition[i] =
                        m_inertialConstraintPosition[i];
                }
                else
                {
                    m_inertialPosition[i] +=
                        m_timestep * m_inertialVelocity[i] +
                        m_timestep * m_timestep *
                            ((0.5 - m_beta) * m_inertialAcceleration[i] +
                             m_beta * m_inertialConstraintAcceleration[i]);
                }
                tmpVel[1][i] = m_inertialConstraintVelocity[i];
                tmpVel[2][i] = m_inertialConstraintAcceleration[i];
                m_inertialVelocity[i]     = tmpVel[1][i];
                m_inertialAcceleration[i] = tmpVel[2][i];
                continue;
            }
            m_inertialPosition[i] +=
                m_timestep * m_inertialVelocity[i] +
                m_timestep * m_timestep *
                    ((0.5 - m_beta) * m_inertialAcceleration[i] +
                     m_beta * tmpVel[2][i]);
            m_inertialVelocity[i]     = tmpVel[1][i];
            m_inertialAcceleration[i] = tmpVel[2][i];
        }
        tmpVel[0][0] = m_inertialPosition[0];
        tmpVel[0][1] = m_inertialPosition[1];
    }
    if (m_isRoot && m_index % m_outputFrequency == 0)
    {
        if (!m_hasRotation)
        {
            tmpVel = m_vel;
        }
        m_outputStream << boost::format("%25.19e") % time << " ";
        for (size_t i = 0; i < tmpVel[0].size(); ++i)
        {
            m_outputStream << boost::format("%25.19e") % tmpVel[0][i] << " "
                           << boost::format("%25.19e") % tmpVel[1][i] << " "
                           << boost::format("%25.19e") % tmpVel[2][i] << " ";
        }
        m_outputStream << std::endl;
    }
    // set displacements at the current time step
    UpdateMRFData(MRFData);
    ++m_index;
}

void RigidSolver::UpdatePrescribedMRFData(const NekDouble &time,
                                          Array<OneD, NekDouble> &MRFData)
{
    using namespace SolverUtils::MovingFrame;

    const auto evaluate = [&](const int idx, const NekDouble defaultValue) {
        auto it = m_frameVelFunction.find(idx);
        return it == m_frameVelFunction.end()
                   ? defaultValue
                   : it->second->Evaluate(0., 0., 0., time);
    };
    const auto hasFunction = [&](const int idx) {
        return m_frameVelFunction.find(idx) != m_frameVelFunction.end();
    };

    const NekDouble dt = m_currentTime >= 0.0 ? time - m_currentTime : 0.0;

    Array<OneD, NekDouble> disp(3, 0.0);
    Array<OneD, NekDouble> velInertial(3, 0.0);
    Array<OneD, NekDouble> omegaInertial(3, 0.0);
    Array<OneD, NekDouble> accInertial(3, 0.0);
    Array<OneD, NekDouble> alphaInertial(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        disp[i] = MRFData[i];
    }

    for (int i = 0; i < 3; ++i)
    {
        velInertial[i]   = evaluate(i, 0.0);
        omegaInertial[i] = evaluate(i + 3, 0.0);
        accInertial[i]   = evaluate(i + 12, 0.0);
        alphaInertial[i] = evaluate(i + 15, 0.0);

        if (hasFunction(i + 6))
        {
            disp[i] = evaluate(i + 6, disp[i]);
        }
        else if (dt > 0.0)
        {
            disp[i] += dt * (velInertial[i] + 0.5 * dt * accInertial[i]);
        }
    }

    const bool hasQuaternion =
        hasFunction(kQuaternionOffset) || hasFunction(kQuaternionOffset + 1) ||
        hasFunction(kQuaternionOffset + 2) ||
        hasFunction(kQuaternionOffset + 3);
    const bool hasEuler = hasFunction(9) || hasFunction(10) || hasFunction(11);
    const bool hasOmega = hasFunction(3) || hasFunction(4) || hasFunction(5);
    const bool hasAlpha = hasFunction(15) || hasFunction(16) || hasFunction(17);

    // XML angular data are lab-frame inputs; MRFData stores body-frame
    // components derived through the quaternion state.
    const Array<OneD, NekDouble> qOld = QuaternionFromFrameData(MRFData);
    Array<OneD, NekDouble> q          = qOld;
    if (hasQuaternion)
    {
        q = IdentityQuaternion();
        for (int i = 0; i < 4; ++i)
        {
            q[i] = evaluate(kQuaternionOffset + i, q[i]);
        }
        q = NormalizeQuaternion(q);
    }
    else if (hasEuler)
    {
        Array<OneD, NekDouble> thetaInput(3, 0.0);
        thetaInput[0] = evaluate(9, 0.0);
        thetaInput[1] = evaluate(10, 0.0);
        thetaInput[2] = evaluate(11, 0.0);
        q             = QuaternionFromConfiguredTheta(thetaInput);
    }
    else if (dt > 0.0 && hasOmega)
    {
        const Array<OneD, NekDouble> omegaBodyForStep =
            RotateInertialToBody(q, omegaInertial);
        q = IntegrateQuaternionBodyOmega(q, omegaBodyForStep, dt);
    }
    q = MakeQuaternionContinuous(qOld, q);

    const Array<OneD, NekDouble> theta = EulerZYXFromQuaternion(q);
    const Array<OneD, NekDouble> velBody =
        RotateInertialToBody(q, velInertial);
    const Array<OneD, NekDouble> accBodyInertial =
        RotateInertialToBody(q, accInertial);
    Array<OneD, NekDouble> omegaBody(3, 0.0);
    Array<OneD, NekDouble> qdot(4, 0.0);
    if (hasOmega)
    {
        omegaBody = RotateInertialToBody(q, omegaInertial);
        qdot      = QuaternionDerivativeFromBodyOmega(q, omegaBody);
    }
    else if (dt > 0.0)
    {
        for (int i = 0; i < 4; ++i)
        {
            qdot[i] = (q[i] - qOld[i]) / dt;
        }
        omegaBody = BodyOmegaFromQuaternionDerivative(q, qdot);
    }

    Array<OneD, NekDouble> alphaBody(3, 0.0);
    if (hasAlpha)
    {
        alphaBody = RotateInertialToBody(q, alphaInertial);
    }
    else if (dt > 0.0)
    {
        Array<OneD, NekDouble> omegaBodyOld(3, 0.0);
        for (int i = 0; i < 3; ++i)
        {
            omegaBodyOld[i] = MRFData[i + 9];
        }
        for (int i = 0; i < 3; ++i)
        {
            alphaBody[i] = (omegaBody[i] - omegaBodyOld[i]) / dt;
        }
    }

    Array<OneD, NekDouble> omegaCrossVel(3, 0.0);
    Cross(omegaBody, velBody, omegaCrossVel);

    for (int i = 0; i < 3; ++i)
    {
        MRFData[i]      = disp[i];
        MRFData[i + 3]  = theta[i];
        MRFData[i + 6]  = velBody[i];
        MRFData[i + 9]  = omegaBody[i];
        MRFData[i + 12] = accBodyInertial[i] - omegaCrossVel[i];
        MRFData[i + 15] = alphaBody[i];
    }
    StoreQuaternionInFrameData(q, MRFData);
}

void RigidSolver::WritePrescribedMRFOutput(
    const NekDouble &time, const Array<OneD, NekDouble> &MRFData)
{
    using namespace SolverUtils::MovingFrame;

    const Array<OneD, NekDouble> q = QuaternionFromFrameData(MRFData);
    Array<OneD, NekDouble> velBody(3, 0.0);
    Array<OneD, NekDouble> accBodyFrame(3, 0.0);
    Array<OneD, NekDouble> omegaBody(3, 0.0);
    Array<OneD, NekDouble> alphaBody(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        velBody[i]      = MRFData[i + 6];
        omegaBody[i]    = MRFData[i + 9];
        accBodyFrame[i] = MRFData[i + 12];
        alphaBody[i]    = MRFData[i + 15];
    }

    Array<OneD, NekDouble> omegaCrossVel(3, 0.0);
    Cross(omegaBody, velBody, omegaCrossVel);
    Array<OneD, NekDouble> accBodyInertial(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        accBodyInertial[i] = accBodyFrame[i] + omegaCrossVel[i];
    }

    const Array<OneD, NekDouble> velInertial =
        RotateBodyToInertial(q, velBody);
    const Array<OneD, NekDouble> accInertial =
        RotateBodyToInertial(q, accBodyInertial);
    const Array<OneD, NekDouble> omegaInertial =
        RotateBodyToInertial(q, omegaBody);
    const Array<OneD, NekDouble> alphaInertial =
        RotateBodyToInertial(q, alphaBody);

    m_outputStream << boost::format("%25.19e") % time << " ";
    for (int i = 0; i < m_spacedim; ++i)
    {
        m_outputStream << boost::format("%25.19e") % MRFData[i] << " "
                       << boost::format("%25.19e") % velInertial[i] << " "
                       << boost::format("%25.19e") % accInertial[i] << " ";
    }
    for (int i = 0; i < 3; ++i)
    {
        m_outputStream << boost::format("%25.19e") % MRFData[i + 3] << " "
                       << boost::format("%25.19e") % omegaInertial[i] << " "
                       << boost::format("%25.19e") % alphaInertial[i] << " ";
    }
    m_outputStream << std::endl;
}

void RigidSolver::UpdateMRFData(Array<OneD, NekDouble> &MRFData)
{
    /// MRFData:
    /// X, Y, Z, Theta_x, Theta_y, Theta_z, [inertial frame 0-5]
    /// U, V, W, Omega_x, Omega_y, Omega_z, [body frame 6-11]
    /// A_x, A_y, A_z, DOmega_x, DOmega_y, DOmega_z, [body frame 12-17]
    /// pivot_x, pivot_y, pivot_z, [body frame]
    /// Q0, Q1, Q2, Q3, [body-to-inertial quaternion]
    for (int i = 0; i < 18; ++i)
    {
        MRFData[i] = 0.0;
    }

    if (m_hasRotation)
    {
        MRFData[0] = m_inertialPosition[0];
        MRFData[1] = m_inertialPosition[1];
    }
    else
    {
        MRFData[0] = m_vel[0][0];
        MRFData[1] = m_vel[0][1];
    }

    if (m_hasRotation)
    {
        MRFData[5] = m_vel[0][m_spacedim];
    }

    for (int i = 0; i < m_spacedim; ++i)
    {
        MRFData[i + 6]  = m_vel[1][i];
        MRFData[i + 12] = m_vel[2][i];
    }

    if (m_hasRotation)
    {
        MRFData[11] = m_vel[1][m_spacedim];
        MRFData[17] = m_vel[2][m_spacedim];
    }
    SolverUtils::MovingFrame::StoreQuaternionInFrameData(
        SolverUtils::MovingFrame::QuaternionFromEulerZYX(
            MRFData[3], MRFData[4], MRFData[5]),
        MRFData);
}

void RigidSolver::UpdateFree3DMRFData(Array<OneD, NekDouble> &MRFData)
{
    using namespace SolverUtils::MovingFrame;

    Array<OneD, NekDouble> velBody(3, 0.0);
    Array<OneD, NekDouble> accBody(3, 0.0);
    Array<OneD, NekDouble> omegaBody(3, 0.0);
    Array<OneD, NekDouble> alphaBody(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        velBody[i]   = m_vel[1][i];
        accBody[i]   = m_vel[2][i];
        omegaBody[i] = m_vel[1][i + 3];
        alphaBody[i] = m_vel[2][i + 3];
    }
    for (const int axis : m_bodyAngularConstraints)
    {
        omegaBody[axis] = 0.0;
        alphaBody[axis] = 0.0;
    }
    Array<OneD, NekDouble> omegaCrossVel(3, 0.0);
    Cross(omegaBody, velBody, omegaCrossVel);
    Array<OneD, NekDouble> accInertialBody(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        accInertialBody[i] = accBody[i] + omegaCrossVel[i];
    }
    Array<OneD, NekDouble> velInertial =
        RotateBodyToInertial(m_quaternion, velBody);
    Array<OneD, NekDouble> accInertial =
        RotateBodyToInertial(m_quaternion, accInertialBody);
    if (!m_inertialTransConstraints.empty())
    {
        for (const int direction : m_inertialTransConstraints)
        {
            velInertial[direction] =
                m_inertialConstraintVelocity[direction];
            accInertial[direction] =
                m_inertialConstraintAcceleration[direction];
        }
        velBody = RotateInertialToBody(m_quaternion, velInertial);
        accInertialBody =
            RotateInertialToBody(m_quaternion, accInertial);
        Cross(omegaBody, velBody, omegaCrossVel);
        for (int i = 0; i < 3; ++i)
        {
            accBody[i] = accInertialBody[i] - omegaCrossVel[i];
        }
    }
    for (int i = 0; i < 3; ++i)
    {
        if (m_inertialTransConstraints.find(i) !=
            m_inertialTransConstraints.end())
        {
            if (m_hasInertialConstraintPosition[i])
            {
                m_inertialPosition[i] =
                    m_inertialConstraintPosition[i];
            }
            else
            {
                m_inertialPosition[i] +=
                    m_timestep * m_inertialVelocity[i] +
                    m_timestep * m_timestep *
                        ((0.5 - m_beta) * m_inertialAcceleration[i] +
                         m_beta * m_inertialConstraintAcceleration[i]);
            }
            m_inertialVelocity[i] = m_inertialConstraintVelocity[i];
            m_inertialAcceleration[i] =
                m_inertialConstraintAcceleration[i];
        }
        else
        {
            m_inertialPosition[i] +=
                m_timestep * m_inertialVelocity[i] +
                m_timestep * m_timestep *
                    ((0.5 - m_beta) * m_inertialAcceleration[i] +
                     m_beta * accInertial[i]);
            m_inertialVelocity[i]     = velInertial[i];
            m_inertialAcceleration[i] = accInertial[i];
        }
        MRFData[i]      = m_inertialPosition[i];
        MRFData[i + 6]  = velBody[i];
        MRFData[i + 9]  = omegaBody[i];
        MRFData[i + 12] = accBody[i];
        MRFData[i + 15] = alphaBody[i];
    }
    const Array<OneD, NekDouble> theta = EulerZYXFromQuaternion(m_quaternion);
    for (int i = 0; i < 3; ++i)
    {
        MRFData[i + 3] = theta[i];
        m_vel[0][i] = m_inertialPosition[i];
        m_vel[0][i + 3] = theta[i];
    }
    StoreQuaternionInFrameData(m_quaternion, MRFData);
}

void RigidSolver::WriteFree3DMRFOutput(
    const NekDouble &time, const Array<OneD, NekDouble> &MRFData)
{
    using namespace SolverUtils::MovingFrame;
    const Array<OneD, NekDouble> q = QuaternionFromFrameData(MRFData);
    Array<OneD, NekDouble> velocityBody(3, 0.0);
    Array<OneD, NekDouble> accelerationBody(3, 0.0);
    Array<OneD, NekDouble> omegaBody(3, 0.0);
    Array<OneD, NekDouble> alphaBody(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        velocityBody[i]     = MRFData[i + 6];
        accelerationBody[i] = MRFData[i + 12];
        omegaBody[i]        = MRFData[i + 9];
        alphaBody[i]        = MRFData[i + 15];
    }
    const Array<OneD, NekDouble> velocity =
        RotateBodyToInertial(q, velocityBody);
    Array<OneD, NekDouble> omegaCrossVelocity(3, 0.0);
    Cross(omegaBody, velocityBody, omegaCrossVelocity);
    Array<OneD, NekDouble> accelerationBodyInertial(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        accelerationBodyInertial[i] =
            accelerationBody[i] + omegaCrossVelocity[i];
    }
    const Array<OneD, NekDouble> acceleration =
        RotateBodyToInertial(q, accelerationBodyInertial);
    const Array<OneD, NekDouble> omega = RotateBodyToInertial(q, omegaBody);
    const Array<OneD, NekDouble> alpha = RotateBodyToInertial(q, alphaBody);
    m_outputStream << boost::format("%25.19e") % time << " ";
    for (int i = 0; i < 3; ++i)
    {
        m_outputStream << boost::format("%25.19e") % MRFData[i] << " "
                       << boost::format("%25.19e") % velocity[i] << " "
                       << boost::format("%25.19e") % acceleration[i] << " ";
    }
    for (int i = 0; i < 3; ++i)
    {
        m_outputStream << boost::format("%25.19e") % MRFData[i + 3] << " "
                       << boost::format("%25.19e") % omega[i] << " "
                       << boost::format("%25.19e") % alpha[i] << " ";
    }
    m_outputStream << std::endl;
}

void RigidSolver::SolveBodyMotion(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                                  const Array<OneD, NekDouble> &forcebody,
                                  std::map<int, NekDouble> &Dirs)
{
    if (0 == m_solveType)
    {
        SolveInertialFrame(bodyVel, forcebody, Dirs);
    }
    else if (1 == m_solveType)
    {
        SolveRotOneFree(bodyVel, forcebody, Dirs);
    }
    else if (2 == m_solveType)
    {
        SolveBodyFrame(bodyVel, forcebody, Dirs);
    }
    else if (eFree3D6DoF == m_solveType)
    {
        SolveFree3D6DoF(bodyVel, forcebody, Dirs);
    }
    else
    {
        ASSERTL0(false, "Unsupported rigid solver type.");
    }
}

void RigidSolver::SolveFree3D6DoF(
    Array<OneD, Array<OneD, NekDouble>> &bodyVel,
    const Array<OneD, NekDouble> &forcebody, std::map<int, NekDouble> &Dirs)
{
    using namespace SolverUtils::MovingFrame;
    ASSERTL0(Dirs.empty(),
             "Full free 3D 6DoF motion cannot prescribe individual DoFs.");
    ASSERTL0(forcebody.size() >= 6,
             "Full free 3D 6DoF motion requires three forces and moments.");

    Array<OneD, NekDouble> force(6, 0.0);
    Array<OneD, NekDouble> externalForce(3, 0.0);
    Array<OneD, NekDouble> externalMoment(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        externalForce[i]  = m_extForceXYZ[i];
        externalMoment[i] = m_extForceXYZ[i + 3];
    }
    externalForce = RotateInertialToBody(m_quaternion, externalForce);
    externalMoment = RotateInertialToBody(m_quaternion, externalMoment);
    for (int i = 0; i < 3; ++i)
    {
        force[i]     = forcebody[i] + externalForce[i];
        force[i + 3] = forcebody[i + 3] + externalMoment[i];
    }

    Array<OneD, Array<OneD, NekDouble>> oldState(3);
    for (int i = 0; i < 3; ++i)
    {
        oldState[i] = Array<OneD, NekDouble>(bodyVel[i].size(), 0.0);
        Vmath::Vcopy(bodyVel[i].size(), bodyVel[i], 1, oldState[i], 1);
    }
    Array<OneD, NekDouble> trialVelocity(6, 0.0);
    Array<OneD, NekDouble> omegaOld(3, 0.0);
    for (int i = 0; i < 6; ++i)
    {
        trialVelocity[i] = oldState[1][i];
        if (i < 3)
        {
            omegaOld[i] = oldState[1][i + 3];
        }
    }

    const auto addCrossMatrix = [](const Array<OneD, NekDouble> &a,
                                   const NekDouble scale,
                                   Array<OneD, NekDouble> &matrix,
                                   const int rowOffset, const int colOffset) {
        constexpr int nDofs = 6;
        matrix[(rowOffset + 0) + (colOffset + 1) * nDofs] +=
            -scale * a[2];
        matrix[(rowOffset + 0) + (colOffset + 2) * nDofs] +=
            scale * a[1];
        matrix[(rowOffset + 1) + (colOffset + 0) * nDofs] +=
            scale * a[2];
        matrix[(rowOffset + 1) + (colOffset + 2) * nDofs] +=
            -scale * a[0];
        matrix[(rowOffset + 2) + (colOffset + 0) * nDofs] +=
            -scale * a[1];
        matrix[(rowOffset + 2) + (colOffset + 1) * nDofs] +=
            scale * a[0];
    };

    bool converged = false;
    for (int iter = 0; iter < m_nonlinearMaxIterations; ++iter)
    {
        Array<OneD, NekDouble> velocity(3, 0.0);
        Array<OneD, NekDouble> omega(3, 0.0);
        Array<OneD, NekDouble> angularMomentum(3, 0.0);
        Array<OneD, NekDouble> nonlinear(6, 0.0);
        Array<OneD, NekDouble> jacobian(36, 0.0);
        for (int i = 0; i < 3; ++i)
        {
            velocity[i]        = trialVelocity[i];
            omega[i]           = trialVelocity[i + 3];
            angularMomentum[i] = m_rotationInertia[i] * omega[i];
        }

        Array<OneD, NekDouble> omegaCrossVelocity(3, 0.0);
        Array<OneD, NekDouble> gyroMoment(3, 0.0);
        Cross(omega, velocity, omegaCrossVelocity);
        Cross(omega, angularMomentum, gyroMoment);
        for (int i = 0; i < 3; ++i)
        {
            nonlinear[i]     = m_mass * omegaCrossVelocity[i];
            nonlinear[i + 3] = gyroMoment[i];
        }

        // Linearise m Omega x u and Omega x (I Omega) about the
        // current Picard/Newton iterate in the body frame.
        addCrossMatrix(omega, m_mass, jacobian, 0, 0);
        addCrossMatrix(velocity, -m_mass, jacobian, 0, 3);
        addCrossMatrix(angularMomentum, -1.0, jacobian, 3, 3);
        // [Omega]_x I: each column of the skew matrix is scaled by the
        // corresponding principal moment of inertia.
        jacobian[3 + 4 * 6] += -omega[2] * m_rotationInertia[1];
        jacobian[3 + 5 * 6] += omega[1] * m_rotationInertia[2];
        jacobian[4 + 3 * 6] += omega[2] * m_rotationInertia[0];
        jacobian[4 + 5 * 6] += -omega[0] * m_rotationInertia[2];
        jacobian[5 + 3 * 6] += -omega[1] * m_rotationInertia[0];
        jacobian[5 + 4 * 6] += omega[0] * m_rotationInertia[1];

        for (int i = 0; i < 3; ++i)
        {
            Vmath::Vcopy(oldState[i].size(), oldState[i], 1, bodyVel[i], 1);
        }
        if (m_inertialTransConstraints.empty() &&
            m_bodyAngularConstraints.empty())
        {
            m_bodySolver.SolveFreeVarMat6DoF(bodyVel, force, nonlinear,
                                              jacobian, trialVelocity);
        }
        else
        {
            const int nConstraints = m_inertialTransConstraints.size() +
                                     m_bodyAngularConstraints.size();
            Array<OneD, NekDouble> constraints(6 * nConstraints, 0.0);
            Array<OneD, NekDouble> values(nConstraints, 0.0);
            Array<OneD, NekDouble> omegaMid(3, 0.0);
            for (int i = 0; i < 3; ++i)
            {
                omegaMid[i] = 0.5 * (omegaOld[i] + omega[i]);
            }
            const Array<OneD, NekDouble> constraintQuaternion =
                NormalizeQuaternion(IntegrateQuaternionBodyOmega(
                    m_quaternion, omegaMid, m_timestep));
            int k = 0;
            for (const int direction : m_inertialTransConstraints)
            {
                Array<OneD, NekDouble> directionInertial(3, 0.0);
                directionInertial[direction] = 1.0;
                const Array<OneD, NekDouble> directionBody =
                    RotateInertialToBody(constraintQuaternion,
                                         directionInertial);
                for (int j = 0; j < 3; ++j)
                {
                    constraints[k * 6 + j] = directionBody[j];
                }
                values[k] = m_inertialConstraintVelocity[direction];
                ++k;
            }
            for (const int axis : m_bodyAngularConstraints)
            {
                constraints[k * 6 + 3 + axis] = 1.0;
                values[k] = 0.0;
                ++k;
            }
            m_bodySolver.SolveFreeVarMatNDofConstrained(
                bodyVel, force, nonlinear, jacobian, trialVelocity,
                constraints, values, 6, nConstraints);
        }

        NekDouble deltaNorm = 0.0;
        NekDouble stateNorm = 0.0;
        for (int i = 0; i < 6; ++i)
        {
            const NekDouble delta = bodyVel[1][i] - trialVelocity[i];
            deltaNorm += delta * delta;
            stateNorm += bodyVel[1][i] * bodyVel[1][i];
            trialVelocity[i] = bodyVel[1][i];
        }
        if (std::sqrt(deltaNorm) <=
            m_nonlinearTolerance * std::max(1.0, std::sqrt(stateNorm)))
        {
            converged = true;
            break;
        }
    }
    ASSERTL0(converged,
             "The nonlinear full free 3D 6DoF rigid-body solve did not "
             "converge; reduce TimeStep or increase "
             "RigidBodyNonlinearMaxIterations.");

    Array<OneD, NekDouble> omegaNew(3, 0.0);
    Array<OneD, NekDouble> omegaMid(3, 0.0);
    for (int i = 0; i < 3; ++i)
    {
        omegaNew[i] = bodyVel[1][i + 3];
        omegaMid[i] = 0.5 * (omegaOld[i] + omegaNew[i]);
    }
    const Array<OneD, NekDouble> qNew = NormalizeQuaternion(
        IntegrateQuaternionBodyOmega(m_quaternion, omegaMid, m_timestep));
    m_quaternion = MakeQuaternionContinuous(m_quaternion, qNew);
}

void RigidSolver::SolveInertialFrame(
    Array<OneD, Array<OneD, NekDouble>> &bodyVel,
    const Array<OneD, NekDouble> &forcebody, std::map<int, NekDouble> &Dirs)
{
    // only translational motion or prescribed tranlation with rotational
    // solve in absolute frame and transform to body frame
    m_bodySolver.SolvePrescribed(bodyVel, Dirs);
    if (m_hasFreeMotion)
    {
        Array<OneD, NekDouble> force(6, 0.);
        for (int i = 0; i < m_spacedim; ++i)
        {
            force[i] = forcebody[i] + m_extForceXYZ[i];
        }
        force[m_spacedim] = forcebody[5] + m_extForceXYZ[5];
        m_bodySolver.SolveFreeFixMat(bodyVel, force);
    }
}

// with rotational and one free tranlation
void RigidSolver::SolveRotOneFree(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                                  const Array<OneD, NekDouble> &forcebody,
                                  std::map<int, NekDouble> &Dirs)
{
    // one direction free
    NekDouble uy = 0.;
    if (Dirs.find(1) != Dirs.end())
    {
        uy = Dirs[1];
        Dirs.erase(1);
    }
    m_bodySolver.SolvePrescribed(bodyVel, Dirs);
    Array<OneD, NekDouble> force(6, 0.), angle(3, 0.);
    angle[2] = bodyVel[0][m_spacedim];
    m_frame.SetAngle(angle);
    m_frame.IneritalToBody(3, m_extForceXYZ, force);
    for (int i = 0; i < m_spacedim; ++i)
    {
        force[i] = forcebody[i] + force[i];
    }
    force[0] += m_mass * bodyVel[1][m_spacedim] * bodyVel[1][m_spacedim] *
                m_pivotdistance;
    force[1] -= m_mass * bodyVel[2][m_spacedim] * m_pivotdistance;
    m_bodySolver.SolveOneFree(bodyVel, force, angle, uy, m_mass);
}

// with rotational and all free tranlation
void RigidSolver::SolveFreeRigidBody2D(
    Array<OneD, Array<OneD, NekDouble>> &bodyVel,
    const Array<OneD, NekDouble> &forcebody, std::map<int, NekDouble> &)
{
    Array<OneD, Array<OneD, NekDouble>> oldState(3);
    for (int i = 0; i < 3; ++i)
    {
        oldState[i] = Array<OneD, NekDouble>(3, 0.0);
        Vmath::Vcopy(3, bodyVel[i], 1, oldState[i], 1);
    }
    Array<OneD, NekDouble> trial(3, 0.0);
    Vmath::Vcopy(3, oldState[1], 1, trial, 1);
    bool converged = false;
    for (int iter = 0; iter < m_nonlinearMaxIterations; ++iter)
    {
        Array<OneD, NekDouble> force(3, 0.0), angle(3, 0.0);
        angle[2] = bodyVel[0][2];
        m_frame.SetAngle(angle);
        m_frame.IneritalToBody(3, m_extForceXYZ, force);
        force[0] += forcebody[0];
        force[1] += forcebody[1];
        force[2] = forcebody[5] + m_extForceXYZ[5];

        Array<OneD, NekDouble> nonlinear(3, 0.0), jacobian(9, 0.0);
        nonlinear[0] = -m_mass * trial[2] * trial[1];
        nonlinear[1] = m_mass * trial[2] * trial[0];
        jacobian[0 + 1 * 3] = -m_mass * trial[2];
        jacobian[0 + 2 * 3] = -m_mass * trial[1];
        jacobian[1 + 0 * 3] = m_mass * trial[2];
        jacobian[1 + 2 * 3] = m_mass * trial[0];

        for (int i = 0; i < 3; ++i)
        {
            Vmath::Vcopy(3, oldState[i], 1, bodyVel[i], 1);
        }
        if (m_inertialTransConstraints.empty())
        {
            m_bodySolver.SolveFreeVarMatNDof(
                bodyVel, force, nonlinear, jacobian, trial,
                m_bodySolver.m_motionDofs);
        }
        else
        {
            const int nConstraints = m_inertialTransConstraints.size();
            Array<OneD, NekDouble> constraints(3 * nConstraints, 0.0);
            Array<OneD, NekDouble> values(nConstraints, 0.0);
            int k = 0;
            const NekDouble c = cos(angle[2]);
            const NekDouble s = sin(angle[2]);
            for (const int direction : m_inertialTransConstraints)
            {
                if (direction == 0)
                {
                    constraints[k * 3]     = c;
                    constraints[k * 3 + 1] = -s;
                }
                else
                {
                    constraints[k * 3]     = s;
                    constraints[k * 3 + 1] = c;
                }
                values[k] = m_inertialConstraintVelocity[direction];
                ++k;
            }
            m_bodySolver.SolveFreeVarMatNDofConstrained(
                bodyVel, force, nonlinear, jacobian, trial, constraints, values,
                m_bodySolver.m_motionDofs, nConstraints);
        }
        NekDouble delta2 = 0.0, state2 = 0.0;
        for (int i = 0; i < m_bodySolver.m_motionDofs; ++i)
        {
            const int i1 = m_bodySolver.m_index[i];
            const NekDouble delta = bodyVel[1][i1] - trial[i1];
            delta2 += delta * delta;
            state2 += bodyVel[1][i1] * bodyVel[1][i1];
            trial[i1] = bodyVel[1][i1];
        }
        if (std::sqrt(delta2) <= m_nonlinearTolerance *
                                     std::max(1.0, std::sqrt(state2)))
        {
            converged = true;
            break;
        }
    }
    ASSERTL0(converged, "The unified 2D rigid-body Newton solve did not converge.");
}

void RigidSolver::SolveBodyFrame(Array<OneD, Array<OneD, NekDouble>> &bodyVel,
                                 const Array<OneD, NekDouble> &forcebody,
                                 std::map<int, NekDouble> &Dirs)
{
    m_bodySolver.SolvePrescribed(bodyVel, Dirs); // at most 1 rotation
    Array<OneD, NekDouble> force(6, 0.), angle(3, 0.);
    if (!m_hasFreeMotion)
    {
        angle[2] = bodyVel[0][m_spacedim];
        m_frame.SetAngle(angle);
        m_frame.IneritalToBody(m_spacedim, bodyVel[1], bodyVel[1]);
        m_frame.IneritalToBody(m_spacedim, bodyVel[2], bodyVel[2]);
        bodyVel[2][0] += bodyVel[1][m_spacedim] * bodyVel[1][1];
        bodyVel[2][1] -= bodyVel[1][m_spacedim] * bodyVel[1][0];
        return;
    }
    if (m_dirDoFs.find(m_spacedim) != m_dirDoFs.end())
    {
        // known rotation
        angle[2] = bodyVel[0][m_spacedim];
        m_frame.SetAngle(angle);
        m_frame.IneritalToBody(m_spacedim, m_extForceXYZ, force);
        for (int i = 0; i < m_spacedim; ++i)
        {
            force[i] = forcebody[i] + force[i];
        }
        force[0] += m_mass * bodyVel[1][m_spacedim] * bodyVel[1][m_spacedim] *
                    m_pivotdistance;
        force[1] -= m_mass * bodyVel[2][m_spacedim] * m_pivotdistance;
        m_bodySolver.SolveFreeVarMat(bodyVel, force, m_mass);
    }
    else
    {
        if (m_useUnifiedFreeRigidBody &&
            fabs(m_pivotdistance) < NekConstants::kNekZeroTol)
        {
            SolveFreeRigidBody2D(bodyVel, forcebody, Dirs);
            return;
        }
        // all free
        Array<OneD, Array<OneD, NekDouble>> tmpbodyVel(bodyVel.size());
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            tmpbodyVel[i] = Array<OneD, NekDouble>(bodyVel[i].size());
            Vmath::Vcopy(bodyVel[i].size(), bodyVel[i], 1, tmpbodyVel[i], 1);
        }
        for (int iter = 0; iter < 2; ++iter)
        {
            if (iter > 0)
            {
                for (size_t i = 0; i < bodyVel.size(); ++i)
                {
                    Vmath::Vcopy(bodyVel[i].size() - 1, bodyVel[i], 1,
                                 tmpbodyVel[i], 1);
                }
            }
            angle[2] = tmpbodyVel[0][m_spacedim];
            m_frame.SetAngle(angle);
            m_frame.IneritalToBody(3, m_extForceXYZ, force);
            for (int i = 0; i < m_spacedim; ++i)
            {
                force[i] = forcebody[i] + force[i];
            }
            force[m_spacedim] = forcebody[5] + m_extForceXYZ[5];
            m_bodySolver.SolveFreeVarMat(tmpbodyVel, force, m_mass);
        }
        // copy final results
        for (size_t i = 0; i < bodyVel.size(); ++i)
        {
            Vmath::Vcopy(bodyVel[i].size(), tmpbodyVel[i], 1, bodyVel[i], 1);
        }
    }
}

void Newmark_BetaSolver::SolveFreeVarMat(Array<OneD, Array<OneD, NekDouble>> u,
                                         Array<OneD, NekDouble> force,
                                         const NekDouble mass)
{
    Array<OneD, NekDouble> bm(m_motionDofs, 0.);
    Array<OneD, NekDouble> bk(m_motionDofs, 0.);
    double *dMatrix = new double[m_motionDofs * m_motionDofs];
    int *ipiv       = new int[m_motionDofs];
    double *drhs    = new double[m_motionDofs];
    int info;

    for (int j = 0; j < m_motionDofs; ++j)
    {
        int j1 = m_index[j];
        bm[j]  = m_coeffs[0] * u[1][j1] + m_coeffs[1] * u[2][j1];
        bk[j]  = u[0][j1] + m_coeffs[3] * u[1][j1] + m_coeffs[4] * u[2][j1];
    }
    Array<OneD, NekDouble> rhs(m_motionDofs, 0.);
    for (int i = 0; i < m_motionDofs; ++i)
    {
        rhs[i] = force[m_index[i]];
        for (int j = 0; j < m_motionDofs; ++j)
        {
            rhs[i] += m_M[i][j] * bm[j] - m_K[i][j] * bk[j];
        }
        for (int j = m_motionDofs; j < m_rows; ++j)
        {
            int j1 = m_index[j];
            rhs[i] -= m_M[i][j] * u[2][j1] + m_C[i][j] * u[1][j1] +
                      m_K[i][j] * u[0][j1];
        }
    }
    for (int i = 0; i < m_motionDofs; ++i)
    {
        drhs[i] = rhs[i];
        for (int j = 0; j < m_motionDofs; ++j)
        {
            dMatrix[j * m_motionDofs + i] = m_Matrix[i][j];
        }
    }
    dMatrix[1] += mass * u[1][m_rows - 1];
    dMatrix[m_motionDofs] += -mass * u[1][m_rows - 1];

    Lapack::DoSgetrf(m_motionDofs, m_motionDofs, dMatrix, m_motionDofs, ipiv,
                     info);
    Lapack::Dgetrs('N', m_motionDofs, 1, dMatrix, m_motionDofs, ipiv, drhs,
                   m_motionDofs, info);
    for (int j = 0; j < m_motionDofs; ++j)
    {
        int j1   = m_index[j];
        u[1][j1] = drhs[j];
    }

    for (int j = 0; j < m_motionDofs; ++j)
    {
        int j1   = m_index[j];
        u[0][j1] = m_coeffs[2] * u[1][j1] + bk[j];
        u[2][j1] = m_coeffs[0] * u[1][j1] - bm[j];
    }

    delete[] dMatrix;
    delete[] drhs;
    delete[] ipiv;
}

void RigidSolver::SetNewmarkBetaSolver(Array<OneD, NekDouble> &AddedMass)
{
    int NumDof = m_free3D6DoF ? 6 : m_spacedim + 1;
    if (AddedMass.size() >= NumDof * NumDof)
    {
        Vmath::Vadd(NumDof * NumDof, AddedMass, 1, m_M, 1, m_M, 1);
    }
    m_bodySolver.SetNewmarkBeta(m_beta, m_gamma, m_timestep, m_M, m_C, m_K,
                                m_dirDoFs, m_solveType);
}

void RigidSolver::SetInitialConditions(
    const LibUtilities::SessionReaderSharedPtr session,
    Array<OneD, NekDouble> MRFData)
{
    NekDouble time = 0.;
    std::map<std::string, std::string> fieldMetaDataMap;
    std::vector<std::string> strFrameData = {
        "X",   "Y",   "Z",   "Theta_x",  "Theta_y",  "Theta_z",
        "U",   "V",   "W",   "Omega_x",  "Omega_y",  "Omega_z",
        "A_x", "A_y", "A_z", "DOmega_x", "DOmega_y", "DOmega_z",
        "Q0",  "Q1",  "Q2",  "Q3"};
    std::map<std::string, NekDouble> fileData;
    const int nLegacyFrameData = 18;
    const auto hasLegacyFrameData =
        [&fileData, &strFrameData, nLegacyFrameData]() {
            for (int i = 0; i < nLegacyFrameData; ++i)
            {
                if (fileData.find(strFrameData[i]) == fileData.end())
                {
                    return false;
                }
            }
            return true;
        };
    const auto hasQuaternionData = [&fileData]() {
        return fileData.find("Q0") != fileData.end() &&
               fileData.find("Q1") != fileData.end() &&
               fileData.find("Q2") != fileData.end() &&
               fileData.find("Q3") != fileData.end();
    };
    if (session->DefinesFunction("InitialConditions"))
    {
        for (int i = 0; i < session->GetVariables().size(); ++i)
        {
            if (session->GetFunctionType("InitialConditions",
                                         session->GetVariable(i)) ==
                LibUtilities::eFunctionTypeFile)
            {
                std::string filename = session->GetFunctionFilename(
                    "InitialConditions", session->GetVariable(i));
                fs::path pfilename(filename);
                // redefine path for parallel file which is in directory
                if (fs::is_directory(pfilename))
                {
                    fs::path metafile("Info.xml");
                    fs::path fullpath = pfilename / metafile;
                    filename          = LibUtilities::PortablePath(fullpath);
                }
                LibUtilities::FieldIOSharedPtr fld =
                    LibUtilities::FieldIO::CreateForFile(session, filename);
                fld->ImportFieldMetaData(filename, fieldMetaDataMap);

                // check to see if time is defined
                if (fieldMetaDataMap != LibUtilities::NullFieldMetaDataMap)
                {
                    if (fieldMetaDataMap.find("Time") != fieldMetaDataMap.end())
                    {
                        time = std::stod(fieldMetaDataMap["Time"]);
                    }
                    fileData.clear();
                    for (auto &var : strFrameData)
                    {
                        if (fieldMetaDataMap.find(var) !=
                            fieldMetaDataMap.end())
                        {
                            fileData[var] = std::stod(fieldMetaDataMap[var]);
                        }
                    }
                    if (hasLegacyFrameData())
                    {
                        break;
                    }
                }
            }
        }
    }
    if (session->DefinesCmdLineArgument("set-start-time"))
    {
        time = std::stod(
            session->GetCmdLineArgument<std::string>("set-start-time"));
    }
    if (hasLegacyFrameData())
    {
        for (int i = 0; i < nLegacyFrameData; ++i)
        {
            MRFData[i] = fileData[strFrameData[i]];
        }
        if (hasQuaternionData())
        {
            Array<OneD, NekDouble> q(4, 0.0);
            for (int i = 0; i < 4; ++i)
            {
                q[i] = fileData[strFrameData[nLegacyFrameData + i]];
            }
            SolverUtils::MovingFrame::StoreQuaternionInFrameData(q, MRFData);
        }
        else
        {
            SolverUtils::MovingFrame::StoreQuaternionInFrameData(
                SolverUtils::MovingFrame::QuaternionFromEulerZYX(
                    MRFData[3], MRFData[4], MRFData[5]),
                MRFData);
        }
    }
    if (m_prescribed3DMRF)
    {
        UpdatePrescribedMRFData(time, MRFData);
        if (m_isRoot)
        {
            WritePrescribedMRFOutput(time, MRFData);
        }
        m_currentTime = time;
        return;
    }
    if (hasLegacyFrameData())
    {
        if (m_free3D6DoF)
        {
            for (int i = 0; i < 3; ++i)
            {
                m_inertialPosition[i] = fileData[strFrameData[i]];
                m_vel[0][i] = m_inertialPosition[i];
                m_vel[1][i] = fileData[strFrameData[i + 6]];
                m_vel[2][i] = fileData[strFrameData[i + 12]];
                m_vel[1][i + 3] = fileData[strFrameData[i + 9]];
                m_vel[2][i + 3] = fileData[strFrameData[i + 15]];
            }
            m_quaternion = hasQuaternionData()
                ? SolverUtils::MovingFrame::QuaternionFromFrameData(MRFData)
                : SolverUtils::MovingFrame::QuaternionFromEulerZYX(
                      fileData[strFrameData[3]], fileData[strFrameData[4]],
                      fileData[strFrameData[5]]);
            const Array<OneD, NekDouble> theta =
                SolverUtils::MovingFrame::EulerZYXFromQuaternion(m_quaternion);
            for (int i = 0; i < 3; ++i)
            {
                m_vel[0][i + 3] = theta[i];
            }
        }
        else
        {
        int NumDofm1 = m_vel[0].size() - 1;
        for (int i = 0; i < m_spacedim; ++i)
        {
            m_vel[0][i] = fileData[strFrameData[i]];
            m_vel[1][i] = fileData[strFrameData[i + 6]];
            m_vel[2][i] = fileData[strFrameData[i + 12]];
        }
        m_vel[0][NumDofm1] = fileData[strFrameData[5]];
        m_vel[1][NumDofm1] = fileData[strFrameData[11]];
        m_vel[2][NumDofm1] = fileData[strFrameData[17]];
        if (m_hasRotation)
        {
            m_inertialPosition[0] = m_vel[0][0];
            m_inertialPosition[1] = m_vel[0][1];
        }
        }
    }
    std::map<int, NekDouble> Dirs;
    UpdatePrescribed(time, Dirs);
    SetInitialConditions(Dirs);
    if (m_spacedim == 2 && m_hasRotation &&
        !m_inertialTransConstraints.empty())
    {
        for (const int direction : m_inertialTransConstraints)
        {
            if (m_hasInertialConstraintPosition[direction])
            {
                m_inertialPosition[direction] =
                    m_inertialConstraintPosition[direction];
            }
        }
        Array<OneD, NekDouble> angle(3, 0.0);
        angle[2] = m_vel[0][m_spacedim];
        m_frame.SetAngle(angle);

        Array<OneD, NekDouble> velocityBody(3, 0.0);
        Array<OneD, NekDouble> velocityInertial(3, 0.0);
        Array<OneD, NekDouble> accelerationBody(3, 0.0);
        Array<OneD, NekDouble> accelerationInertial(3, 0.0);
        for (int i = 0; i < 2; ++i)
        {
            velocityBody[i]     = m_vel[1][i];
            accelerationBody[i] = m_vel[2][i];
        }
        m_frame.BodyToInerital(2, velocityBody, velocityInertial);
        for (const int direction : m_inertialTransConstraints)
        {
            velocityInertial[direction] =
                m_inertialConstraintVelocity[direction];
        }
        m_frame.IneritalToBody(2, velocityInertial, velocityBody);
        for (int i = 0; i < 2; ++i)
        {
            m_vel[1][i] = velocityBody[i];
        }

        const NekDouble omega = m_vel[1][m_spacedim];
        accelerationBody[0] -= omega * velocityBody[1];
        accelerationBody[1] += omega * velocityBody[0];
        m_frame.BodyToInerital(2, accelerationBody, accelerationInertial);
        for (const int direction : m_inertialTransConstraints)
        {
            accelerationInertial[direction] =
                m_inertialConstraintAcceleration[direction];
        }
        m_frame.IneritalToBody(2, accelerationInertial, accelerationBody);
        accelerationBody[0] += omega * velocityBody[1];
        accelerationBody[1] -= omega * velocityBody[0];
        for (int i = 0; i < 2; ++i)
        {
            m_vel[2][i] = accelerationBody[i];
        }
    }
    if (m_spacedim == 2 && m_hasRotation)
    {
        Array<OneD, NekDouble> angle(3, 0.0);
        Array<OneD, NekDouble> velocityBody(3, 0.0);
        Array<OneD, NekDouble> accelerationBody(3, 0.0);
        angle[2] = m_vel[0][m_spacedim];
        for (int i = 0; i < 2; ++i)
        {
            velocityBody[i]     = m_vel[1][i];
            accelerationBody[i] = m_vel[2][i];
        }
        m_frame.SetAngle(angle);
        m_frame.BodyToInerital(2, velocityBody, m_inertialVelocity);
        m_frame.BodyToInerital(2, accelerationBody,
                                m_inertialAcceleration);
        const NekDouble omega = m_vel[1][m_spacedim];
        m_inertialAcceleration[0] -= omega * m_inertialVelocity[1];
        m_inertialAcceleration[1] += omega * m_inertialVelocity[0];
        for (const int direction : m_inertialTransConstraints)
        {
            m_inertialVelocity[direction] =
                m_inertialConstraintVelocity[direction];
            m_inertialAcceleration[direction] =
                m_inertialConstraintAcceleration[direction];
        }
    }
    if (m_free3D6DoF)
    {
        Array<OneD, NekDouble> velocityBody(3, 0.0);
        Array<OneD, NekDouble> accelerationBody(3, 0.0);
        Array<OneD, NekDouble> omegaBody(3, 0.0);
        for (int i = 0; i < 3; ++i)
        {
            velocityBody[i]     = m_vel[1][i];
            accelerationBody[i] = m_vel[2][i];
            omegaBody[i]        = m_vel[1][i + 3];
        }
        Array<OneD, NekDouble> omegaCrossVelocity(3, 0.0);
        SolverUtils::MovingFrame::Cross(omegaBody, velocityBody,
                                         omegaCrossVelocity);
        for (int i = 0; i < 3; ++i)
        {
            accelerationBody[i] += omegaCrossVelocity[i];
        }
        m_inertialVelocity = SolverUtils::MovingFrame::RotateBodyToInertial(
            m_quaternion, velocityBody);
        m_inertialAcceleration =
            SolverUtils::MovingFrame::RotateBodyToInertial(
                m_quaternion, accelerationBody);
        for (const int direction : m_inertialTransConstraints)
        {
            if (m_hasInertialConstraintPosition[direction])
            {
                m_inertialPosition[direction] =
                    m_inertialConstraintPosition[direction];
            }
            m_inertialVelocity[direction] =
                m_inertialConstraintVelocity[direction];
            m_inertialAcceleration[direction] =
                m_inertialConstraintAcceleration[direction];
        }
        velocityBody = SolverUtils::MovingFrame::RotateInertialToBody(
            m_quaternion, m_inertialVelocity);
        accelerationBody = SolverUtils::MovingFrame::RotateInertialToBody(
            m_quaternion, m_inertialAcceleration);
        SolverUtils::MovingFrame::Cross(omegaBody, velocityBody,
                                         omegaCrossVelocity);
        for (int i = 0; i < 3; ++i)
        {
            m_vel[1][i] = velocityBody[i];
            m_vel[2][i] = accelerationBody[i] - omegaCrossVelocity[i];
        }
        const Array<OneD, NekDouble> theta =
            SolverUtils::MovingFrame::EulerZYXFromQuaternion(m_quaternion);
        for (int i = 0; i < 3; ++i)
        {
            MRFData[i]      = m_inertialPosition[i];
            MRFData[i + 3]  = theta[i];
            MRFData[i + 6]  = m_vel[1][i];
            MRFData[i + 9]  = m_vel[1][i + 3];
            MRFData[i + 12] = m_vel[2][i];
            MRFData[i + 15] = m_vel[2][i + 3];
        }
        SolverUtils::MovingFrame::StoreQuaternionInFrameData(m_quaternion,
                                                               MRFData);
        if (m_isRoot)
        {
            WriteFree3DMRFOutput(time, MRFData);
        }
        return;
    }
    UpdateMRFData(MRFData);
    // output initial status for rigid body
    if (m_isRoot)
    {
        Array<OneD, Array<OneD, NekDouble>> tmpVel;
        if (m_hasRotation)
        {
            tmpVel = Array<OneD, Array<OneD, NekDouble>>(m_vel.size());
            for (size_t i = 0; i < m_vel.size(); ++i)
            {
                tmpVel[i] = Array<OneD, NekDouble>(m_vel[i].size());
                Vmath::Vcopy(m_vel[i].size(), m_vel[i], 1, tmpVel[i], 1);
            }
            // transform u, du to body frame
            Array<OneD, NekDouble> angle(3, 0.);
            angle[2] = tmpVel[0][m_spacedim];
            m_frame.SetAngle(angle);
            m_frame.BodyToInerital(m_spacedim, tmpVel[1], tmpVel[1]);
            m_frame.BodyToInerital(m_spacedim, tmpVel[2], tmpVel[2]);
            tmpVel[2][0] -= tmpVel[1][m_spacedim] * tmpVel[1][1];
            tmpVel[2][1] += tmpVel[1][m_spacedim] * tmpVel[1][0];
            tmpVel[0][0] = m_inertialPosition[0];
            tmpVel[0][1] = m_inertialPosition[1];
        }
        else
        {
            tmpVel = m_vel;
        }
        m_outputStream << boost::format("%25.19e") % time << " ";
        for (size_t i = 0; i < tmpVel[0].size(); ++i)
        {
            m_outputStream << boost::format("%25.19e") % tmpVel[0][i] << " "
                           << boost::format("%25.19e") % tmpVel[1][i] << " "
                           << boost::format("%25.19e") % tmpVel[2][i] << " ";
        }
        m_outputStream << std::endl;
    }
}

void RigidSolver::GetFilterInfo(
    const LibUtilities::SessionReaderSharedPtr session,
    std::map<std::string, std::string> &vParams)
{
    TiXmlElement *pForce      = session->GetElement("Nektar/RIGIDSOLVER");
    const TiXmlElement *param = pForce->FirstChildElement("BOUNDARY");
    ASSERTL0(param, "Body surface should be assigned");

    vParams["Boundary"]           = param->GetText();
    const TiXmlElement *pivotElmt = pForce->FirstChildElement("PIVOTPOINT");
    if (pivotElmt)
    {
        std::string pstr = pivotElmt->GetText();
        std::replace(pstr.begin(), pstr.end(), ',', ' ');
        vParams["MomentPoint"] = pstr;
    }
}

void RigidSolver::SetInitialConditions(std::map<int, NekDouble> &Dirs)
{
    for (auto it : Dirs)
    {
        int NumDof  = m_vel[0].size();
        int NumDof2 = NumDof << 1;
        if (it.first < m_vel[0].size())
        {
            m_vel[1][it.first] = it.second;
        }
        else if (it.first < NumDof2)
        {
            m_vel[0][it.first - NumDof] = it.second;
        }
        else
        {
            m_vel[2][it.first - NumDof2] = it.second;
        }
    }
}

void Newmark_BetaSolver::SetNewmarkBeta(NekDouble beta, NekDouble gamma,
                                        NekDouble dt, Array<OneD, NekDouble> M,
                                        Array<OneD, NekDouble> C,
                                        Array<OneD, NekDouble> K,
                                        std::set<int> DirDoFs, int solveType)
{
    m_coeffs    = Array<OneD, NekDouble>(5, 0.);
    m_coeffs[0] = 1. / (gamma * dt);
    m_coeffs[1] = 1. / gamma - 1.;
    m_coeffs[2] = beta * dt / gamma;
    m_coeffs[3] = dt * (1. - beta / gamma);
    m_coeffs[4] = (0.5 - beta / gamma) * dt * dt;

    m_rows = sqrt(M.size());
    m_index.resize(m_rows, -1);
    m_motionDofs = 0;
    for (int i = 0; i < m_rows; ++i)
    {
        if (DirDoFs.find(i) == DirDoFs.end())
        {
            m_index[m_motionDofs++] = i;
        }
    }
    for (int i = 0, count = m_motionDofs; i < m_rows; ++i)
    {
        if (DirDoFs.find(i) != DirDoFs.end())
        {
            m_index[count++] = i;
        }
    }
    if (1 == solveType)
    {
        ASSERTL0(m_motionDofs == 2,
                 "2 Dofs if body is free only in x direction.");
    }
    if (m_motionDofs)
    {
        Array<OneD, NekDouble> temp;
        m_M = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        m_C = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        m_K = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        DNekMatSharedPtr inverseMatrix =
            MemoryManager<DNekMat>::AllocateSharedPtr(m_motionDofs,
                                                      m_motionDofs, 0.0, eFULL);
        for (int i = 0; i < m_motionDofs; ++i)
        {
            m_M[i]     = Array<OneD, NekDouble>(m_rows, 0.);
            m_C[i]     = Array<OneD, NekDouble>(m_rows, 0.);
            m_K[i]     = Array<OneD, NekDouble>(m_rows, 0.);
            int offset = m_index[i] * m_rows;
            for (int j = 0; j < m_rows; ++j)
            {
                int ind   = offset + m_index[j];
                m_M[i][j] = M[ind];
                m_C[i][j] = C[ind];
                m_K[i][j] = K[ind];
                NekDouble value =
                    m_coeffs[0] * M[ind] + C[ind] + m_coeffs[2] * K[ind];
                if (j < m_motionDofs)
                {
                    inverseMatrix->SetValue(i, j, value);
                }
            }
        }

        m_Matrix = Array<OneD, Array<OneD, NekDouble>>(m_motionDofs);
        if (solveType == 0)
        {
            inverseMatrix->Invert();
        }
        for (int i = 0; i < m_motionDofs; ++i)
        {
            m_Matrix[i] = Array<OneD, NekDouble>(m_motionDofs, 0.);
            for (int j = 0; j < m_rows; ++j)
            {
                if (j < m_motionDofs)
                {
                    m_Matrix[i][j] = inverseMatrix->GetValue(i, j);
                }
            }
        }
    }
}

void Newmark_BetaSolver::SolvePrescribed(
    Array<OneD, Array<OneD, NekDouble>> u,
    std::map<int, NekDouble> motionPrescribed)
{
    for (int i = 0; i < m_rows; ++i)
    {
        if (motionPrescribed.find(i) != motionPrescribed.end())
        {
            int i0       = i + m_rows;
            int i2       = i0 + m_rows;
            NekDouble bm = 0., bk = 0.;
            if (motionPrescribed.find(i2) == motionPrescribed.end())
            {
                bm = m_coeffs[0] * u[1][i] + m_coeffs[1] * u[2][i];
            }
            if (motionPrescribed.find(i0) == motionPrescribed.end())
            {
                bk = u[0][i] + m_coeffs[3] * u[1][i] + m_coeffs[4] * u[2][i];
            }

            u[1][i] = motionPrescribed[i];
            if (motionPrescribed.find(i2) == motionPrescribed.end())
            {
                u[2][i] = m_coeffs[0] * u[1][i] - bm;
            }
            else
            {
                u[2][i] = motionPrescribed[i2];
            }
            if (motionPrescribed.find(i0) == motionPrescribed.end())
            {
                u[0][i] = m_coeffs[2] * u[1][i] + bk;
            }
            else
            {
                u[0][i] = motionPrescribed[i0];
            }
        }
    }
}

void Newmark_BetaSolver::SolveFreeFixMat(Array<OneD, Array<OneD, NekDouble>> u,
                                         Array<OneD, NekDouble> force)
{
    if (m_motionDofs)
    {
        Array<OneD, NekDouble> bm(m_motionDofs, 0.);
        Array<OneD, NekDouble> bk(m_motionDofs, 0.);
        for (int j = 0; j < m_motionDofs; ++j)
        {
            int j1 = m_index[j];
            bm[j]  = m_coeffs[0] * u[1][j1] + m_coeffs[1] * u[2][j1];
            bk[j]  = u[0][j1] + m_coeffs[3] * u[1][j1] + m_coeffs[4] * u[2][j1];
        }
        Array<OneD, NekDouble> rhs(m_motionDofs, 0.);
        for (int i = 0; i < m_motionDofs; ++i)
        {
            rhs[i] = force[m_index[i]];
            for (int j = 0; j < m_motionDofs; ++j)
            {
                rhs[i] += m_M[i][j] * bm[j] - m_K[i][j] * bk[j];
            }
            for (int j = m_motionDofs; j < m_rows; ++j)
            {
                int j1 = m_index[j];
                rhs[i] -= m_M[i][j] * u[2][j1] + m_C[i][j] * u[1][j1] +
                          m_K[i][j] * u[0][j1];
            }
        }
        for (int j = 0; j < m_motionDofs; ++j)
        {
            int j1   = m_index[j];
            u[1][j1] = Vmath::Dot(m_motionDofs, m_Matrix[j], 1, rhs, 1);
            u[0][j1] = m_coeffs[2] * u[1][j1] + bk[j];
            u[2][j1] = m_coeffs[0] * u[1][j1] - bm[j];
        }
    }
}

void Newmark_BetaSolver::SolveFreeVarMat6DoF(
    Array<OneD, Array<OneD, NekDouble>> u,
    const Array<OneD, NekDouble> &force,
    const Array<OneD, NekDouble> &nonlinearTerm,
    const Array<OneD, NekDouble> &nonlinearJacobian,
    const Array<OneD, NekDouble> &linearisationVelocity)
{
    SolveFreeVarMatNDof(u, force, nonlinearTerm, nonlinearJacobian,
                         linearisationVelocity, 6);
}

void Newmark_BetaSolver::SolveFreeVarMatNDof(
    Array<OneD, Array<OneD, NekDouble>> u,
    const Array<OneD, NekDouble> &force,
    const Array<OneD, NekDouble> &nonlinearTerm,
    const Array<OneD, NekDouble> &nonlinearJacobian,
    const Array<OneD, NekDouble> &linearisationVelocity, const int nDofs)
{
    ASSERTL0(m_motionDofs == nDofs,
             "The nonlinear moving-frame solver has inconsistent DoFs.");
    ASSERTL0(force.size() >= m_rows && nonlinearTerm.size() >= m_rows &&
                 nonlinearJacobian.size() >= m_rows * m_rows &&
                 linearisationVelocity.size() >= m_rows,
             "Invalid nonlinear rigid-body system size.");

    Array<OneD, NekDouble> bm(nDofs, 0.0);
    Array<OneD, NekDouble> bk(nDofs, 0.0);
    Array<OneD, NekDouble> rhs(nDofs, 0.0);
    std::vector<double> matrix(nDofs * nDofs, 0.0);
    std::vector<double> drhs(nDofs, 0.0);
    std::vector<int> ipiv(nDofs, 0);

    for (int j = 0; j < nDofs; ++j)
    {
        const int j1 = m_index[j];
        bm[j] = m_coeffs[0] * u[1][j1] + m_coeffs[1] * u[2][j1];
        bk[j] = u[0][j1] + m_coeffs[3] * u[1][j1] + m_coeffs[4] * u[2][j1];
    }

    for (int i = 0; i < nDofs; ++i)
    {
        const int i1 = m_index[i];
        rhs[i]       = force[i1] - nonlinearTerm[i1];
        for (int j = 0; j < nDofs; ++j)
        {
            const int j1 = m_index[j];
            rhs[i] += nonlinearJacobian[i1 + j1 * m_rows] *
                      linearisationVelocity[j1];
            rhs[i] += m_M[i][j] * bm[j] - m_K[i][j] * bk[j];
            matrix[j * nDofs + i] =
                m_coeffs[0] * m_M[i][j] + m_C[i][j] +
                m_coeffs[2] * m_K[i][j] +
                nonlinearJacobian[i1 + j1 * m_rows];
        }
        for (int j = nDofs; j < m_rows; ++j)
        {
            const int j1 = m_index[j];
            rhs[i] -= m_M[i][j] * u[2][j1] + m_C[i][j] * u[1][j1] +
                      m_K[i][j] * u[0][j1];
        }
        drhs[i] = rhs[i];
    }

    int info = 0;
    Lapack::DoSgetrf(nDofs, nDofs, matrix.data(), nDofs, ipiv.data(), info);
    ASSERTL0(info == 0,
             "Singular nonlinear rigid-body Newmark effective matrix.");
    Lapack::Dgetrs('N', nDofs, 1, matrix.data(), nDofs, ipiv.data(),
                   drhs.data(), nDofs, info);
    ASSERTL0(info == 0, "Failed to solve nonlinear rigid-body Newmark system.");

    for (int j = 0; j < nDofs; ++j)
    {
        const int j1 = m_index[j];
        u[1][j1]     = drhs[j];
        u[0][j1]     = m_coeffs[2] * u[1][j1] + bk[j];
        u[2][j1]     = m_coeffs[0] * u[1][j1] - bm[j];
    }
}

void Newmark_BetaSolver::SolveFreeVarMatNDofConstrained(
    Array<OneD, Array<OneD, NekDouble>> u,
    const Array<OneD, NekDouble> &force,
    const Array<OneD, NekDouble> &nonlinearTerm,
    const Array<OneD, NekDouble> &nonlinearJacobian,
    const Array<OneD, NekDouble> &linearisationVelocity,
    const Array<OneD, NekDouble> &velocityConstraints,
    const Array<OneD, NekDouble> &constraintVelocity, const int nDofs,
    const int nConstraints)
{
    ASSERTL0(m_motionDofs == nDofs && nConstraints > 0,
             "Invalid constrained nonlinear rigid-body system size.");
    ASSERTL0(force.size() >= m_rows && nonlinearTerm.size() >= m_rows &&
                 nonlinearJacobian.size() >= m_rows * m_rows &&
                 linearisationVelocity.size() >= m_rows &&
                 velocityConstraints.size() >= nConstraints * m_rows &&
                 constraintVelocity.size() >= nConstraints,
             "Invalid constrained nonlinear rigid-body system data.");

    Array<OneD, NekDouble> bm(nDofs, 0.0);
    Array<OneD, NekDouble> bk(nDofs, 0.0);
    const int systemSize = nDofs + nConstraints;
    std::vector<double> matrix(systemSize * systemSize, 0.0);
    std::vector<double> rhs(systemSize, 0.0);
    std::vector<int> ipiv(systemSize, 0);

    for (int j = 0; j < nDofs; ++j)
    {
        const int j1 = m_index[j];
        bm[j] = m_coeffs[0] * u[1][j1] + m_coeffs[1] * u[2][j1];
        bk[j] = u[0][j1] + m_coeffs[3] * u[1][j1] + m_coeffs[4] * u[2][j1];
    }

    for (int i = 0; i < nDofs; ++i)
    {
        const int i1 = m_index[i];
        rhs[i]       = force[i1] - nonlinearTerm[i1];
        for (int j = 0; j < nDofs; ++j)
        {
            const int j1 = m_index[j];
            rhs[i] += nonlinearJacobian[i1 + j1 * m_rows] *
                      linearisationVelocity[j1];
            rhs[i] += m_M[i][j] * bm[j] - m_K[i][j] * bk[j];
            matrix[j * systemSize + i] =
                m_coeffs[0] * m_M[i][j] + m_C[i][j] +
                m_coeffs[2] * m_K[i][j] +
                nonlinearJacobian[i1 + j1 * m_rows];
        }
        for (int j = nDofs; j < m_rows; ++j)
        {
            const int j1 = m_index[j];
            rhs[i] -= m_M[i][j] * u[2][j1] + m_C[i][j] * u[1][j1] +
                      m_K[i][j] * u[0][j1];
        }
        for (int k = 0; k < nConstraints; ++k)
        {
            const NekDouble value = velocityConstraints[k * m_rows + i1];
            matrix[(nDofs + k) * systemSize + i] = value;
            matrix[i * systemSize + nDofs + k] = value;
        }
    }
    for (int k = 0; k < nConstraints; ++k)
    {
        rhs[nDofs + k] = constraintVelocity[k];
    }

    int info = 0;
    Lapack::DoSgetrf(systemSize, systemSize, matrix.data(), systemSize,
                     ipiv.data(), info);
    ASSERTL0(info == 0,
             "Singular constrained nonlinear rigid-body Newmark matrix.");
    Lapack::Dgetrs('N', systemSize, 1, matrix.data(), systemSize, ipiv.data(),
                   rhs.data(), systemSize, info);
    ASSERTL0(info == 0,
             "Failed to solve constrained nonlinear rigid-body Newmark system.");

    for (int j = 0; j < nDofs; ++j)
    {
        const int j1 = m_index[j];
        u[1][j1]     = rhs[j];
        u[0][j1]     = m_coeffs[2] * u[1][j1] + bk[j];
        u[2][j1]     = m_coeffs[0] * u[1][j1] - bm[j];
    }
}

/**
 * e_x \cdot M [du0, du1, du2]  + m (e_x \times Omega) \cdot u = F \cdot e_x
 * e_x = (c, -s, 0)
 * (c M00 - s M10, c M01 - s M11, c M02 - s M 12) [du0, du1, du2]^T +
 * (-s Omega, -c Omega, 0) [u0, u1, u2]^T =
 * c F0 - s F1
 * equation 2, e_y \cdot (u0, u1) = uy = s * u0 + c * u1
 **/
void Newmark_BetaSolver::SolveOneFree(Array<OneD, Array<OneD, NekDouble>> u,
                                      Array<OneD, NekDouble> force,
                                      const Array<OneD, NekDouble> theta,
                                      const NekDouble uy, const NekDouble mass)
{
    int iOmega  = m_rows - 1;
    NekDouble c = cos(theta[2]), s = sin(theta[2]);
    NekDouble C00, C01, C10, C11, M0, M1, M2, F0, F1;
    M0  = c * m_M[0][0] - s * m_M[1][0];
    M1  = c * m_M[0][1] - s * m_M[1][1];
    M2  = c * m_M[0][iOmega] - s * m_M[1][iOmega];
    C00 = m_coeffs[0] * M0 - s * u[1][iOmega] * mass;
    C01 = m_coeffs[0] * M1 - c * u[1][iOmega] * mass;
    C10 = s;
    C11 = c;
    F1  = uy;
    // F0
    Array<OneD, NekDouble> bm(m_motionDofs, 0.);
    Array<OneD, NekDouble> bk(m_motionDofs, 0.);
    for (int j = 0; j < m_motionDofs; ++j)
    {
        bm[j] = m_coeffs[0] * u[1][j] + m_coeffs[1] * u[2][j];
        bk[j] = u[0][j] + m_coeffs[3] * u[1][j] + m_coeffs[4] * u[2][j];
    }
    F0 = c * force[0] - s * force[1];
    F0 = F0 + M0 * bm[0] + M1 * bm[1] - M2 * u[2][iOmega];
    // solve
    NekDouble det = 1. / (C00 * C11 - C01 * C10);
    u[1][0]       = det * (C11 * F0 - C01 * F1);
    u[1][1]       = det * (-C10 * F0 + C00 * F1);
    for (int j = 0; j < m_motionDofs; ++j)
    {
        u[0][j] = m_coeffs[2] * u[1][j] + bk[j];
        u[2][j] = m_coeffs[0] * u[1][j] - bm[j];
    }
}

FrameTransform::FrameTransform()
{
    m_matrix = Array<OneD, NekDouble>(2, 0.);
}

void FrameTransform::SetAngle(const Array<OneD, NekDouble> theta)
{
    m_matrix[0] = cos(theta[2]);
    m_matrix[1] = sin(theta[2]);
}

void FrameTransform::BodyToInerital(const int dim,
                                    const Array<OneD, NekDouble> &body,
                                    Array<OneD, NekDouble> &inertial)
{
    NekDouble xi = body[0], eta = body[1];
    inertial[0] = xi * m_matrix[0] - eta * m_matrix[1];
    inertial[1] = eta * m_matrix[0] + xi * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            inertial[i] = body[i];
        }
    }
}

void FrameTransform::IneritalToBody(const int dim,
                                    const Array<OneD, NekDouble> &inertial,
                                    Array<OneD, NekDouble> &body)
{
    NekDouble x = inertial[0], y = inertial[1];
    body[0] = x * m_matrix[0] + y * m_matrix[1];
    body[1] = y * m_matrix[0] - x * m_matrix[1];
    if (body != inertial && dim > 2)
    {
        for (int i = 2; i < dim; ++i)
        {
            body[i] = inertial[i];
        }
    }
}

} // namespace Nektar
