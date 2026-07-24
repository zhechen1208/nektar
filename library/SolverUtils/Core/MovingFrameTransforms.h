///////////////////////////////////////////////////////////////////////////////
//
// File: MovingFrameTransforms.h
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
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERUTILS_CORE_MOVINGFRAMETRANSFORMS_H
#define NEKTAR_SOLVERUTILS_CORE_MOVINGFRAMETRANSFORMS_H

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <SolverUtils/SolverUtilsDeclspec.h>

#include <map>
#include <string>

namespace Nektar::SolverUtils::MovingFrame
{

static constexpr int kLegacyFrameDataSize     = 21;
static constexpr int kQuaternionOffset        = 21;
static constexpr int kFrameDataSizeWithQuat   = 25;
static constexpr int kQuaternionCompactOffset = 6;

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> IdentityQuaternion();

SOLVER_UTILS_EXPORT bool TryNormalizeQuaternion(
    const Array<OneD, NekDouble> &q, Array<OneD, NekDouble> &qout);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> NormalizeQuaternion(
    const Array<OneD, NekDouble> &q);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> QuaternionFromEulerZYX(
    NekDouble thetaX, NekDouble thetaY, NekDouble thetaZ);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> QuaternionFromAxisAngle(
    int axis, NekDouble theta);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> MultiplyQuaternions(
    const Array<OneD, NekDouble> &qa, const Array<OneD, NekDouble> &qb);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> EulerZYXFromQuaternion(
    const Array<OneD, NekDouble> &q);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> IntegrateQuaternionBodyOmega(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &omegaBody,
    NekDouble dt);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> MakeQuaternionContinuous(
    const Array<OneD, NekDouble> &qref, const Array<OneD, NekDouble> &q);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> QuaternionDerivativeFromBodyOmega(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &omegaBody);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble>
QuaternionSecondDerivativeFromBodyOmegaAlpha(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &omegaBody,
    const Array<OneD, NekDouble> &alphaBody);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> BodyOmegaFromQuaternionDerivative(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &qdot);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble>
BodyAlphaFromQuaternionSecondDerivative(const Array<OneD, NekDouble> &q,
                                        const Array<OneD, NekDouble> &qddot);

SOLVER_UTILS_EXPORT void Cross(const Array<OneD, NekDouble> &a,
                               const Array<OneD, NekDouble> &b,
                               Array<OneD, NekDouble> &out);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> RotateBodyToInertial(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &v);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> RotateInertialToBody(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &v);

SOLVER_UTILS_EXPORT bool QuaternionFromParams(
    const std::map<std::string, NekDouble> &params, Array<OneD, NekDouble> &q);

SOLVER_UTILS_EXPORT Array<OneD, NekDouble> QuaternionFromFrameData(
    const Array<OneD, NekDouble> &frameData);

SOLVER_UTILS_EXPORT void StoreQuaternionInFrameData(
    const Array<OneD, NekDouble> &q, Array<OneD, NekDouble> &frameData);

SOLVER_UTILS_EXPORT void StoreQuaternionInCompactDisp(
    const Array<OneD, NekDouble> &q, Array<OneD, NekDouble> &frameDisp);

} // namespace Nektar::SolverUtils::MovingFrame

#endif
