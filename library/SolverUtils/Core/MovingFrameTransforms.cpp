///////////////////////////////////////////////////////////////////////////////
//
// File: MovingFrameTransforms.cpp
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

#include <SolverUtils/Core/MovingFrameTransforms.h>

#include <LibUtilities/BasicUtils/ErrorUtil.hpp>

#include <algorithm>
#include <cmath>

namespace Nektar::SolverUtils::MovingFrame
{

namespace
{

void MultiplyByPureQuaternion(const Array<OneD, NekDouble> &q,
                              const Array<OneD, NekDouble> &v,
                              Array<OneD, NekDouble> &out)
{
    out[0] = -(q[1] * v[0] + q[2] * v[1] + q[3] * v[2]);
    out[1] = q[0] * v[0] + q[2] * v[2] - q[3] * v[1];
    out[2] = q[0] * v[1] + q[3] * v[0] - q[1] * v[2];
    out[3] = q[0] * v[2] + q[1] * v[1] - q[2] * v[0];
}

Array<OneD, NekDouble> BodyVectorFromQuaternionDerivative(
    const Array<OneD, NekDouble> &qin,
    const Array<OneD, NekDouble> &dq)
{
    const Array<OneD, NekDouble> q = NormalizeQuaternion(qin);
    Array<OneD, NekDouble> out(3, 0.0);
    out[0] = 2.0 * (q[0] * dq[1] - q[1] * dq[0] - q[2] * dq[3] +
                    q[3] * dq[2]);
    out[1] = 2.0 * (q[0] * dq[2] + q[1] * dq[3] - q[2] * dq[0] -
                    q[3] * dq[1]);
    out[2] = 2.0 * (q[0] * dq[3] - q[1] * dq[2] + q[2] * dq[1] -
                    q[3] * dq[0]);
    return out;
}

} // namespace

Array<OneD, NekDouble> IdentityQuaternion()
{
    Array<OneD, NekDouble> q(4, 0.0);
    q[0] = 1.0;
    return q;
}

bool TryNormalizeQuaternion(const Array<OneD, NekDouble> &q,
                             Array<OneD, NekDouble> &qout)
{
    const NekDouble q0    = q[0];
    const NekDouble q1    = q[1];
    const NekDouble q2    = q[2];
    const NekDouble q3    = q[3];
    const NekDouble norm2 = q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3;
    if (!std::isfinite(norm2) ||
        norm2 <= NekConstants::kNekZeroTol * NekConstants::kNekZeroTol)
    {
        return false;
    }

    const NekDouble norm = std::sqrt(norm2);
    if (!std::isfinite(norm))
    {
        return false;
    }

    qout = Array<OneD, NekDouble>(4, 0.0);
    qout[0] = q0 / norm;
    qout[1] = q1 / norm;
    qout[2] = q2 / norm;
    qout[3] = q3 / norm;
    return true;
}

Array<OneD, NekDouble> NormalizeQuaternion(const Array<OneD, NekDouble> &q)
{
    Array<OneD, NekDouble> qout;
    ASSERTL0(TryNormalizeQuaternion(q, qout),
             "Moving reference frame quaternion has invalid norm.");

    return qout;
}

Array<OneD, NekDouble> QuaternionFromEulerZYX(NekDouble thetaX,
                                               NekDouble thetaY,
                                               NekDouble thetaZ)
{
    const NekDouble cr = std::cos(0.5 * thetaX);
    const NekDouble sr = std::sin(0.5 * thetaX);
    const NekDouble cp = std::cos(0.5 * thetaY);
    const NekDouble sp = std::sin(0.5 * thetaY);
    const NekDouble cy = std::cos(0.5 * thetaZ);
    const NekDouble sy = std::sin(0.5 * thetaZ);

    Array<OneD, NekDouble> q(4, 0.0);
    q[0] = cy * cp * cr + sy * sp * sr;
    q[1] = cy * cp * sr - sy * sp * cr;
    q[2] = sy * cp * sr + cy * sp * cr;
    q[3] = sy * cp * cr - cy * sp * sr;
    return NormalizeQuaternion(q);
}

Array<OneD, NekDouble> QuaternionFromAxisAngle(int axis, NekDouble theta)
{
    ASSERTL0(axis >= 0 && axis < 3,
             "Moving reference frame rotation axis must be x, y or z.");

    Array<OneD, NekDouble> q(4, 0.0);
    q[0]          = std::cos(0.5 * theta);
    q[axis + 1]   = std::sin(0.5 * theta);
    return NormalizeQuaternion(q);
}

Array<OneD, NekDouble> MultiplyQuaternions(
    const Array<OneD, NekDouble> &qain,
    const Array<OneD, NekDouble> &qbin)
{
    const Array<OneD, NekDouble> qa = NormalizeQuaternion(qain);
    const Array<OneD, NekDouble> qb = NormalizeQuaternion(qbin);

    Array<OneD, NekDouble> q(4, 0.0);
    q[0] = qa[0] * qb[0] - qa[1] * qb[1] - qa[2] * qb[2] - qa[3] * qb[3];
    q[1] = qa[0] * qb[1] + qa[1] * qb[0] + qa[2] * qb[3] - qa[3] * qb[2];
    q[2] = qa[0] * qb[2] - qa[1] * qb[3] + qa[2] * qb[0] + qa[3] * qb[1];
    q[3] = qa[0] * qb[3] + qa[1] * qb[2] - qa[2] * qb[1] + qa[3] * qb[0];
    return NormalizeQuaternion(q);
}

Array<OneD, NekDouble> EulerZYXFromQuaternion(
    const Array<OneD, NekDouble> &qin)
{
    const Array<OneD, NekDouble> q = NormalizeQuaternion(qin);

    const NekDouble sinr =
        2.0 * (q[0] * q[1] + q[2] * q[3]);
    const NekDouble cosr =
        1.0 - 2.0 * (q[1] * q[1] + q[2] * q[2]);

    const NekDouble sinp =
        2.0 * (q[0] * q[2] - q[3] * q[1]);

    const NekDouble siny =
        2.0 * (q[0] * q[3] + q[1] * q[2]);
    const NekDouble cosy =
        1.0 - 2.0 * (q[2] * q[2] + q[3] * q[3]);

    Array<OneD, NekDouble> theta(3, 0.0);
    theta[0] = std::atan2(sinr, cosr);
    theta[1] = std::asin(std::max(-1.0, std::min(1.0, sinp)));
    theta[2] = std::atan2(siny, cosy);
    return theta;
}

Array<OneD, NekDouble> IntegrateQuaternionBodyOmega(
    const Array<OneD, NekDouble> &qin,
    const Array<OneD, NekDouble> &omegaBody, NekDouble dt)
{
    const Array<OneD, NekDouble> q = NormalizeQuaternion(qin);
    const NekDouble omegaNorm = std::sqrt(
        omegaBody[0] * omegaBody[0] + omegaBody[1] * omegaBody[1] +
        omegaBody[2] * omegaBody[2]);
    if (omegaNorm <= NekConstants::kNekZeroTol || dt == 0.0)
    {
        return q;
    }

    // Exact Lie-group update for a constant body-frame angular velocity over
    // this step. When omegaBody is evaluated at the time-step midpoint, this
    // is the second-order midpoint update of qdot = 0.5 q otimes omegaBody.
    const NekDouble halfAngle = 0.5 * omegaNorm * dt;
    Array<OneD, NekDouble> delta(4, 0.0);
    delta[0] = std::cos(halfAngle);
    const NekDouble scale = std::sin(halfAngle) / omegaNorm;
    for (int i = 0; i < 3; ++i)
    {
        delta[i + 1] = scale * omegaBody[i];
    }
    return MultiplyQuaternions(q, delta);
}

Array<OneD, NekDouble> MakeQuaternionContinuous(
    const Array<OneD, NekDouble> &qref, const Array<OneD, NekDouble> &qin)
{
    const Array<OneD, NekDouble> q0 = NormalizeQuaternion(qref);
    Array<OneD, NekDouble> q        = NormalizeQuaternion(qin);
    NekDouble dot                  = 0.0;
    for (int i = 0; i < 4; ++i)
    {
        dot += q0[i] * q[i];
    }
    if (dot < 0.0)
    {
        for (int i = 0; i < 4; ++i)
        {
            q[i] = -q[i];
        }
    }
    return q;
}

Array<OneD, NekDouble> QuaternionDerivativeFromBodyOmega(
    const Array<OneD, NekDouble> &qin,
    const Array<OneD, NekDouble> &omegaBody)
{
    const Array<OneD, NekDouble> q = NormalizeQuaternion(qin);
    Array<OneD, NekDouble> qdot(4, 0.0);
    MultiplyByPureQuaternion(q, omegaBody, qdot);
    for (int i = 0; i < 4; ++i)
    {
        qdot[i] *= 0.5;
    }
    return qdot;
}

Array<OneD, NekDouble> QuaternionSecondDerivativeFromBodyOmegaAlpha(
    const Array<OneD, NekDouble> &qin,
    const Array<OneD, NekDouble> &omegaBody,
    const Array<OneD, NekDouble> &alphaBody)
{
    const Array<OneD, NekDouble> q     = NormalizeQuaternion(qin);
    const Array<OneD, NekDouble> qdot =
        QuaternionDerivativeFromBodyOmega(q, omegaBody);
    Array<OneD, NekDouble> qdotOmega(4, 0.0);
    Array<OneD, NekDouble> qAlpha(4, 0.0);
    Array<OneD, NekDouble> qddot(4, 0.0);
    MultiplyByPureQuaternion(qdot, omegaBody, qdotOmega);
    MultiplyByPureQuaternion(q, alphaBody, qAlpha);
    for (int i = 0; i < 4; ++i)
    {
        qddot[i] = 0.5 * (qdotOmega[i] + qAlpha[i]);
    }
    return qddot;
}

Array<OneD, NekDouble> BodyOmegaFromQuaternionDerivative(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &qdot)
{
    return BodyVectorFromQuaternionDerivative(q, qdot);
}

Array<OneD, NekDouble> BodyAlphaFromQuaternionSecondDerivative(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &qddot)
{
    return BodyVectorFromQuaternionDerivative(q, qddot);
}

void Cross(const Array<OneD, NekDouble> &a, const Array<OneD, NekDouble> &b,
           Array<OneD, NekDouble> &out)
{
    const NekDouble x = a[1] * b[2] - a[2] * b[1];
    const NekDouble y = a[2] * b[0] - a[0] * b[2];
    const NekDouble z = a[0] * b[1] - a[1] * b[0];
    out[0]            = x;
    out[1]            = y;
    out[2]            = z;
}

Array<OneD, NekDouble> RotateBodyToInertial(
    const Array<OneD, NekDouble> &qin, const Array<OneD, NekDouble> &v)
{
    const Array<OneD, NekDouble> q = NormalizeQuaternion(qin);
    Array<OneD, NekDouble> qv(3, 0.0);
    Array<OneD, NekDouble> t(3, 0.0);
    Array<OneD, NekDouble> u(3, 0.0);
    Array<OneD, NekDouble> out(3, 0.0);

    qv[0] = q[1];
    qv[1] = q[2];
    qv[2] = q[3];
    Cross(qv, v, t);
    Cross(qv, t, u);

    for (int i = 0; i < 3; ++i)
    {
        out[i] = v[i] + 2.0 * (q[0] * t[i] + u[i]);
    }
    return out;
}

Array<OneD, NekDouble> RotateInertialToBody(
    const Array<OneD, NekDouble> &q, const Array<OneD, NekDouble> &v)
{
    Array<OneD, NekDouble> qconj(4, 0.0);
    qconj[0] = q[0];
    qconj[1] = -q[1];
    qconj[2] = -q[2];
    qconj[3] = -q[3];
    return RotateBodyToInertial(qconj, v);
}

bool QuaternionFromParams(const std::map<std::string, NekDouble> &params,
                           Array<OneD, NekDouble> &q)
{
    const bool hasQ = params.find("Q0") != params.end() ||
                      params.find("Q1") != params.end() ||
                      params.find("Q2") != params.end() ||
                      params.find("Q3") != params.end();
    if (hasQ)
    {
        q = IdentityQuaternion();
        const std::string names[4] = {"Q0", "Q1", "Q2", "Q3"};
        for (int i = 0; i < 4; ++i)
        {
            auto it = params.find(names[i]);
            if (it != params.end())
            {
                q[i] = it->second;
            }
        }
        q = NormalizeQuaternion(q);
        return true;
    }

    const bool hasTheta = params.find("Theta_x") != params.end() ||
                          params.find("Theta_y") != params.end() ||
                          params.find("Theta_z") != params.end();
    if (hasTheta)
    {
        Array<OneD, NekDouble> theta(3, 0.0);
        const std::string names[3] = {"Theta_x", "Theta_y", "Theta_z"};
        for (int i = 0; i < 3; ++i)
        {
            auto it = params.find(names[i]);
            if (it != params.end())
            {
                theta[i] = it->second;
            }
        }
        q = QuaternionFromEulerZYX(theta[0], theta[1], theta[2]);
        return true;
    }
    return false;
}

Array<OneD, NekDouble> QuaternionFromFrameData(
    const Array<OneD, NekDouble> &frameData)
{
    if (frameData.size() >= kFrameDataSizeWithQuat)
    {
        Array<OneD, NekDouble> qraw(4, 0.0);
        Array<OneD, NekDouble> qout;
        for (int i = 0; i < 4; ++i)
        {
            qraw[i] = frameData[kQuaternionOffset + i];
        }
        if (TryNormalizeQuaternion(qraw, qout))
        {
            return qout;
        }
    }

    if (frameData.size() >= 6)
    {
        return QuaternionFromEulerZYX(frameData[3], frameData[4],
                                      frameData[5]);
    }

    return IdentityQuaternion();
}

void StoreQuaternionInFrameData(const Array<OneD, NekDouble> &q,
                                 Array<OneD, NekDouble> &frameData)
{
    if (frameData.size() >= kFrameDataSizeWithQuat)
    {
        const Array<OneD, NekDouble> qn = NormalizeQuaternion(q);
        for (int i = 0; i < 4; ++i)
        {
            frameData[kQuaternionOffset + i] = qn[i];
        }
    }
}

void StoreQuaternionInCompactDisp(const Array<OneD, NekDouble> &q,
                                  Array<OneD, NekDouble> &frameDisp)
{
    if (frameDisp.size() >= kQuaternionCompactOffset + 4)
    {
        const Array<OneD, NekDouble> qn = NormalizeQuaternion(q);
        for (int i = 0; i < 4; ++i)
        {
            frameDisp[kQuaternionCompactOffset + i] = qn[i];
        }
    }
}

} // namespace Nektar::SolverUtils::MovingFrame
