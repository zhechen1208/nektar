////////////////////////////////////////////////////////////////////////////////
//
//  File: ProcessBVF.cpp
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: Computes boundary vorticity flux.
//  input: u, v, w, u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z
//  input: p, p_x, p_y
//  input: u_xx, u_yy, v_xx, v_yy
////////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <string>
using namespace std;

#include <GlobalMapping/Mapping.h>
#include <LibUtilities/BasicUtils/ParseUtils.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>

#include "ProcessBVF.h"
#include "ProcessMapping.h"

namespace Nektar::FieldUtils
{

ModuleKey ProcessBVF::className = GetModuleFactory().RegisterCreatorFunction(
    ModuleKey(eProcessModule, "BVF"), ProcessBVF::create,
    "Computes boundary vorticity flux.");

ProcessBVF::ProcessBVF(FieldSharedPtr f) : ProcessModule(f)
{
}

ProcessBVF::~ProcessBVF()
{
}

void ProcessBVF::v_Process(po::variables_map &vm)
{
    m_f->SetUpExp(vm);

    int nfields   = m_f->m_variables.size();
    int addfields = 1;
    m_f->m_exp.resize(nfields + addfields);
    m_f->m_variables.push_back("uu_x");

    // Skip in case of empty partition
    if (m_f->m_exp[0]->GetNumElmts() == 0)
    {
        return;
    }

    int npoints = m_f->m_exp[0]->GetNpoints();
    Array<OneD, NekDouble> grad(npoints, 0.);
    set<int> uset, uxset;
    ParseUtils::GenerateVariableSet("u", m_f->m_variables, uset);
    ParseUtils::GenerateVariableSet("u_x", m_f->m_variables, uxset);
    ASSERTL0(!uset.empty() && !uxset.empty(), "u or u_x not found.");
    int uid  = *uset.begin();
    int uxid = *uxset.begin();
    Vmath::Vmul(npoints, m_f->m_exp[uid]->GetPhys(), 1,
                m_f->m_exp[uxid]->GetPhys(), 1, grad, 1);

    m_f->m_exp[nfields] = m_f->AppendExpList(m_f->m_numHomogeneousDir);
    Vmath::Vcopy(npoints, grad, 1, m_f->m_exp[nfields]->UpdatePhys(), 1);
    m_f->m_exp[nfields]->FwdTransLocalElmt(grad,
                                           m_f->m_exp[nfields]->UpdateCoeffs());
}
} // namespace Nektar::FieldUtils
