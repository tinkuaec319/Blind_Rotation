//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2022, NJIT, Duality Technologies Inc. and other contributors
//
// All rights reserved.
//
// Author TPOC: contact@openfhe.org
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================

/*
 * Custom Modifications:
 * - [This code is the implementation of the algorithm in the paper https://eprint.iacr.org/2023/1564]
 * 
 * This modified section follows the terms of the original BSD 2-Clause License.
 * Other modifications are provided under the terms of the BSD 2-Clause License.
 * See the BSD 2-Clause License text below:
 */

//==================================================================================
// Additional BSD License for Custom Modifications:
//
// Copyright (c) 2023 Binwu Xiang,Kaixing Wang and other contributors
//
// All rights reserved.
//
// Author TPOC: wangkaixing22@mails.ucas.ac.cn
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//==================================================================================
#include <ostream>

#include "binfhe-constants.h"

namespace lbcrypto {

std::ostream& operator<<(std::ostream& s, BINFHE_PARAMSET f) {
    switch (f) {
        case TOY:
            s << "TOY";
            break;
        case MEDIUM:
            s << "MEDIUM";
            break;
        case STD128_LMKCDEY:
            s << "STD128_LMKCDEY";
            break;
        case STD128_AP:
            s << "STD128_AP";
            break;
        case STD128:
            s << "STD128";
            break;
        case STD192:
            s << "STD192";
            break;
        case STD256:
            s << "STD256";
            break;
        case STD128Q:
            s << "STD128Q";
            break;
        case STD128Q_LMKCDEY:
            s << "STD128Q_LMKCDEY";
            break;
        case STD192Q:
            s << "STD192Q";
            break;
        case STD256Q:
            s << "STD256Q";
            break;
        case STD128_3:
            s << "STD128_3";
            break;
        case STD128_3_LMKCDEY:
            s << "STD128_3_LMKCDEY";
            break;
        case STD128Q_3:
            s << "STD128Q_3";
            break;
        case STD128Q_3_LMKCDEY:
            s << "STD128Q_3_LMKCDEY";
            break;
        case STD192Q_3:
            s << "STD192Q_3";
            break;
        case STD256Q_3:
            s << "STD256Q_3";
            break;
        case STD128_4:
            s << "STD128_4";
            break;
        case STD128_4_LMKCDEY:
            s << "STD128_4_LMKCDEY";
            break;
        case STD128Q_4:
            s << "STD128Q_4";
            break;
        case STD128Q_4_LMKCDEY:
            s << "STD128Q_4_LMKCDEY";
            break;
        case STD192Q_4:
            s << "STD192Q_4";
            break;
        case STD256Q_4:
            s << "STD256Q_4";
            break;
        case SIGNED_MOD_TEST:
            s << "SIGNED_MOD_TEST";
            break;
        case P128T:
            s << "P128T";
            break;
        case P128G:
            s << "P128G";
            break;
        case P128T_2:
            s << "P128T_2";
            break;
        case P128G_2:
            s << "P128G_2";
            break;
        case STD128_LMKCDEY_New:
            s << "STD128_LMKCDEY_New";
            break;
        case P192G:
            s << "P192G";
            break;
        case P192T:
            s << "P192T";
            break;

        case STD128_LMKCDEY_LWR:
            s << "STD128_LMKCDEY_LWR";
            break;
        case STD128_LMKCDEY_LWE:
            s << "STD128_LMKCDEY_LWE";
            break;
        case P192G_LWE:
            s << "P192G_LWE";
            break;
        case P192G_LWR:
            s << "P192G_LWR";
            break;                                    

        default:
            s << "UNKNOWN";
            break;
    }
    return s;
}

std::ostream& operator<<(std::ostream& s, BINFHE_OUTPUT f) {
    switch (f) {
        case FRESH:
            s << "FRESH";
            break;
        case BOOTSTRAPPED:
            s << "BOOTSTRAPPED";
            break;
        default:
            s << "UNKNOWN";
            break;
    }
    return s;
}

std::ostream& operator<<(std::ostream& s, BINFHE_METHOD f) {
    switch (f) {
        case AP:
            s << "DM";
            break;
        case GINX:
            s << "CGGI";
            break;
        case LMKCDEY:
            s << "LMKCDEY";
            break;
        case XZDDF:
            s << "XZDDF";
            break;
        default:
            s << "UNKNOWN";
            break;
    }
    return s;
}

std::ostream& operator<<(std::ostream& s, BINGATE f) {
    switch (f) {
        case OR:
            s << "OR";
            break;
        case AND:
            s << "AND";
            break;
        case NOR:
            s << "NOR";
            break;
        case NAND:
            s << "NAND";
            break;
        case XOR_FAST:
            s << "XOR_FAST";
            break;
        case XNOR_FAST:
            s << "XNOR_FAST";
            break;
        case XOR:
            s << "XOR";
            break;
        case XNOR:
            s << "XNOR";
            break;
        case AND3:
            s << "AND3";
            break;
        case OR3:
            s << "OR3";
            break;
        case AND4:
            s << "AND4";
            break;
        case OR4:
            s << "OR4";
            break;
        case MAJORITY:
            s << "MAJORITY";
            break;
        case CMUX:
            s << "CMUX";
            break;
        default:
            s << "UNKNOWN";
            break;
    }
    return s;
}

};  // namespace lbcrypto
/*
    NativeInteger temmp, temmpn,diff=0;
    for (size_t i=0;i<n;i++){
        temmp=b_vec[i].first;
        if (temmp%2==0){
            if(i+1<n){
                temmpn = b_vec[i+1].first;
                diff = temmpn - temmp;
                b_vec[i].first =temmp + floor(diff/2) ;
            }
            if(b_vec[i-1].first==temmp-1){
                b_vec[i].first=temmp-1;
            }
            else{
                b_vec[i].first=temmp+1;
            }            
        }
        else{
            if(b_vec[i-1].first==temmp){ 
                    b_vec[i].first=temmp;
                }
            else if (b_vec[i+1].first==temmp+1) {
                    b_vec[i].first=temmp;
            }
            else if(b_vec[i-1].first==temmp-2){ 
                    b_vec[i].first=temmp-2;
                }
            else if (b_vec[i+1].first==temmp+2) {
                    b_vec[i].first=temmp+2;
            }
            else
                b_vec[i].first=temmp;
        }        
    }

    NativeInteger temmp;
    for (size_t i=0;i<n;i++){
        temmp=b_vec[i].first;
        if (temmp%2==0){
            if(b_vec[i-1].first==temmp-1){
                b_vec[i].first=temmp-1;
            }
            else{
                b_vec[i].first=temmp+1;
            }            
        }
        else{
            if(b_vec[i-1].first==temmp){ 
                    b_vec[i].first=temmp;
                }
            else if (b_vec[i+1].first==temmp+1) {
                    b_vec[i].first=temmp;
            }
            else if(b_vec[i-1].first==temmp-2){ 
                    b_vec[i].first=temmp-2;
                }
            else if (b_vec[i+1].first==temmp+2) {
                    b_vec[i].first=temmp+2;
            }
            else
                b_vec[i].first=temmp;
        }        
    }
    */