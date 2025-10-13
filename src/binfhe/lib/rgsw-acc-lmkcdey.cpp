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

#include "rgsw-acc-lmkcdey.h"
#include <bits/stdc++.h>
#include <vector>
#include <string>

namespace lbcrypto {

// Key generation as described in https://eprint.iacr.org/2022/198
RingGSWACCKey RingGSWAccumulatorLMKCDEY::KeyGenAcc(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                   const NativePoly& skNTT, ConstLWEPrivateKey& LWEsk) const {
    auto sv{LWEsk->GetElement()};
    auto mod{sv.GetModulus().ConvertToInt<int32_t>()};
    auto modHalf{mod >> 1};
    uint32_t N{params->GetN()};
    size_t n{sv.GetLength()};
    uint32_t numAutoKeys{params->GetNumAutoKeys()};
    uint32_t total_rgsw_keys=(numAutoKeys+1)*n;
    NativeInteger gen = NativeInteger(5);

    // std::cout<<"Here is the secret key:"<<sv<<std::endl;
    // dim2, 0: for RGSW(X^si), 1: for automorphism keys
    // only w automorphism keys required
    // allocates (n - w) more memory for pointer (not critical for performance)
    RingGSWACCKey ek = std::make_shared<RingGSWACCKeyImpl>(1, 3, total_rgsw_keys);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (uint32_t j = 0; j <=numAutoKeys; ++j){
    //     std::cout<<"gen.ModExp(j, 2 * N).ConvertToInt<LWEPlaintext>():"<<gen.ModExp(j, 2 * N).ConvertToInt<LWEPlaintext>()<<std::endl;
        for (uint32_t i = 0; i < n; ++i) {
            auto s{sv[i].ConvertToInt<int32_t>()};
            (*ek)[0][0][j*n+i] = KeyGenLMKCDEYC(params, skNTT, s > modHalf ? s - mod : s, j, gen.ModExp(j, 2 * N).ConvertToInt<LWEPlaintext>());
            // (*ek)[0][0][i] = KeyGenLMKCDEY(params, skNTT, s > modHalf ? s - mod : s);
            //如果s大于modHalf，则返回s - mod，否则返回s
        }
    }

#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (uint32_t j = 0; j <=numAutoKeys; ++j){
        for (uint32_t i = 0; i < n; ++i) {
            auto s{sv[i].ConvertToInt<int32_t>()};
            (*ek)[0][1][j*n+i] = KeyGenLMKCDEYC(params, skNTT, s > modHalf ? mod-s : -s, j, gen.ModExp(j, 2 * N).ConvertToInt<LWEPlaintext>());            
            // (*ek)[0][1][i] = KeyGenLMKCDEY(params, skNTT, s > modHalf ? mod - s : -s);
            //如果s大于modHalf，则返回s - mod，否则返回s
            // (*ek)[0][1][i] = KeyGenLMKCDEYC(params, skNTT, s > modHalf ? mod - s : -s, 0);
        }
    }    
    // NativeInteger gen = NativeInteger(5);

    (*ek)[0][2][0] = KeyGenAuto(params, skNTT, 2 * N - gen.ConvertToInt());
    //create and store the secret key square requred for creating RGSW ciphertext
    // (*ek)[0][2][0] = KeyGenSecSqr(params, skNTT);

    // m_window: window size, consider parameterization in the future
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(numAutoKeys))
    for (uint32_t i = 1; i <= numAutoKeys; ++i){
        (*ek)[0][2][i] = KeyGenAuto(params, skNTT, gen.ModExp(i, 2 * N).ConvertToInt<LWEPlaintext>());
    }
    return ek;
}

void RingGSWAccumulatorLMKCDEY::EvalAccTS(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek,
                                        RLWECiphertext& acc, const NativeVector& a) const {
    // assume a is all-odd ciphertext (using round-to-odd technique)
    size_t n             = a.GetLength();
    uint32_t Nh          = params->GetN() / 2;
    uint32_t M           = 2 * params->GetN();
    uint32_t numAutoKeys = params->GetNumAutoKeys();
    NativeInteger MNative(M);
    auto logGen = params->GetLogGen();

    std::vector<std::pair<NativeInteger,int>> b_vec;
    b_vec.reserve(n);
    
    for (std::size_t i = 0; i < n; ++i)
    {
        b_vec.emplace_back(std::pair<NativeInteger, int> (a[i],i));
    }    
    // NativeInteger temp;
    // std::cout<<"a vector is before sorting:"<<std::endl;
    // for (size_t i = 0; i < n; i++)
    // {   
    //     std::cout<<a[i]<<"\t";
    // }

    // std::cout<<"\n\n";
    // std::cout<<"b vector is:"<<std::endl;
    // for (NativeInteger i : b_vec)
    //     std::cout<< i << "\t";
    
    
    sort(b_vec.begin(), b_vec.end());
    // for (size_t i=0; i<b_vec.size();i++)
    // {
    //     std::cout<<b_vec[i].first<<"  "<<b_vec[i].second<<"\t";
    // }
    // std::cout<<"\n\n\n";
    // for (auto i: b_vec)
    //     std::cout<< i.first << " "<<i.second<<"\t"; //<<std::endl;

    // std::cout<<"\n\n\n";

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
    int index;
    NativeInteger c[n];
    // NativeInteger temp;
    for (size_t i=0;i<n; i++){
        index = b_vec[i].second;
        c[index] = b_vec[i].first;
    }

    std::unordered_map<int32_t, std::vector<int32_t>> permuteMap;
    
    for (size_t i = 0; i < n; i++) {  // put ail a_i in the permuteMap
        // make it odd; round-to-odd(https://eprint.iacr.org/2022/198) will improve error.
        int32_t aIOdd = NativeInteger(0).ModSubFast(a[i], MNative).ConvertToInt<uint32_t>() | 0x1;
        int32_t index = logGen[aIOdd];
        // std::cout<<" aIOdd and index are:"<<aIOdd<<"\t"<<index<<std::endl;
        if (permuteMap.find(index) == permuteMap.end()) {
            std::vector<int32_t> indexVec;
            permuteMap[index] = indexVec;
        }
        auto& indexVec = permuteMap[index];
        indexVec.push_back(i);
    }

    int count4=0, count5=0, count6=0;
    for (uint32_t i = Nh - 1; i > 0; i--) {
        int flag6=0;
        if (permuteMap.find(-i) != permuteMap.end()) {
            int flag=1;
            //Multiply the positive secret key
            auto& indexVec = permuteMap[-i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                if (flag)
                {
                    count4++;
                    flag--;    
                }                
            }        
            count6++;
            flag6=1;
        }
        if (permuteMap.find(i) != permuteMap.end()) {
            int flag=1;            
            //Multiply the positive secret key
            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                if (flag)
                {
                    count5++;
                    flag--;    
                }
            }       
            if (!flag6)
            {
                count6++;
            }
        }
    }
    int flag7=1;
    if (permuteMap.find(M) != permuteMap.end()) {
        int flag=1;
        auto& indexVec = permuteMap[M];
        for (size_t j = 0; j < indexVec.size(); j++) {
            if (flag)
            {
                count4++;
                flag--;    
            }
        }    
        if(flag7){
            count6++;    
            flag7--;
        }
    }
    if (permuteMap.find(0) != permuteMap.end()) {
        int flag=1;
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            if (flag)
            {
                count5++;
                flag--;    
            }            
        }        
        if(flag7)
            count6++;    
    } 
    

    
    // std::cout<<"\n";
    // std::cout<<"Total unique, negative and positive sets are equal to ------------------------>"<<count6<<"\t"<<count4<<"\t"<<count5<<std::endl;
    
    // std::cout<<"For positive and negative numbers:------------------------>"<<std::endl;
    /*
    for (uint32_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(i) != permuteMap.end()) {
            //Multiply the positive secret key
            auto& indexVec = permuteMap[i];
            // std::cout<<"indexVec+["<<i<<"]::";
            for (size_t j = 0; j < indexVec.size(); j++) {
                // std::cout<<indexVec[j]<<"\t";
                count5++;
            }        
            // std::cout<<"\n";
        }
        if (permuteMap.find(-i) != permuteMap.end()) {
            //Multiply the positive secret key
            auto& indexVec = permuteMap[i];
            // std::cout<<"indexVec-["<<i<<"]::";
            for (size_t j = 0; j < indexVec.size(); j++) {
                // std::cout<<indexVec[j]<<"\t";
                count5++;
            }        
            // std::cout<<"\n";
        }        
    }

    if (permuteMap.find(0) != permuteMap.end()) {
        auto& indexVec = permuteMap[0];
        // std::cout<<"indexVec["<<0<<"]::";
        for (size_t j = 0; j < indexVec.size(); j++) {
            // std::cout<<indexVec[j]<<"\t";
            count5++;
        }        
        // std::cout<<"\n";
    }    
    */

    // std::cout<<"\nValue of count5 is equal to----------------------------------------------->>>::"<<count5<<std::endl;
    int count=0;
    NativeInteger gen(5);
    uint32_t genInt       = 5;
    uint32_t nSkips       = 0;
    
    acc->GetElements()[1] = (acc->GetElements()[1]).AutomorphismTransform(M - genInt);
    
    // for a_j = -5^i
    for (uint32_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(-i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
            }
            // int flag=1;
            //Multiply the negative secret key
            auto& indexVec = permuteMap[-i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc); //check keys are correct
                // if (flag){
                //     std::cout<<"From negative:["<<i<<"]:"<<indexVec[j]<<"\t";
                //     flag=0;
                // }
                // else{
                //     std::cout<<indexVec[j]<<"\t";
                // }
            }
            // std::cout<<"\n";
        }
        if (permuteMap.find(i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
            }
            // int flag=1;
            //Multiply the positive secret key
            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccLMKCDEY(params, (*ek)[0][1][indexVec[j]], acc);
                // if (flag){
                //     std::cout<<"From positive:["<<i<<"]:"<<indexVec[j]<<"\t";
                //     flag=0;
                // }
                // else{
                //     std::cout<<indexVec[j]<<"\t";
                // }
            }
            // std::cout<<"\n";
        } 
        nSkips++;

        if (nSkips == numAutoKeys || i == 1) {
            // std::cout<<"I am from numAutoKeys || i == 1 and nSkips="<<i<<"\t"<<nSkips<<std::endl;
            Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
            nSkips = 0;
            count++;
        }
    }
    // std::cout<<"Value of count is:-------------------->"<<"\t"<<count<<std::endl;
    // if(count>1)
        // {for (int tk=0;tk<5;tk++) std::cout<<"\n\n";}
    // std::cout<<"nSkips is:"<<"\t"<<nSkips<<std::endl;
    // for -1
    if (permuteMap.find(M) != permuteMap.end()) {
        // int flag=1;
        auto& indexVec = permuteMap[M];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
            // if (flag){
            //     std::cout<<"From negative check ckeck:"<<indexVec[j]<<"\t";
            //     flag=0;
            // }
            // else{
            //     std::cout<<indexVec[j]<<"\t";
            // }
        }
        // std::cout<<"\n";
    }
    // for 1
    if (permuteMap.find(0) != permuteMap.end()) {
        // int flag=1;
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            // std::cout<<"I am from 0"<<std::endl;
            AddToAccLMKCDEY(params, (*ek)[0][1][indexVec[j]], acc);
            // if (flag){
            //     std::cout<<"From positive check ckeck:"<<indexVec[j]<<"\t";
            //     flag=0;
            // }
            // else{
            //     std::cout<<indexVec[j]<<"\t";
            // }
        }
        // std::cout<<"\n";
    }

    /*
    //perform keyswitch to (2*N-g)
    Automorphism(params, NativeInteger(M - genInt), (*ek)[0][1][0], acc);

    // for a_j = 5^i
    for (size_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][1][nSkips], acc);
                nSkips = 0;
            }

            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
            }
        }
        nSkips++;

        if (nSkips == numAutoKeys || i == 1) {
            Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][1][nSkips], acc);
            nSkips = 0;
        }
    }

    // for 0
    if (permuteMap.find(0) != permuteMap.end()) {
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
        }
    }
    */
}

void RingGSWAccumulatorLMKCDEY::EvalAccTSS(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek,
                                        RLWECiphertext& acc, const NativeVector& a) const {                                      
    // assume a is all-odd ciphertext (using round-to-odd technique)
    size_t n             = a.GetLength();
    uint32_t Nh          = params->GetN() / 2;
    uint32_t M           = 2 * params->GetN();
    uint32_t numAutoKeys = params->GetNumAutoKeys();
    NativeInteger MNative(M);
    auto logGen = params->GetLogGen();

    std::unordered_map<int32_t, std::vector<int32_t>> permuteMap;
    for (size_t i = 0; i < n; i++) {  // put ail a_i in the permuteMap
        // make it odd; round-to-odd(https://eprint.iacr.org/2022/198) will improve error.
        int32_t aIOdd1 = NativeInteger(0).ModSubFast(a[i], MNative).ConvertToInt<uint32_t>() | 0x1;
        int32_t index = logGen[aIOdd1];
        if (permuteMap.find(index) == permuteMap.end()) {
            std::vector<int32_t> indexVec;
            permuteMap[index] = indexVec;
        }
        auto& indexVec = permuteMap[index];
        indexVec.push_back(i);
    }

    NativeInteger gen(5);
    uint32_t genInt       = 5;
    uint32_t nSkips       = 0;
    uint32_t flag         = 1;
    acc->GetElements()[1] = (acc->GetElements()[1]).AutomorphismTransform(M - genInt);
    
    for (uint32_t i = Nh - 1; i > 0; i--) {

        if( (flag==1 && nSkips>0) && (permuteMap.find(-i) != permuteMap.end() || permuteMap.find(i) != permuteMap.end()) ) {
            for (size_t i = 0; i < nSkips; i++)
            {
                Automorphism(params, gen.ModExp(1, M), (*ek)[0][2][1], acc);                       
            }                
            flag=0,nSkips=0;
        }
        if (permuteMap.find(-i) != permuteMap.end()) {
            uint32_t Set_nextN=i-1, Total_jumpN;
            while(permuteMap.find(-Set_nextN) == permuteMap.end()) // && permuteMap.find(-Set_nextN) != permuteMap.find(-1))
            {
                if (Set_nextN==1){
                    Set_nextN--;             
                    break;                          
                }   
                if(Set_nextN==0){
                    break;                    
                }     
                Set_nextN--;                
            }
            Total_jumpN=(i-Set_nextN);

            uint32_t Set_nextP=i-1, Total_jumpP;
            while(permuteMap.find(Set_nextP) == permuteMap.end()) // && permuteMap.find(Set_nextN) != permuteMap.find(1))
            {
                if (Set_nextP==1){
                    Set_nextP--;   
                    break;    
                }
                if(Set_nextP==0){
                    break;                    
                }
                Set_nextP--;                
            }
            Total_jumpP=(i-Set_nextP);


            uint32_t Switch_key;
            if (Total_jumpP > Total_jumpN)
            {
                Switch_key=Total_jumpN;
            }
            else{      
                Switch_key=Total_jumpP;
            }           

            //Multiply the negative secret key
            auto& indexVec = permuteMap[-i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                if (j==indexVec.size()-1){
                    if(Switch_key>numAutoKeys){
                        AddToAccLMKCDEY_A(params, gen.ModExp(numAutoKeys, M), (*ek)[0][0][n*numAutoKeys + indexVec[j]], acc);                        
                        for (size_t i = 0; i < Switch_key-numAutoKeys; i++){
                            Automorphism(params, gen.ModExp(1, M), (*ek)[0][2][1], acc);  
                        }
                    }
                    else{                         
                        AddToAccLMKCDEY_A(params, gen.ModExp(Switch_key, M), (*ek)[0][0][n*Switch_key + indexVec[j]], acc);                      
                    }                                             
                }                
                else{
                    AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc); //check keys are correct                
                }
            }
            flag=0;
        }
        if (permuteMap.find(i) != permuteMap.end()) {        
            uint32_t Set_nextN=i-1, Total_jumpN;
            while(permuteMap.find(-Set_nextN) == permuteMap.end()) // && permuteMap.find(-Set_nextN) != permuteMap.find(-1))
            {
                if (Set_nextN==1){
                    Set_nextN--;         
                    break;        
                }     
                if (Set_nextN==0){
                    break;        
                }     
                Set_nextN--;                
            }
            Total_jumpN=(i-Set_nextN);

            uint32_t Set_nextP=i-1,Total_jumpP;
            while(permuteMap.find(Set_nextP) == permuteMap.end()) // && permuteMap.find(Set_nextN) != permuteMap.find(1))
            {
                if (Set_nextP==1){
                    Set_nextP--;
                    break;        
                }     
                if (Set_nextP==0){
                    break;        
                }             
                Set_nextP--;
            }
            
            Total_jumpP=(i-Set_nextP);
    
            uint32_t Switch_key;
            if (Total_jumpP > Total_jumpN)
            {
                Switch_key=Total_jumpN;
            }
            else{       
                Switch_key=Total_jumpP;
            }

            //Multiply the positive secret key
            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {                
                if(j==indexVec.size()-1){ 
                    if(Switch_key>numAutoKeys){                
                        AddToAccLMKCDEY_A(params, gen.ModExp(numAutoKeys, M), (*ek)[0][1][n*numAutoKeys + indexVec[j]], acc);
                        for (size_t i = 0; i <Switch_key-numAutoKeys; i++){
                            Automorphism(params, gen.ModExp(1, M), (*ek)[0][2][1], acc);                       
                        }
                    }
                    else{                                                 
                        AddToAccLMKCDEY_A(params, gen.ModExp(Switch_key, M), (*ek)[0][1][n*Switch_key + indexVec[j]], acc);                                                                                          
                    }                                                       
                    nSkips=0;
                }
                else{
                    AddToAccLMKCDEY(params, (*ek)[0][1][indexVec[j]], acc);
                }
            }
            flag=0;
        }       
        nSkips++;    
    }
    // for -1
    if (permuteMap.find(M) != permuteMap.end()) {
        auto& indexVec = permuteMap[M];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
        }
    }
    // for 1
    if (permuteMap.find(0) != permuteMap.end()) {
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccLMKCDEY(params, (*ek)[0][1][indexVec[j]], acc);
        }
    }
}

//Implementation of LEE et al. paper
void RingGSWAccumulatorLMKCDEY::EvalAcc(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek,
                                        RLWECiphertext& acc, const NativeVector& a) const {
    //please check keys before running this program...where are the aut keys are stored.                                            
    // assume a is all-odd ciphertext (using round-to-odd technique)
    size_t n             = a.GetLength();
    uint32_t Nh          = params->GetN() / 2;
    uint32_t M           = 2 * params->GetN();
    uint32_t numAutoKeys = params->GetNumAutoKeys();

    NativeInteger MNative(M);

    auto logGen = params->GetLogGen();
   std::unordered_map<int32_t, std::vector<int32_t>> permuteMap;

    for (size_t i = 0; i < n; i++) {  // put ail a_i in the permuteMap
        // make it odd; round-to-odd(https://eprint.iacr.org/2022/198) will improve error.
        int32_t aIOdd = NativeInteger(0).ModSubFast(a[i], MNative).ConvertToInt<uint32_t>() | 0x1;
        int32_t index = logGen[aIOdd];
        // std::cout<<" aIOdd and index are:"<<aIOdd<<"\t"<<index<<std::endl;
        if (permuteMap.find(index) == permuteMap.end()) {
            std::vector<int32_t> indexVec;
            permuteMap[index] = indexVec;
        }
        auto& indexVec = permuteMap[index];
        indexVec.push_back(i);
    }

    NativeInteger gen(5);
    uint32_t genInt       = 5;
    uint32_t nSkips       = 0;
    
    acc->GetElements()[1] = (acc->GetElements()[1]).AutomorphismTransform(M - genInt);
    
    // for a_j = -5^i
    for (uint32_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(-i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
            }
            auto& indexVec = permuteMap[-i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
            }
        }
        nSkips++;

        if (nSkips == numAutoKeys || i == 1) {
            Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
            nSkips = 0;
        }
    }

    // for -1
    if (permuteMap.find(M) != permuteMap.end()) {
        auto& indexVec = permuteMap[M];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
        }
    }

    //perform keyswitch to (2*N-g)
    Automorphism(params, NativeInteger(M - genInt), (*ek)[0][2][0], acc);

    // for a_j = 5^i
    for (size_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
            }

            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
            }
        }
        nSkips++;

        if (nSkips == numAutoKeys || i == 1) {
            Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
            nSkips = 0;
        }
    }

    // for 0
    if (permuteMap.find(0) != permuteMap.end()) {
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccLMKCDEY(params, (*ek)[0][0][indexVec[j]], acc);
        }
    }
}

// Encryption as described in Section 5 of https://eprint.iacr.org/2022/198
// Same as KeyGenAP, but only for X^{s_i}
// skNTT corresponds to the secret key z
RingGSWEvalKey RingGSWAccumulatorLMKCDEY::KeyGenLMKCDEY(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                        const NativePoly& skNTT, LWEPlaintext m) const {
    auto polyParams = params->GetPolyParams();//(Q,2N)
    auto Gpow       = params->GetGPower();

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);//确保dug的模数是Q
    // Reduce mod q (dealing with negative number as well)
    int64_t q  = params->Getq().ConvertToInt<int64_t>();//
    int64_t N  = params->GetN();
    int64_t mm = (((m % q) + q) % q) * (2 * N / q);
    bool isReducedMM{false};
    if (mm >= N) {
        mm -= N;
        isReducedMM = true;
    }

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    std::vector<NativePoly> tempA(digitsG2, NativePoly(dug, polyParams, Format::COEFFICIENT));
    RingGSWEvalKeyImpl result(digitsG2, 2);

    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i][0]=NativePoly(polyParams, Format::COEFFICIENT, true);
        result[i][1]=NativePoly(polyParams, Format::COEFFICIENT, true);
    }
    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i][0] = tempA[i];
        tempA[i].SetFormat(Format::EVALUATION);
        // result[i][1].SetFormat(Format::EVALUATION);    
        if (!isReducedMM){
            // result[i][i & 0x1][mm].ModAddFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G            
            if( i & 0x1){
                result[i][1][mm].ModAddFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G
            }
            else{
                result[i][1][mm].ModAddFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G
                result[i][1].SetFormat(Format::EVALUATION);  
                result[i][1] = result[i][1] * (-skNTT);
            }
        }
        else{
            // result[i][i & 0x1][mm].ModSubFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G
            if( i & 0x1){     
                result[i][1][mm].ModSubFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G
            }
            else{
                result[i][1][mm].ModSubFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G                
                result[i][1].SetFormat(Format::EVALUATION);                
                result[i][1] = result[i][1] * (-skNTT);
            }
        }
        result[i][1].SetFormat(Format::COEFFICIENT);                
        result[i][1] = result[i][1] + NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);                
        result[i][0].SetFormat(Format::EVALUATION);
        result[i][1].SetFormat(Format::EVALUATION);       
        result[i][1] += (tempA[i] *= skNTT);
    }

    return std::make_shared<RingGSWEvalKeyImpl>(result);
}

RingGSWEvalKey RingGSWAccumulatorLMKCDEY::KeyGenLMKCDEYC(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                        const NativePoly& skNTT, LWEPlaintext m, uint32_t j, LWEPlaintext t) const {
    auto polyParams = params->GetPolyParams();//(Q,2N)
    auto Gpow       = params->GetGPower();

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);//确保dug的模数是Q
    // Reduce mod q (dealing with negative number as well)
    int64_t q  = params->Getq().ConvertToInt<int64_t>();//
    int64_t N  = params->GetN();
    // int64_t mm = ((((m*t) % q) + q) % q) * (2 * N / q);
    int64_t mm = (((m % q) + q) % q) * (2 * N / q);    
    bool isReducedMM{false};
    if (mm >= N) {
        mm -= N;
        isReducedMM = true;
    }
    // std::cout<<"Here is the value of t:"<<t<<std::endl;
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    std::vector<NativePoly> tempA(digitsG2, NativePoly(dug, polyParams, Format::COEFFICIENT));
    RingGSWEvalKeyImpl result(digitsG2, 2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i][0]=NativePoly(polyParams, Format::COEFFICIENT, true);
        result[i][1]=NativePoly(polyParams, Format::COEFFICIENT, true);
    }
    for (uint32_t i = 0; i < digitsG2; ++i) {
        if (j==0){
            result[i][0] = tempA[i];
            tempA[i].SetFormat(Format::EVALUATION);            
            if (!isReducedMM){
                // result[i][i & 0x1][mm].ModAddFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G            
                if( i & 0x1){
                    result[i][1][mm].ModAddFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G
                }
                else{
                    result[i][1][mm].ModAddFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G
                    result[i][1].SetFormat(Format::EVALUATION);  
                    result[i][1] = result[i][1] * (-skNTT);
                }
            }
            else{
                // result[i][i & 0x1][mm].ModSubFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G
                if( i & 0x1){     
                    result[i][1][mm].ModSubFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G
                }
                else{
                    result[i][1][mm].ModSubFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G                
                    result[i][1].SetFormat(Format::EVALUATION);                
                    result[i][1] = result[i][1] * (-skNTT);
                }
            }
            result[i][1].SetFormat(Format::COEFFICIENT);                
            result[i][1] = result[i][1] + NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);                
            result[i][0].SetFormat(Format::EVALUATION);
            result[i][1].SetFormat(Format::EVALUATION);       
            result[i][1] += (tempA[i] *= skNTT);
        }
        else{
            auto skAuto{skNTT.AutomorphismTransform(t)};
            result[i][0] = tempA[i];
            tempA[i].SetFormat(Format::EVALUATION);            
            if (!isReducedMM){
                // result[i][i & 0x1][mm].ModAddFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G            
                if( i & 0x1){
                    result[i][1][mm].ModAddFastEq((Gpow[(i >> 1) + 1]),Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G
                    result[i][1]=result[i][1].AutomorphismTransform(t);
                }
                else{
                    result[i][1][mm].ModAddFastEq((Gpow[(i >> 1) + 1]), Q);  // (i even) Add G Multiple, (i odd) [a,as+e] + X^m*G
                    result[i][1]=result[i][1].AutomorphismTransform(t);                    
                    result[i][1].SetFormat(Format::EVALUATION);  
                    result[i][1] = result[i][1] * (-skAuto);
                }
            }
            else{
                // result[i][i & 0x1][mm].ModSubFastEq(Gpow[(i >> 1) + 1],Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G
                if( i & 0x1){     
                    result[i][1][mm].ModSubFastEq((Gpow[(i >> 1) + 1]),Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G
                    result[i][1]=result[i][1].AutomorphismTransform(t);                    
                }
                else{
                    result[i][1][mm].ModSubFastEq((Gpow[(i >> 1) + 1]),Q);  // (i even) Sub G Multiple, (i odd) [a,as+e] - X^m*G     
                    result[i][1]=result[i][1].AutomorphismTransform(t);                               
                    result[i][1].SetFormat(Format::EVALUATION);                
                    result[i][1] = result[i][1] * (-skAuto);
                }
            }
            result[i][1].SetFormat(Format::COEFFICIENT);                
            result[i][1] = result[i][1] + NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);                
            result[i][0].SetFormat(Format::EVALUATION);
            result[i][1].SetFormat(Format::EVALUATION);       
            result[i][1] += (tempA[i] *= skNTT);
        }
    }

    return std::make_shared<RingGSWEvalKeyImpl>(result);
}

// Generation of an autormorphism key
RingGSWEvalKey RingGSWAccumulatorLMKCDEY::KeyGenAuto(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                     const NativePoly& skNTT, LWEPlaintext k) const {
    auto polyParams{params->GetPolyParams()};
    // m_polyParams{std::make_shared<ILNativeParams>(2 * N, Q)},
    auto Gpow{params->GetGPower()};//m_Gpower,是一个3长度vector (0,1024,1048576)

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);
            
    auto skAuto{skNTT.AutomorphismTransform(k)};

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    RingGSWEvalKeyImpl result(digitsG, 2);

    for (uint32_t i = 0; i < digitsG; ++i) {
        result[i][0] = NativePoly(dug, polyParams, EVALUATION);
        result[i][1] = NativePoly(params->GetDgg(), polyParams, EVALUATION) - skAuto * Gpow[i + 1];
        result[i][1] += result[i][0] * skNTT;
    }
    return std::make_shared<RingGSWEvalKeyImpl>(result);
}

// Generation of encryption of seceret key square
RingGSWEvalKey RingGSWAccumulatorLMKCDEY::KeyGenSecSqr(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                     const NativePoly& skNTT) const {
    auto polyParams{params->GetPolyParams()};
    // m_polyParams{std::make_shared<ILNativeParams>(2 * N, Q)},
    auto Gpow{params->GetGPower()};//m_Gpower,是一个3长度vector (0,1024,1048576)

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);
    auto skNNTSq{skNTT * skNTT};
    // std::cout<<"RLWE secret key is:\n"<<skNTT;
    // std::cout<<"RLWE secret key is:\n"<<skNNTSq;

    // approximate gadget decomposition is used; the first digit is ignored
    // uint32_t digitsG{params->GetDigitsG() - 1};
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    RingGSWEvalKeyImpl result(digitsG2, 2);

    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i][0] = NativePoly(dug, polyParams, EVALUATION);
        result[i][1] = NativePoly(params->GetDgg(), polyParams, EVALUATION) - skNNTSq * Gpow[i + 1];
        result[i][1] += result[i][0] * skNTT;
    }
    return std::make_shared<RingGSWEvalKeyImpl>(result);
}

// Create RGSW ciphertext of the given RLWE' ciphertext
RingGSWACCKey RingGSWAccumulatorLMKCDEY::RLWEtoRGSW(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                     ConstRingGSWEvalKey& SkeySq, ConstRingGSWEvalKey& key) const {    
    auto polyParams{params->GetPolyParams()};
    RingGSWACCKey Com_Key = std::make_shared<RingGSWACCKeyImpl>(1, 1, 2);
    (*Com_Key)[0][0][1] = std::make_shared<RingGSWEvalKeyImpl>(*key);
    auto Gpow{params->GetGPower()};//m_Gpower,是一个3长度vector (0,1024,1048576)
    
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1}; 
    RingGSWEvalKeyImpl result(digitsG2, 2);
    
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG2))    
    for (uint32_t k = 0; k < digitsG2; ++k){
        std::vector<NativePoly> ct(key->GetElements()[k]);
        ct[0].SetFormat(Format::COEFFICIENT);
        ct[1].SetFormat(Format::COEFFICIENT);
    
        // approximate gadget decomposition is used; the first digit is ignored
        uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
        std::vector<NativePoly> dct(digitsG2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));
        SignedDigitDecompose(params, ct, dct);

    // calls digitsG2 NTTs
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG2))
        for (uint32_t d = 0; d < digitsG2; ++d)
            dct[d].SetFormat(Format::EVALUATION);

        // acc = dct * ek (matrix product);
        const std::vector<std::vector<NativePoly>>& ev = SkeySq->GetElements();
        result[k][0]                          = (dct[0] * ev[0][0]);
        for (uint32_t d = 1; d < digitsG2; ++d)
            result[k][0] += (dct[d] * ev[d][0]);
        //Add (b,0) to the ciphertext 
        result[k][0] += ct[1];

        result[k][1] = (dct[0] *= ev[0][1]);
        for (uint32_t d = 1; d < digitsG2; ++d)
            result[k][1] += (dct[d] *= ev[d][1]);
    }

    (*Com_Key)[0][0][0]=std::make_shared<RingGSWEvalKeyImpl>(result);
    return Com_Key;
}

//RLWE' and RGSW key multplication-----------------------------------------------------------------------------------------------
RingGSWEvalKey RingGSWAccumulatorLMKCDEY::MulRLWEtoRGSW(const std::shared_ptr<RingGSWCryptoParams>& params,
                                        ConstRingGSWEvalKey& rgswkey0, ConstRingGSWEvalKey& rgswkey1, ConstRingGSWEvalKey& ek) const {

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
    RingGSWEvalKeyImpl result(digitsG2, 2);

    for (uint32_t k = 0; k < digitsG2; ++k){
        // const std::vector<NativePoly> ev = ek->GetElements()[k];
        std::vector<NativePoly> ct(ek->GetElements()[k]);
        ct[0].SetFormat(Format::COEFFICIENT);
        ct[1].SetFormat(Format::COEFFICIENT);
        
        // approximate gadget decomposition is used; the first digit is ignored
        uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};
        std::vector<NativePoly> dct(digitsG2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));

        SignedDigitDecompose(params, ct, dct);

        // calls digitsG2 NTTs
        #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG2))
            for (uint32_t d = 0; d < digitsG2; ++d)
                dct[d].SetFormat(Format::EVALUATION);

            // acc = dct * ek (matrix product);
            const std::vector<std::vector<NativePoly>>& ev0 = rgswkey0->GetElements();
            const std::vector<std::vector<NativePoly>>& ev1 = rgswkey1->GetElements();

            result[k][0]                          = (dct[0] * ev0[0][0]);
            for (uint32_t d = 1; d < digitsG2; ++d)
                result[k][0] += (dct[d] * ev0[d][0]);

            result[k][1] = (dct[0] *= ev1[0][1]);
            for (uint32_t d = 1; d < digitsG2; ++d)
                result[k][1] += (dct[d] *= ev1[d][1]);
    }

    // std::cout<<"Here is the value of result[0][0]:"<<result[0][0]<<std::endl;
    // std::cout<<"Here is the value of result[0][1]:"<<result[0][1]<<std::endl;
    return std::make_shared<RingGSWEvalKeyImpl>(result);
}

//RLWE' and RGSW key multplication-----------------------------------------------------------------------------------------------


// LMKCDEY Accumulation as described in https://eprint.iacr.org/2022/198
// Same as AP, but multiplied once
void RingGSWAccumulatorLMKCDEY::AddToAccLMKCDEY(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                ConstRingGSWEvalKey& ek, RLWECiphertext& acc) const {
    std::vector<NativePoly> ct(acc->GetElements());
    ct[0].SetFormat(Format::COEFFICIENT);
    ct[1].SetFormat(Format::COEFFICIENT);
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};

    std::vector<NativePoly> dct(digitsG2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));

    SignedDigitDecompose(params, ct, dct);

    // calls digitsG2 NTTs
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG2))
    for (uint32_t d = 0; d < digitsG2; ++d)
        dct[d].SetFormat(Format::EVALUATION);

    // acc = dct * ek (matrix product);
    const std::vector<std::vector<NativePoly>>& ev = ek->GetElements();
    acc->GetElements()[0]                          = (dct[0] * ev[0][0]);
    for (uint32_t d = 1; d < digitsG2; ++d)
        acc->GetElements()[0] += (dct[d] * ev[d][0]);
    acc->GetElements()[1] = (dct[0] *= ev[0][1]);
    for (uint32_t d = 1; d < digitsG2; ++d)
        acc->GetElements()[1] += (dct[d] *= ev[d][1]);
}

void RingGSWAccumulatorLMKCDEY::AddToAccLMKCDEY_A(const std::shared_ptr<RingGSWCryptoParams>& params,
                                const NativeInteger& a, ConstRingGSWEvalKey& ek, RLWECiphertext& acc) const {

    // precompute bit reversal for the automorphism into vec
    uint32_t N{params->GetN()};
    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, a.ConvertToInt<usint>(), &vec);//

    acc->GetElements()[1] = acc->GetElements()[1].AutomorphismTransform(a.ConvertToInt<usint>(), vec);
    acc->GetElements()[0] = acc->GetElements()[0].AutomorphismTransform(a.ConvertToInt<usint>(), vec);

    std::vector<NativePoly> ct(acc->GetElements());
    ct[0].SetFormat(Format::COEFFICIENT);
    ct[1].SetFormat(Format::COEFFICIENT);
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1) << 1};

    std::vector<NativePoly> dct(digitsG2, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));

    SignedDigitDecompose(params, ct, dct);

    // calls digitsG2 NTTs
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG2))
    for (uint32_t d = 0; d < digitsG2; ++d)
        dct[d].SetFormat(Format::EVALUATION);

    // acc = dct * ek (matrix product);
    const std::vector<std::vector<NativePoly>>& ev = ek->GetElements();
    acc->GetElements()[0]                          = (dct[0] * ev[0][0]);
    for (uint32_t d = 1; d < digitsG2; ++d)
        acc->GetElements()[0] += (dct[d] * ev[d][0]);
    acc->GetElements()[1] = (dct[0] *= ev[0][1]);
    for (uint32_t d = 1; d < digitsG2; ++d)
        acc->GetElements()[1] += (dct[d] *= ev[d][1]);
}


// Automorphism
void RingGSWAccumulatorLMKCDEY::Automorphism(const std::shared_ptr<RingGSWCryptoParams>& params, const NativeInteger& a,
                                             ConstRingGSWEvalKey& ak, RLWECiphertext& acc) const {
    // precompute bit reversal for the automorphism into vec
    uint32_t N{params->GetN()};
    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, a.ConvertToInt<usint>(), &vec);//

    acc->GetElements()[1] = acc->GetElements()[1].AutomorphismTransform(a.ConvertToInt<usint>(), vec);

    NativePoly cta(acc->GetElements()[0]);
    acc->GetElements()[0].SetValuesToZero();
    cta = cta.AutomorphismTransform(a.ConvertToInt<usint>(), vec);
    
    cta.SetFormat(COEFFICIENT);

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    std::vector<NativePoly> dcta(digitsG, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));

    SignedDigitDecompose(params, cta, dcta);

#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dcta[d].SetFormat(Format::EVALUATION);

    // acc = dct * input (matrix product);
    const std::vector<std::vector<NativePoly>>& ev = ak->GetElements();
    for (uint32_t d = 0; d < digitsG; ++d)
        acc->GetElements()[0] += (dcta[d] * ev[d][0]);
    for (uint32_t d = 0; d < digitsG; ++d)
        acc->GetElements()[1] += (dcta[d] *= ev[d][1]);
}

};  // namespace lbcrypto
