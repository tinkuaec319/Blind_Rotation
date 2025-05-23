//==================================================================================
// BSD 2-Clause License
//
// Copyright (c) 2014-2023, NJIT, Duality Technologies Inc. and other contributors
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

#include "vntru-acc-xzddf.h"

#include <string>
#include <vector>
namespace lbcrypto {

VectorNTRUACCKey VectorNTRUAccumulatorXZDDF::KeyGenAcc(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                       const NativePoly& skNTT, const NativePoly& invskNTT,
                                                       ConstLWEPrivateKey& LWEsk) const {
    auto sv{LWEsk->GetElement()};
    auto mod{sv.GetModulus().ConvertToInt<int32_t>()};  //q_ks
    auto modHalf{mod >> 1};
    size_t n{sv.GetLength()};
    auto q{params->Getq().ConvertToInt<size_t>()};
    VectorNTRUACCKey ek = std::make_shared<VectorNTRUACCKeyImpl>(1, 2, q - 1 > n + 1 ? q - 1 : n + 1);
    //生成评估秘钥
    auto s{sv[0].ConvertToInt<int32_t>()};                                          // 0 +-1
    (*ek)[0][0][0] = KDMKeyGenXZDDF(params, invskNTT, s > modHalf ? mod - s : -s);  //第一个evk(KDM-form)

#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))
    for (size_t i = 1; i < n; ++i) {
        auto s{sv[i].ConvertToInt<int32_t>()};
        (*ek)[0][0][i] = KeyGenXZDDF(params, invskNTT, s > modHalf ? mod - s : -s);
        //如果s大于modHalf，则返回s - mod，否则返回s
    }
    auto sums = 0;
    for (size_t i = 0; i < n; ++i) {
        auto s{sv[i].ConvertToInt<int32_t>()};
        sums = sums +s;
    }
    sums %= mod;
    if (sums > modHalf) {
        sums -= mod;
    }
    (*ek)[0][0][n] = KeyGenXZDDF(params, invskNTT, sums);

    //生成自同构秘钥
    int64_t intq = params->Getq().ConvertToInt<int64_t>();  
    int64_t N    = params->GetN();
    for (auto i = 0; i < intq - 1; ++i) {
        (*ek)[0][1][i] = KeyGenAuto(params, skNTT, invskNTT, (2 * N / intq) * (i + 1) + 1);
    }
    return ek;
}

VectorNTRUACCKey VectorNTRUAccumulatorXZDDF::KeyGenAccS(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                       const NativePoly& skNTT, const NativePoly& invskNTT,
                                                       ConstLWEPrivateKey& LWEsk) const {
    auto sv{LWEsk->GetElement()};
    auto mod{sv.GetModulus().ConvertToInt<int32_t>()};  //q_ks
    auto modHalf{mod >> 1};
    int64_t N    = params->GetN();
    size_t n{sv.GetLength()};
    int64_t numAutoKeys{params->GetNumAutoKeys()};
    uint32_t total_rgsw_keys=(numAutoKeys+1)*n;
    NativeInteger gen = NativeInteger(5);

    // auto q{params->Getq().ConvertToInt<size_t>()};

    // auto qp{params->Getqp()};
    // uint32_t qp{params->Getqp()}; //.ConvertToInt<size_t>()};
    // dim2, 0: for RGSW(X^si), 1: for automorphism keys
    // only w automorphism keys required
    // allocates (n - w) more memory for pointer (not critical for performance)
    VectorNTRUACCKey ek = std::make_shared<VectorNTRUACCKeyImpl>(1, 3, total_rgsw_keys+1); //(1, 2, q - 1 > n + 1 ? q - 1 : n + 1);
    
    //生成评估秘钥
    //create the first ev_key with all secret elements addition and use it to conver ciphertext to NTRU ciphertext 
    // std::cout<<"mod value is:"<<q<<std::endl;
    // auto sums = 0;
    // for (size_t i = 0; i < n; ++i) {
    //     auto s{sv[i].ConvertToInt<int32_t>()};
    //     // if (s>modHalf)
    //     //     s=s-mod;
    //     // std::cout<<"Secret key is:"<<s<<std::endl;
    //     // sums = sums + s;
    // }
    // std::cout<<"Secret key sum is:"<<sums<<std::endl;
    // sums = (-sums)% mod;
    // std::cout<<"After negation and mod secret key sum is:"<<sums<<std::endl;
    // if (sums > modHalf) {
    //     sums -= mod;
    // }
    // (*ek)[0][0][0] = KDMKeyGenXZDDFS(params, invskNTT, s > modHalf ? s-mod : s);

    // (*ek)[0][0][2*n] = KDMKeyGenXZDDFS(params, invskNTT, sums);

    /*This is should be enable when window style bootstrapping is performed without merging with that of the 
     *bootstrapping key.
    */
    // (*ek)[0][0][n] = KDMKeyGenXZDDFS(params, invskNTT, 0);

    //This is should be enable when the secret key is merged with that of the automorphism key
    (*ek)[0][0][total_rgsw_keys] = KDMKeyGenXZDDFS(params, invskNTT, 0);


    // auto s{sv[0].ConvertToInt<int32_t>()};                                          // 0 +-1
    // (*ek)[0][0][0] = KDMKeyGenXZDDF(params, invskNTT, s > modHalf ? mod - s : -s);  //第一个evk(KDM-form)

    //compute vector ntru ciphertext for all secret elements of the LWE secret
// #pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(2*n))
    // std::cout<<"Value of n is equal to :"<<n<<std::endl;        

    // compute vector ntru ciphertext for all secret elements of the LWE secret     
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n)) 
    for (uint32_t j = 0; j <=numAutoKeys; ++j){       
        for (size_t i = 0; i<n; ++i) {
            auto s{sv[i].ConvertToInt<int32_t>()};   
            // (*ek)[0][0][i] = KeyGenXZDDFS(params, invskNTT,   s > modHalf ? s-mod : s);
            (*ek)[0][0][j*n+i] = KeyGenXZDDFSC(params, skNTT, invskNTT, s > modHalf ? s - mod : s, j, gen.ModExp(j, 2 * N).ConvertToInt<LWEPlaintext>());            
        }
    }

    // compute vector ntru ciphertext for all -secret elements of the LWE secret
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(n))        
    for (uint32_t j = 0; j <=numAutoKeys; ++j){
        for (size_t i = 0; i < n; ++i) {
                auto s{sv[i].ConvertToInt<int32_t>()};
                // std::cout<<"Here is the secret components:(i,s)"<<i<<" , "<<s<<std::endl;        
                // (*ek)[0][1][i] = KeyGenXZDDFS(params, invskNTT,   s > modHalf ? mod-s : -s);
                (*ek)[0][1][j*n+i] = KeyGenXZDDFSC(params, skNTT, invskNTT, s > modHalf ? mod-s : -s, j, gen.ModExp(j, 2 * N).ConvertToInt<LWEPlaintext>());                            
        }
    }


    (*ek)[0][2][0] = KeyGenAutoS(params, skNTT, invskNTT, 2 * N - gen.ConvertToInt());

    // m_window: window size, consider parameterization in the future
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(numAutoKeys))    
    for (int64_t i = 1; i <=numAutoKeys ; ++i) {
        (*ek)[0][2][i] = KeyGenAutoS(params, skNTT, invskNTT, gen.ModExp(i, 2 * N).ConvertToInt()); //invskNTT, (2 * N / intq) * (i + 1) + 1);
    }
    // (*ek)[0][1][numAutoKeys+1] = KeyGenAutoS(params, skNTT, invskNTT, sums);
    // std::cout<<"Generating keys"<<std::endl;
    return ek;
}

//This is for keeping q number of automorphism keys.
void VectorNTRUAccumulatorXZDDF::EvalAcc(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a) const {
    size_t n   = a.GetLength();
    uint32_t N = params->GetN();
    int32_t q  = params->Getq().ConvertToInt<int32_t>();
    // int32_t qp = params->Getqp();
    // int32_t qp = params->Getqp(); //.ConvertToInt<int32_t>();
    std::vector<uint32_t> ua(n);
    std::vector<uint32_t> w(n);
    std::vector<uint32_t> invw(n + 1);
    invw[n] = 1;
    std::vector<NativeInteger> NATIVEw(n);  //自同构的次数
    std::vector<uint32_t> invindex(n);      //对应到autk 的index
    // std::cout<<"Test 5"<<std::endl;
    for (size_t i = 0; i < n; i++) {
        ua[i]   = a[i].ConvertToInt<int32_t>();       //a
        w[i]    = (2 * N / q) * ua[i] + 1;            //w_i
        invw[i] = ModInverse(w[i], 2 * N) % (2 * N);  //w_inv
    }
    
    for (size_t i = 0; i < n; i++) {
        NATIVEw[i] = NativeVector::Integer((w[i] * invw[i + 1]) % (2 * N));
        invindex[i] = (NATIVEw[i].ConvertToInt<int32_t>() - 2*N/q -1) / (2*N/q);
        // invindex[i] = (NATIVEw[i].ConvertToInt<int32_t>() - 2*N/q -1) / (2*N/q);
        // std::cout<<"invindex[i] is:"<<invindex[i]<<std::endl;
    }
    // std::cout<<"Test 7"<<std::endl;
    for (size_t i = 0; i < n; i++) {
        AddToAccXZDDF(params, (*ek)[0][0][i], acc);  ///evk_{0 ~ n-1}
        if (NATIVEw[i].ConvertToInt<int32_t>() != 1) {
            Automorphism(params, NATIVEw[i], (*ek)[0][1][invindex[i]], acc);
        }
    }
//    std::cout<<"Test 6"<<std::endl;
    AddToAccXZDDF(params, (*ek)[0][0][n], acc);
}

//This is for the window style automorphism computation.
void VectorNTRUAccumulatorXZDDF::EvalAccS(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a) const {
                                                                        
    size_t   n = a.GetLength();
    uint32_t N = params->GetN();
    uint32_t Nh= N / 2;
    uint32_t M = 2 * N;
    uint32_t numAutoKeys = params->GetNumAutoKeys();

    NativeInteger MNative(M);
    // int32_t q  = params->Getq().ConvertToInt<int32_t>();

    int count=0,count1=0,count2=0;

    auto logGen = params->GetLogGen();
    std::unordered_map<int32_t, std::vector<int32_t>> permuteMap;

    //convert all coefficients of ciphertext first components into odd terms to make them elements of Z*_2N 
    //so that they becomes isomorphic to Z*_2N =Z_N/2 * Z_2 
    // std::cout<<"Value of c0 is:"<<a<<std::endl;
    for (size_t i = 0; i < n; i++) {
        // aiOdd = NativeInteger(0).ModSubFast(a[i]*(2*N/q)+ 1, MNative).ConvertToInt<int32_t>();            //w_i        
        // std::cout<<"Value of a[i] and aiOdd is:"<<a[i]<<"\t"<<aiOdd<<std::endl;
        // 2n/q*s[i] is encrypted in the secret key encryption.
        int32_t aiOdd = NativeInteger(0).ModSubFast(a[i], MNative).ConvertToInt<uint32_t>() | 0x1;
        // aiOdd =a[i].ConvertToInt<uint32_t>() | 0x1;

        // put all aiOdd in the permuteMap
        int32_t index = logGen[aiOdd];
        if (permuteMap.find(index) == permuteMap.end()) {
            std::vector<int32_t> indexVec;
            permuteMap[index] = indexVec;
            count++;
        }
        auto& indexVec = permuteMap[index];
        indexVec.push_back(i);
    }
    
    NativeInteger gen(5);
    uint32_t genInt       = 5;
    uint32_t nSkips       = 0;


    acc->GetElements() = (acc->GetElements()).AutomorphismTransform(M - genInt);
    
    /*acc is just a polynomial and need is to convert it into an NTRU polynomial and is converted
    to NTRU polynomial by multilying it with encrypted 1/f  */
    AddToAccXZDDFS(params, (*ek)[0][0][n], acc);
    

    // for a_j = -5^i
    for (uint32_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(-i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
                count1++;
            }
            auto& indexVec = permuteMap[-i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}
            }
        }
        nSkips++;

        if (nSkips == numAutoKeys || i == 1) {
            Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
            nSkips = 0;
            count1++;
        }
    }
    // for -1
    if (permuteMap.find(M) != permuteMap.end()) {
        auto& indexVec = permuteMap[M];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}            
        }
    }

    Automorphism(params, NativeInteger(M - genInt), (*ek)[0][2][0], acc);
    // for a_j = 5^i
    for (size_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
                count2++;
            }

            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}
            }
        }
        nSkips++;

        if (nSkips == numAutoKeys || i == 1) {
            Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
            nSkips = 0;
            count2++;
        }
    }
    // for 0
    if (permuteMap.find(0) != permuteMap.end()) {
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}            
        }
    } 

    std::cout<<"(count,count1,count2) = ("<<count<<","<<count1<<","<<count2<<")"<<std::endl;
}

//This is to combine sets corresponding to same number and reduce the total number of automorphisms.
void VectorNTRUAccumulatorXZDDF::EvalAccTS(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a) const {
                                                                        
    size_t   n = a.GetLength();
    uint32_t N = params->GetN();
    uint32_t Nh= N / 2;
    uint32_t M = 2 * N;
    uint32_t numAutoKeys = params->GetNumAutoKeys();
    NativeInteger MNative(M);
    // int32_t q  = params->Getq().ConvertToInt<int32_t>();

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
        // temp=b_vec[i].first;
        c[index] = b_vec[i].first;
    }

    auto logGen = params->GetLogGen();
    std::unordered_map<int32_t, std::vector<int32_t>> permuteMap;

    //convert all coefficients of ciphertext first components into odd terms to make them elements of Z*_2N 
    //so that they becomes isomorphic to Z*_2N =Z_N/2 * Z_2 
    // std::cout<<"Value of c0 is:"<<a<<std::endl;
    for (size_t i = 0; i < n; i++) {
        // aiOdd = NativeInteger(0).ModSubFast(a[i]*(2*N/q)+ 1, MNative).ConvertToInt<int32_t>();            //w_i        
        // std::cout<<"Value of a[i] and aiOdd is:"<<a[i]<<"\t"<<aiOdd<<std::endl;
        // 2n/q*s[i] is encrypted in the secret key encryption.
        int32_t aiOdd = NativeInteger(0).ModSubFast(a[i], MNative).ConvertToInt<uint32_t>() | 0x1;
        // aiOdd =a[i].ConvertToInt<uint32_t>() | 0x1;

        // put all aiOdd in the permuteMap
        int32_t index = logGen[aiOdd];
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


    acc->GetElements() = (acc->GetElements()).AutomorphismTransform(M - genInt);
    
    /*acc is just a polynomial and need is to convert it into an NTRU polynomial and is converted
    to NTRU polynomial by multilying it with encrypted 1/f  */
    AddToAccXZDDFS(params, (*ek)[0][0][n], acc);
    

    // for a_j = -5^i
    for (uint32_t i = Nh - 1; i > 0; i--) {
        if (permuteMap.find(-i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
            }
            auto& indexVec = permuteMap[-i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}
            }
        }
        if (permuteMap.find(i) != permuteMap.end()) {
            if (nSkips != 0) {  // Rotation by 5^nSkips
                Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][2][nSkips], acc);
                nSkips = 0;
            }

            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                AddToAccXZDDFS(params, (*ek)[0][1][indexVec[j]], acc);  ///evk_{0 ~ n-1}
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
            AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}            
        }
    }

    // for 0
    if (permuteMap.find(0) != permuteMap.end()) {
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccXZDDFS(params, (*ek)[0][1][indexVec[j]], acc);  ///evk_{0 ~ n-1}            
        }
    } 

    // Automorphism(params, NativeInteger(M - genInt), (*ek)[0][1][0], acc);
    // // for a_j = 5^i
    // for (size_t i = Nh - 1; i > 0; i--) {


    //     if (nSkips == numAutoKeys || i == 1) {
    //         Automorphism(params, gen.ModExp(nSkips, M), (*ek)[0][1][nSkips], acc);
    //         nSkips = 0;
    //         count2++;
    //     }
    // }


    // std::cout<<"(count,count1,count2) = ("<<count<<","<<count1<<","<<count2<<")"<<std::endl;
}

//This is to combine sets corresponding to same number and reduce the total number of automorphisms.
void VectorNTRUAccumulatorXZDDF::EvalAccTSC(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                         ConstVectorNTRUACCKey& ek, NTRUCiphertext& acc, const NativeVector& a) const {
                                                                        
    size_t   n = a.GetLength();
    uint32_t N = params->GetN();
    uint32_t Nh= N / 2;
    uint32_t M = 2 * N;
    uint32_t numAutoKeys = params->GetNumAutoKeys();
    NativeInteger MNative(M);
    uint32_t total_rgsw_keys=(numAutoKeys+1)*n;

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
        // temp=b_vec[i].first;
        c[index] = b_vec[i].first;
    }

    auto logGen = params->GetLogGen();
    std::unordered_map<int32_t, std::vector<int32_t>> permuteMap;

    //convert all coefficients of ciphertext first components into odd terms to make them elements of Z*_2N 
    //so that they becomes isomorphic to Z*_2N =Z_N/2 * Z_2 
    // std::cout<<"Value of c0 is:"<<a<<std::endl;
    for (size_t i = 0; i < n; i++) {
        // aiOdd = NativeInteger(0).ModSubFast(a[i]*(2*N/q)+ 1, MNative).ConvertToInt<int32_t>();            //w_i        
        // std::cout<<"Value of a[i] and aiOdd is:"<<a[i]<<"\t"<<aiOdd<<std::endl;
        // 2n/q*s[i] is encrypted in the secret key encryption.
        int32_t aiOdd = NativeInteger(0).ModSubFast(a[i], MNative).ConvertToInt<uint32_t>() | 0x1;
        // aiOdd =a[i].ConvertToInt<uint32_t>() | 0x1;

        // put all aiOdd in the permuteMap
        int32_t index = logGen[aiOdd];
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
    acc->GetElements() = (acc->GetElements()).AutomorphismTransform(M - genInt);
    
    /*acc is just a polynomial and need is to convert it into an NTRU polynomial and is converted
    to NTRU polynomial by multilying it with encrypted 1/f  */
    AddToAccXZDDFS(params, (*ek)[0][0][total_rgsw_keys], acc);
    

    // for a_j = -5^i
    for (uint32_t i = Nh - 1; i > 0; i--) {

        if( (flag==1 && nSkips>0) && (permuteMap.find(-i) != permuteMap.end() || permuteMap.find(i) != permuteMap.end()) ) {
            for (size_t i = 0; i < nSkips; i++)
            {
                Automorphism(params, gen.ModExp(1, M), (*ek)[0][2][1], acc);                       
            }                     
            flag=0,nSkips=0;
        }

        if (permuteMap.find(-i) != permuteMap.end()) {
            // std::cout<<"------ Value of i is:-------"<<i<<std::endl;
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
            // std::cout<<"Total_jump positive is:"<<Total_jumpP<<std::endl;


            uint32_t Switch_key;
            if (Total_jumpP > Total_jumpN)
            {
                // std::cout<<"We need to key switch to negative:"<<Total_jumpN<<std::endl;
                Switch_key=Total_jumpN;
            }
            else{
                // std::cout<<"We need to key switch to positive:"<<Total_jumpP<<std::endl;                
                Switch_key=Total_jumpP;
            }
            //Multiply the negative secret key
            auto& indexVec = permuteMap[-i];
            for (size_t j = 0; j < indexVec.size(); j++) {
                if (j==indexVec.size()-1){
                    if(Switch_key>numAutoKeys){
                        // std::cout<<"I entered in negative:"<<Switch_key<<"\t"<<Switch_key-numAutoKeys<<std::endl;                        
                        // AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}    
                        // Automorphism(params, gen.ModExp(numAutoKeys, M), (*ek)[0][2][numAutoKeys], acc);  
                        AddToAccXZDDFS_A(params, gen.ModExp(numAutoKeys, M), (*ek)[0][0][n*numAutoKeys + indexVec[j]], acc);                        
                        for (size_t i = 0; i < Switch_key-numAutoKeys; i++){
                            Automorphism(params, gen.ModExp(1, M), (*ek)[0][2][1], acc);  
                        }
                    }
                    else{                         
                        AddToAccXZDDFS_A(params, gen.ModExp(Switch_key, M), (*ek)[0][0][n*Switch_key + indexVec[j]], acc);
                        // AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}    
                        // Automorphism(params, gen.ModExp(Switch_key, M), (*ek)[0][2][Switch_key], acc);                          
                    }
                    // std::cout<<Switch_key<<"\t";                                                  
                }                
                else{
                    AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}    
                }
            }                
            flag=0;
        }
        if (permuteMap.find(i) != permuteMap.end()) {
            // std::cout<<"++++++ Value of i is:++++++"<<i<<std::endl;            
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
                // std::cout<<"We need to key switch to negative:"<<Total_jumpN<<std::endl;
                Switch_key=Total_jumpN;
            }
            else{
                // std::cout<<"We need to key switch to positive:"<<Total_jumpP<<std::endl;                
                Switch_key=Total_jumpP;
            }  
            //Multiply the positive secret key
            auto& indexVec = permuteMap[i];
            for (size_t j = 0; j < indexVec.size(); j++) {          
                if(j==indexVec.size()-1){ 
                    if(Switch_key>numAutoKeys){
                        // std::cout<<"I entered in positive:"<<Switch_key<<"\t"<<Switch_key-numAutoKeys<<std::endl;                        
                        // AddToAccXZDDFS(params,  (*ek)[0][1][indexVec[j]], acc);  ///evk_{0 ~ n-1}                            
                        // Automorphism(params, gen.ModExp(numAutoKeys, M), (*ek)[0][2][numAutoKeys], acc);      
                        AddToAccXZDDFS_A(params, gen.ModExp(numAutoKeys, M), (*ek)[0][1][n*numAutoKeys + indexVec[j]], acc);                        
                        for (size_t i = 0; i <Switch_key-numAutoKeys; i++){
                            Automorphism(params, gen.ModExp(1, M), (*ek)[0][2][1], acc);                       
                        }
                    }
                    else{                                                 
                        AddToAccXZDDFS_A(params, gen.ModExp(Switch_key, M), (*ek)[0][1][n*Switch_key + indexVec[j]], acc);                        
                        // AddToAccXZDDFS(params,  (*ek)[0][1][indexVec[j]], acc);  ///evk_{0 ~ n-1}      
                        // Automorphism(params, gen.ModExp(Switch_key, M), (*ek)[0][2][Switch_key], acc);                                                                      
                    }
                    // std::cout<<"nSkips is:"<<nSkips<<std::endl;       
                    // std::cout<<Switch_key<<"\t";                                                                 
                    nSkips=0;
                }
                else{
                    AddToAccXZDDFS(params,  (*ek)[0][1][indexVec[j]], acc);  ///evk_{0 ~ n-1}                    
                }

                // if (flag){
                //     std::cout<<"From positive:["<<i<<"]:"<<indexVec[j]<<"\t";
                //     flag=0;
                // }
                // else{
                //     std::cout<<indexVec[j]<<"\t";
                // }
            }
            flag=0;
        }         
        nSkips++;                   
    }

    // for -1
    if (permuteMap.find(M) != permuteMap.end()) {
        auto& indexVec = permuteMap[M];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccXZDDFS(params, (*ek)[0][0][indexVec[j]], acc);  ///evk_{0 ~ n-1}            
        }
    }

    // for 0
    if (permuteMap.find(0) != permuteMap.end()) {
        auto& indexVec = permuteMap[0];
        for (size_t j = 0; j < indexVec.size(); j++) {
            AddToAccXZDDFS(params, (*ek)[0][1][indexVec[j]], acc);  ///evk_{0 ~ n-1}            
        }
    } 
}

void VectorNTRUAccumulatorXZDDF::AddToAccXZDDFS(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                               ConstVectorNTRUEvalKey& ek, NTRUCiphertext& acc) const {
    NativePoly ct(acc->GetElements());
    ct.SetFormat(Format::COEFFICIENT);
    
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{(params->GetDigitsG() - 1)};

    std::vector<NativePoly> dct(digitsG, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));  // d-1维N长多项式
    SignedDigitDecompose(params, ct, dct);                                                        //分解acc
    // calls digitsG2 NTTs
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);

#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);

    // acc = dct * ek (matrix product);
    const std::vector<NativePoly>& ev = ek->GetElements();
    for (uint32_t d = 0; d < digitsG; ++d){
        sum += (dct[d] *= ev[d]);
    }
    acc->GetElements() = sum;
}

void VectorNTRUAccumulatorXZDDF::AddToAccXZDDFS_A(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                            const NativeInteger& a, ConstVectorNTRUEvalKey& ek, NTRUCiphertext& acc) const {
    // precompute bit reversal for the automorphism into vec
    uint32_t N{params->GetN()};
    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, a.ConvertToInt<usint>(), &vec);  //

    NativePoly ct(acc->GetElements());
    acc->GetElements().SetValuesToZero();
    ct = ct.AutomorphismTransform(a.ConvertToInt<usint>(), vec);
    ct.SetFormat(COEFFICIENT);

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    std::vector<NativePoly> dct(digitsG, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));
    SignedDigitDecompose(params, ct, dct);

    // calls digitsG2 NTTs
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);

    // acc = dct * ek (matrix product);
    const std::vector<NativePoly>& ev = ek->GetElements();
    for (uint32_t d = 0; d < digitsG; ++d){
        sum += (dct[d] *= ev[d]);
    }
    acc->GetElements() = sum;
}

// KDM-form
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KDMKeyGenXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                             const NativePoly& invskNTT, LWEPlaintext m) const {
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);  //确保dug的模数是Q
                        //Reduce mod q (dealing with negative number as well)
    int64_t N  = params->GetN();
    int64_t mm = (((m % N) + N) % N);  // 0 1 N-1
    // int64_t mm = (((m % q) + q) % q) * (2*N/q);  // 0 1 N-1
    bool isReducedMM{false};
    if (m < 0) {
        isReducedMM = true;
    }

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1)};
    std::vector<NativePoly> tempA(digitsG2, NativePoly(dug, polyParams, Format::COEFFICIENT));
    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        if (!isReducedMM)
            result[i][mm].ModAddFastEq(Gpow[i + 1],Q);  // g+X^m*G
        else
            result[i][mm].ModSubFastEq(Gpow[i + 1],Q);  // g-X^m*G
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

// KDM-form
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KDMKeyGenXZDDFS(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                             const NativePoly& invskNTT, LWEPlaintext m) const {
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //
    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);  //确保dug的模数是Q
    
    //Reduce mod q (dealing with negative number as well)
    int64_t N  = params->GetN();
    int64_t mm = (((m % N) + N) % N);  // 0 1 N-1
    bool isReducedMM{false};
    if (m < 0) {
        isReducedMM = true;
    }

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1)};
    std::vector<NativePoly> tempA(digitsG2, NativePoly(dug, polyParams, Format::COEFFICIENT));
    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        if (!isReducedMM){
            result[i][mm].ModAddFastEq(Gpow[i + 1],Q);  // g+X^m*G
        }
        else{
            result[i][mm].ModSubFastEq(Gpow[i + 1],Q);  // g-X^m*G
        }
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;      // (g+X^m*G)*(1/f)   or  (g-X^m*G)*(1/f)
        // std::cout<<"Encrypted (g+X^m*G)*(1/f) or (g-X^m*G)*(1/f) is:"<<result[i]<<std::endl;
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}


//NO KDM-form
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KeyGenXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                          const NativePoly& invskNTT, LWEPlaintext m) const {
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //
    NativeInteger Q{params->GetQ()};
    int64_t N  = params->GetN();
    int64_t mm = (((m % N) + N) % N);  // 0 1 q-1
    bool isReducedMM{false};
    if (m < 0) {
        isReducedMM = true;
    }
    uint32_t digitsG2{(params->GetDigitsG() - 1)};  //2
    NativePoly zeroPoly(polyParams, Format::COEFFICIENT);
    zeroPoly.SetValuesToZero();
    std::vector<NativePoly> tempA(digitsG2, zeroPoly);

    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        // result[i][0] = tempA[i];
        tempA[i].SetFormat(Format::COEFFICIENT);
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;  // g/f
        if (!isReducedMM)
            tempA[i][mm].ModAddFastEq(Gpow[i + 1], Q);  // X^m*G
        else
            tempA[i][mm].ModSubFastEq(Gpow[i + 1], Q);  // X^m*G
        tempA[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] + tempA[i];               // g/f+X^m*G   / g/f-X^m*G
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

//NO KDM-form
//bootstrapping key generation
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KeyGenXZDDFS(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                          const NativePoly& invskNTT, LWEPlaintext m) const {

    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //     
    NativeInteger Q{params->GetQ()};
    int64_t N  = params->GetN();
    int64_t q  = params->Getq().ConvertToInt<int64_t>();//
    int64_t mm = (((m % N) + N) % N) * (2*N/q);  // 0 1 q-1
    bool isReducedMM{false};
    if (mm >= N) {
        mm -=N;
        isReducedMM = true;
    }

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1)};  //2
    NativePoly zeroPoly(polyParams, Format::COEFFICIENT);
    zeroPoly.SetValuesToZero();
    std::vector<NativePoly> tempA(digitsG2, zeroPoly);

    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        // result[i][0] = tempA[i];
        tempA[i].SetFormat(Format::COEFFICIENT);
        result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
        result[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] * invskNTT;  // g/f
        if (!isReducedMM)
            tempA[i][mm].ModAddFastEq(Gpow[i + 1], Q);  //  X^m*G
        else
            tempA[i][mm].ModSubFastEq(Gpow[i + 1], Q);  // -X^m*G
        tempA[i].SetFormat(Format::EVALUATION);
        result[i] = result[i] + tempA[i];               //g/f+X^m*G  and/or  g/f-X^m*G
    }
    // std::cout<<"Complted at KeyGenXZDDFS"<<std::endl;    
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

//bootstrapping key generation
VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KeyGenXZDDFSC(const std::shared_ptr<VectorNTRUCryptoParams>& params, const NativePoly& skNTT,
                                                          const NativePoly& invskNTT, LWEPlaintext m, uint32_t j, LWEPlaintext t) const {

    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //     
    NativeInteger Q{params->GetQ()};
    int64_t N  = params->GetN();
    int64_t q  = params->Getq().ConvertToInt<int64_t>();//
    int64_t mm = (((m % N) + N) % N) * (2*N/q);  // 0 1 q-1
    bool isReducedMM{false};
    if (mm >= N) {
        mm -=N;
        isReducedMM = true;
    }

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG2{(params->GetDigitsG() - 1)};  //2
    NativePoly zeroPoly(polyParams, Format::COEFFICIENT);
    zeroPoly.SetValuesToZero();
    std::vector<NativePoly> tempA(digitsG2, zeroPoly);

    VectorNTRUEvalKeyImpl result(digitsG2);
    for (uint32_t i = 0; i < digitsG2; ++i) {
        if(j==0){    
            tempA[i].SetFormat(Format::COEFFICIENT);
            result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
            result[i].SetFormat(Format::EVALUATION);
            result[i] = result[i] * invskNTT;  // g/f
            if (!isReducedMM)
                tempA[i][mm].ModAddFastEq(Gpow[i + 1], Q);  //  X^m*G
            else
                tempA[i][mm].ModSubFastEq(Gpow[i + 1], Q);  // -X^m*G
            tempA[i].SetFormat(Format::EVALUATION);
            result[i] = result[i] + tempA[i];               //g/f+X^m*G  and/or  g/f-X^m*G
        }
        else{
            auto skNTTAuto{skNTT.AutomorphismTransform(t)};
            tempA[i].SetFormat(Format::COEFFICIENT);
            result[i] = NativePoly(params->GetDgg(), polyParams, Format::COEFFICIENT);  //采样g
            result[i].SetFormat(Format::EVALUATION);
            if (!isReducedMM){
                tempA[i][mm].ModAddFastEq(Gpow[i + 1], Q);  //  X^m*G             
                tempA[i]=tempA[i].AutomorphismTransform(t);
                tempA[i].SetFormat(Format::EVALUATION);
                tempA[i]=tempA[i] * skNTTAuto;                
            }
            else{
                tempA[i][mm].ModSubFastEq(Gpow[i + 1], Q);  // -X^m*G         
                tempA[i]=tempA[i].AutomorphismTransform(t);
                tempA[i].SetFormat(Format::EVALUATION);
                tempA[i]=tempA[i] * skNTTAuto;                
            }
            result[i] = result[i] + tempA[i];               //g/f+X^m*G*f^t  and/or  g/f-X^m*G*f^t                        
            result[i] = result[i] * invskNTT;  // g/f^t             //this needs to be verified            
            // tempA[i].SetFormat(Format::EVALUATION);
            

        }
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KeyGenAuto(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                         const NativePoly& skNTT, const NativePoly& invskNTT,
                                                         LWEPlaintext k) const {
    //auto polyParams{params->GetPolyParams()};
    // m_polyParams{std::make_shared<ILNativeParams>(2 * N, Q)},
    // auto Gpow{params->GetGPower()};//m_Gpower,是一个3长度vector (0,1024,1048576)
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);
    auto skAuto{skNTT.AutomorphismTransform(k)};  //生成f(X^k)

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    VectorNTRUEvalKeyImpl result(digitsG);

    for (uint32_t i = 0; i < digitsG; ++i) {
        result[i] = NativePoly(params->GetDgg(), polyParams, EVALUATION) + skAuto * Gpow[i + 1];
        result[i] = result[i] * invskNTT;
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}

VectorNTRUEvalKey VectorNTRUAccumulatorXZDDF::KeyGenAutoS(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                                         const NativePoly& skNTT, const NativePoly& invskNTT,
                                                         LWEPlaintext k) const {
    //auto polyParams{params->GetPolyParams()};
    // m_polyParams{std::make_shared<ILNativeParams>(2 * N, Q)},
    // auto Gpow{params->GetGPower()};//m_Gpower,是一个3长度vector (0,1024,1048576)
    auto polyParams = params->GetPolyParams();  //(Q,2N)
    auto Gpow       = params->GetGPower();      //

    DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeInteger Q{params->GetQ()};
    dug.SetModulus(Q);
    auto skAuto{skNTT.AutomorphismTransform(k)};  //生成f(X^k)

    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    VectorNTRUEvalKeyImpl result(digitsG);

    for (uint32_t i = 0; i < digitsG; ++i) {
        //in origianl NTRU encryption it was negative.
        // result[i] = NativePoly(params->GetDgg(), polyParams, EVALUATION) + skAuto * Gpow[i + 1]; 
        result[i] = NativePoly(params->GetDgg(), polyParams, EVALUATION) + skAuto * Gpow[i + 1];  // g+s^t*G
        result[i] = result[i] * invskNTT;  //(g+s^t*G)/f
    }
    return std::make_shared<VectorNTRUEvalKeyImpl>(result);
}


void VectorNTRUAccumulatorXZDDF::AddToAccXZDDF(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                               ConstVectorNTRUEvalKey& ek, NTRUCiphertext& acc) const {
    NativePoly ct(acc->GetElements());
    ct.SetFormat(Format::COEFFICIENT);
    
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{(params->GetDigitsG() - 1)};
    std::vector<NativePoly> dct(digitsG,
                                NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));  // d-1维N长多项式
    SignedDigitDecompose(params, ct, dct);     
                                                       //分解acc
    // calls digitsG2 NTTs
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);
    // acc = dct * ek (matrix product);
    const std::vector<NativePoly>& ev = ek->GetElements();
    for (uint32_t d = 0; d < digitsG; ++d)
        sum += (dct[d] *= ev[d]);

    acc->GetElements() = sum;
}


void VectorNTRUAccumulatorXZDDF::Automorphism(const std::shared_ptr<VectorNTRUCryptoParams>& params,
                                              const NativeInteger& a, ConstVectorNTRUEvalKey& ak,
                                              NTRUCiphertext& acc) const {
    // precompute bit reversal for the automorphism into vec
    uint32_t N{params->GetN()};
    std::vector<usint> vec(N);
    PrecomputeAutoMap(N, a.ConvertToInt<usint>(), &vec);  //

    NativePoly ct(acc->GetElements());
    acc->GetElements().SetValuesToZero();
    ct = ct.AutomorphismTransform(a.ConvertToInt<usint>(), vec);
    ct.SetFormat(COEFFICIENT);
    
    // approximate gadget decomposition is used; the first digit is ignored
    uint32_t digitsG{params->GetDigitsG() - 1};
    std::vector<NativePoly> dct(digitsG, NativePoly(params->GetPolyParams(), Format::COEFFICIENT, true));
    SignedDigitDecompose(params, ct, dct);
    
    NativePoly sum(params->GetPolyParams(), Format::EVALUATION, true);
#pragma omp parallel for num_threads(OpenFHEParallelControls.GetThreadLimit(digitsG))
    for (uint32_t d = 0; d < digitsG; ++d)
        dct[d].SetFormat(Format::EVALUATION);

    // acc = dct * input (matrix product);
    const std::vector<NativePoly>& ev = ak->GetElements();
    for (uint32_t d = 0; d < digitsG; ++d)
        sum += (dct[d] * ev[d]);

    acc->GetElements() = sum;
}


};  // namespace lbcrypto