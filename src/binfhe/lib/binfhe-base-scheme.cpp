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


#include "binfhe-base-scheme.h"
#include <iostream>  // Include the input/output stream library
#include <fstream>
#include <string>
#include <typeinfo>
namespace lbcrypto {

// wrapper for KeyGen methods
RingGSWBTKey BinFHEScheme::KeyGen(const std::shared_ptr<BinFHECryptoParams>& params, ConstLWEPrivateKey& LWEsk,
                                  KEYGEN_MODE keygenMode = SYM_ENCRYPT) const {
    const auto& LWEParams = params->GetLWEParams();

    RingGSWBTKey ek;  
    LWEPrivateKey skN;  
    if (keygenMode == SYM_ENCRYPT) {
        skN = LWEscheme->KeyGen(LWEParams->GetN(), LWEParams->GetQ());
    }
    else if (keygenMode == PUB_ENCRYPT) {
        ConstLWEKeyPair kpN = LWEscheme->KeyGenPair(LWEParams);
        skN                 = kpN->secretKey;
        ek.Pkey             = kpN->publicKey;
    }
    else {
        OPENFHE_THROW(config_error, "Invalid KeyGen mode");
    }

    uint32_t N{LWEParams->GetN()};
    //save the created secret key
    std::ofstream outputFile("LWE_Secret_key1.txt");  // Open/create a file named "test.txt" for writing
    if (outputFile.is_open()) {  // Check if the file was successfully opened
        // Write some text into the file
        for (uint32_t i=0; i<N; i++)
        {
            outputFile<<skN->GetElement()[i] <<std::endl;
        }
    }
    // Close the file
    outputFile.close();  // Close the file after writing



    /*-----------KSkey-------------*/
    ek.KSkey = LWEscheme->KeySwitchGen(LWEParams, LWEsk, skN);
    const auto& RGSWParams = params->GetRingGSWParams();
    const auto& polyParams = RGSWParams->GetPolyParams(); 
    NativePoly skNPoly(polyParams);
    skNPoly.SetValues(skN->GetElement(), Format::COEFFICIENT);
    skNPoly.SetFormat(Format::EVALUATION);
    /*-----------BSkey-------------*/
    ek.BSkey = ACCscheme->KeyGenAcc(RGSWParams, skNPoly, LWEsk);

    return ek;
}

VectorNTRUBTKey BinFHEScheme::NKeyGen(const std::shared_ptr<BinFHECryptoParams>& params, ConstLWEPrivateKey& LWEsk,
                                      KEYGEN_MODE keygenMode = SYM_ENCRYPT) const {
    const auto& LWEParams = params->GetLWEParams();
    VectorNTRUBTKey ek; 

    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    uint32_t Q{VNTRUParams->GetQ().ConvertToInt<uint32_t>()};
    uint32_t N{VNTRUParams->GetN()};

    NativeVector NatVec(N, Q);
    NativeVector NatVec_inv(N, Q);
    Get_invertible_NativeVector(NatVec, NatVec_inv, Q, N);
    LWEPrivateKey LWEskN = std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(NatVec));

    ek.KSkey = LWEscheme->KeySwitchGen(LWEParams, LWEsk, LWEskN);

    NativePoly skNPoly(polyParams);
    skNPoly.SetValues(NatVec, Format::COEFFICIENT);
    NativePoly invskNPoly(polyParams);
    invskNPoly.SetValues(NatVec_inv, Format::COEFFICIENT);
    skNPoly.SetFormat(Format::EVALUATION);
    invskNPoly.SetFormat(Format::EVALUATION);
    
    //write NTRU secret key into a file
    std::ofstream outputFile("Self_created_NTRU_Secret_key.txt");  // Open/create a file named "test.txt" for writing
    if (outputFile.is_open()) {  // Check if the file was successfully opened
        // Write some text into the file
        for (uint32_t i=0; i<N; i++)
        {
            outputFile << NatVec[i] <<std::endl;
        }
    }
    // Close the file
    outputFile.close();  // Close the file after writing

    ek.BSkey = NACCscheme->KeyGenAcc(VNTRUParams, skNPoly, invskNPoly, LWEsk);
    return ek;
}
VectorNTRUBTKey BinFHEScheme::NKeyGenS(const std::shared_ptr<BinFHECryptoParams>& params, ConstLWEPrivateKey& LWEsk,
                                      KEYGEN_MODE keygenMode = SYM_ENCRYPT) const {

    const auto& LWEParams = params->GetLWEParams();
    VectorNTRUBTKey ek; 

    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    uint32_t Q{VNTRUParams->GetQ().ConvertToInt<uint32_t>()};
    uint32_t N{VNTRUParams->GetN()};

    NativeVector NatVec(N, Q);
    NativeVector NatVec_inv(N, Q);
    Get_invertible_NativeVector(NatVec, NatVec_inv, Q, N);
    LWEPrivateKey LWEskN = std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(NatVec));
    
    /*-----------KSkey-------------*/
    ek.KSkey = LWEscheme->KeySwitchGenS(LWEParams, LWEsk, LWEskN);

    NativePoly skNPoly(polyParams);
    skNPoly.SetValues(NatVec, Format::COEFFICIENT);
    NativePoly invskNPoly(polyParams);
    invskNPoly.SetValues(NatVec_inv, Format::COEFFICIENT);
    skNPoly.SetFormat(Format::EVALUATION);
    invskNPoly.SetFormat(Format::EVALUATION);
    
    //write NTRU secret key into a file
    std::ofstream outputFile("Self_created_NTRU_Secret_key.txt");  // Open/create a file named "test.txt" for writing
    if (outputFile.is_open()) {  // Check if the file was successfully opened
        // Write some text into the file
        for (uint32_t i=0; i<N; i++)
        {
            outputFile << NatVec[i] <<std::endl;
        }
    }
    // Close the file
    outputFile.close();  // Close the file after writing
    /*-----------BSkey-------------*/
    ek.BSkey = NACCscheme->KeyGenAccS(VNTRUParams, skNPoly, invskNPoly, LWEsk);
    return ek;
}

void Get_invertible_NativeVector(NativeVector& NatVec, NativeVector& NatVec_inv, uint32_t q_boot, uint32_t N) {
    // 三值集合的均匀分布
    // uniform_int_distribution<int> ternary_sampler(-1,1);
    //正态分布
    normal_distribution<double> gaussian_sampler(0.0, 1);
    // 随机引擎
    default_random_engine rand_engine(std::chrono::system_clock::now().time_since_epoch().count());

    std::vector<int> vec     = std::vector<int>(N, 0);
    std::vector<int> vec_inv = std::vector<int>(N, 0);

    uint32_t half_q_boot = q_boot / 2;
    //polynomial with the coefficient vector vec (will be generated later)
    ZZ_pX poly;
    //element of Z_(q_boot)
    ZZ_p coef;
    coef.init(ZZ(q_boot));
    //the inverse of poly modulo poly_mod (will be generated later)
    ZZ_pX inv_poly;
    //random sampling
    //int sum = 0;
    while (true) {
        //create the polynomial with the coefficient vector of the desired form
        SetCoeff(poly, 0, gaussian_sampler(rand_engine));
        for (uint32_t i = 1; i < N; i++) {
            coef = gaussian_sampler(rand_engine);
            // if(coef == 0)
            // {
            //     sum++;
            // }
            SetCoeff(poly, i, coef);
        }
        //cout<<double(sum)/N<<endl;
        //test invertibility
        try {
            // static ZZ_pX get_def_poly()
            ZZ_pX def_poly;
            ZZ_p coef;
            coef.init(ZZ(q_boot));
            coef = 1;
            SetCoeff(def_poly, 0, coef);
            SetCoeff(def_poly, N, coef);

            InvMod(inv_poly, poly, def_poly);
            break;
        }
        catch (...) {
            //cout << "Polynomial " << poly << " isn't a unit" << endl;
            continue;
        }
    }
    uint32_t tmp_coef;
    for (uint32_t i = 0; i <= deg(poly); i++) {
        tmp_coef = conv<long>(poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec[i] = tmp_coef;
    }

    for (uint32_t i = 0; i <= deg(inv_poly); i++) {
        tmp_coef = conv<long>(inv_poly[i]);
        if (tmp_coef > half_q_boot)
            tmp_coef -= q_boot;
        vec_inv[i] = tmp_coef;
    }
    // vector<int> to NativePoly
    for (uint32_t i = 0; i < N; i++) {
        int32_t v     = vec[i];
        int32_t v_inv = vec_inv[i];
        if (v < 0){
            NatVec[i] = q_boot - typename NativeVector::Integer(-v);
        }
        else
            NatVec[i] = typename NativeVector::Integer(v);
        if (v_inv < 0)
            NatVec_inv[i] = q_boot - typename NativeVector::Integer(-v_inv);
        else
            NatVec_inv[i] = typename NativeVector::Integer(v_inv);
    }
}

LWECiphertext BinFHEScheme::EvalBinGate(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const VectorNTRUBTKey& EK, ConstLWECiphertext& ct1,
                                        ConstLWECiphertext& ct2) const {

    // std::cout<<"Reached to EvalBinGate"<<std::endl;
    if (ct1 == ct2)
        OPENFHE_THROW(config_error, "Input ciphertexts should be independant");

    // By default, we compute XOR/XNOR using a combination of AND, OR, and NOT gates 迭代计算
    if ((gate == XOR) || (gate == XNOR)) {
        const auto& ctAND1 = EvalBinGate(params, AND, EK, ct1, EvalNOT(params, ct2));
        const auto& ctAND2 = EvalBinGate(params, AND, EK, EvalNOT(params, ct1), ct2);
        const auto& ctOR   = EvalBinGate(params, OR, EK, ctAND1, ctAND2);

        // NOT is free so there is not cost to do it an extra time for XNOR
        return (gate == XOR) ? ctOR : EvalNOT(params, ctOR);
    }

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ct1);
    //构造(5q/8,0)
    auto n = params->GetLWEParams()->Getn();
    uint32_t N{params->GetLWEParams()->GetN()};

    NativeVector zero(n,0);
    uint32_t q  = params->GetLWEParams()->Getq().ConvertToInt<uint32_t>();
    uint32_t qp = params->GetLWEParams()->Getqp(); //.ConvertToInt<uint32_t>();
    // uint32_t qpmodKS = params->GetLWEParams()->GetqpKS(); //.ConvertToInt<uint32_t>();
    zero.SetModulus(q);

    //For LWE or lwr 
    NativeInteger temp_b= 5*q/(8*qp);  // 5*p/(8) This had been changed due to new cipher mod p
    // std::cout<<"temp_b is="<<temp_b.Mod(q/qp)<<std::endl;
    // LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q/qp)));
    LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q/qp)));

    // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
    // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
    if ((gate == XOR_FAST) || (gate == XNOR_FAST)) {
        LWEscheme->EvalSubEq(qp, ctprep, ct2);
        LWEscheme->EvalAddEq(qp, ctprep, ctprep);
    }
    else {
        LWEscheme->EvalAddEq(qp, ctprep, ct2);
        LWEscheme->EvalSubEq(qp, ct_temp, ctprep);
    }

    //This is for normal Bootstrapping using automorphisms
    auto acc{BootstrapGateCoreS(params, gate, EK.BSkey, ct_temp)};  //一个NativePoly，多项式
    
    // This is to check the computed output ..................................................................................................................
    /*
    NativePoly ct(acc->GetElements());
    ct.SetFormat(Format::COEFFICIENT);
    
    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    NativePoly sum(polyParams, Format::EVALUATION, true);

    uint32_t Q1{VNTRUParams->GetQ().ConvertToInt<uint32_t>()};
    NativeVector NtruSecO(N,0);
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<N; i++)
    {
        NtruSecO[i]=0;
    }
    ifstream inputFile1("Self_created_NTRU_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFile1.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line_1;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int tt_t=0;
    int temp045=0;
    while (getline(inputFile1, line_1)) {
        temp045=stoi(line_1) % Q1;
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        NtruSecO[tt_t]=temp045;
        tt_t++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile1.close();
    

    // NativeVector NatVec(N, Q);
    // NativeVector NatVec_inv(N, Q);
    // Get_invertible_NativeVector(NatVec, NatVec_inv, Q, N);
    // LWEPrivateKey LWEskN = std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(NatVec));

    // ek.KSkey = LWEscheme->KeySwitchGen(LWEParams, LWEsk, LWEskN);

    NativePoly skNPoly(polyParams, Format::COEFFICIENT);
    skNPoly.SetValuesToZero();
    for (uint32_t i=0; i<N; i++){
        skNPoly[i]= NtruSecO[i];
        // std::cout<<"skNPoly[i]"<<zeroPoly[i]<<std::endl;
    }
    // NativePoly skNPoly(polyParams);
    // std::cout<<"working fine"<<std::endl;
    // skNPoly.SetValues(zeroPoly, Format::COEFFICIENT);
    // std::cout<<"working fine"<<std::endl;
    // NativePoly invskNPoly(polyParams);
    // invskNPoly.SetValues(NatVec_inv, Format::COEFFICIENT);
    skNPoly.SetFormat(Format::EVALUATION);
    ct.SetFormat(Format::EVALUATION);
    // invskNPoly.SetFormat(Format::EVALUATION);
    // std::cout<<"working fine"<<std::endl;
    sum=skNPoly*ct;
    sum.SetFormat(Format::COEFFICIENT);

    // std::cout<<"Content of the accumulator polynomial after NTRU decryption is:..........................."<<((sum[0]+(Q1 >> 3) + 1)*4)/Q1<<std::endl;
    // std::cout<<"Content of the accumulator polynomial after NTRU decryption is:..........................."<<(sum[0]*4)/Q1<<std::endl;
    std::cout<<"\nHere is the content of the NTRU decrypted ciphertext:"<<std::endl;
    for (uint32_t i=0; i<N; i++){
        std::cout<<sum[i]<<" , ";
    }
    std::cout<<" \n ";
    */
    // This is to check the computed output ends here ..................................................................................................................





    NativePoly& accVec{acc->GetElements()};
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    // std::cout<<"Content of accumulator vector is:"<<std::endl;
    // for (uint32_t i=0;i<N;i++)
    // {
    //     std::cout<<accVec[i]<<",";
    // }
    // std::cout<<"Before transpose accmulator content is:"<<std::endl;
    // for (u_int32_t i=0;i<N;i++)
    //     std::cout<<(acc->GetElements()[i]*4)/Q<<" , ";
    // std::cout<<"\n";
    accVec = accVec.Transpose();

    // std::cout<<"After transpose"<<std::endl;
    // for (u_int32_t i=0;i<N;i++)
    //     std::cout<<"accVec["<<i<<"]"<<acc->GetElements()[i]<<std::endl;
    // std::cout<<"\n\n\nContent of accumulator vector transpose is:"<<std::endl;
    // for (uint32_t i=0;i<N;i++)
    // {
    //     std::cout<<accVec[i]<<",";
    // }
    accVec.SetFormat(Format::COEFFICIENT);
    // std::cout<<"\n\n\nContent of accumulator vector coefficient is:"<<std::endl;
    // for (uint32_t i=0;i<N;i++)
    // {
    //     std::cout<<accVec[i]<<",";
    // }
    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    // const auto& LWEParams = params->GetLWEParams();
    // NativeInteger Q{LWEParams->GetQ()};

    // uint32_t dim{LWEParams->Getn()};
    NativeInteger b{(Q >> 3) + 1};
    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec.GetValues()), std::move(b));
    
    //Check the returned ciphertext is correct ...........................................................................................
    // NativePoly ct(acc->GetElements());
    // ct.SetFormat(Format::COEFFICIENT);
    
    // const auto& VNTRUParams = params->GetVectorNTRUParams();
    // const auto& polyParams  = VNTRUParams->GetPolyParams();  
    // NativePoly sum(polyParams, Format::EVALUATION, true);

    // uint32_t Q1{VNTRUParams->GetQ().ConvertToInt<uint32_t>()};

    NativeVector NtruSec(N,0);
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<N; i++)
    {
        NtruSec[i]=0;
    }
    ifstream inputFile("Self_created_NTRU_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFile.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    // Read each line of the file and print it 
    string line;
    int ttt=0;
    int temp45=0;
    while (getline(inputFile, line)) {
        temp45=stoll(line) % Q.ConvertToInt();
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        NtruSec[ttt]=temp45;
        ttt++;
    }
    inputFile.close();
    NativeInteger temp{0};
     for (uint32_t i=0; i<N; i++){
        temp +=ctExt->GetA()[i]*NtruSec[i];   
    }
    temp=(temp+Q)%Q;
    NativeInteger bb=ctExt->GetB();
    bb=((bb-temp)+Q)%Q;
    std::cout<<"Decrypted result of NTRU/LWE ciphertext  is:"<<Q<<"\t\t"<<bb<<"\t\t"<<(((bb+Q/8)*4)/Q) % 4<<std::endl;
    //Checking completed ..................... ...........................................................................................

    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    
    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    //-----------------------------------------------------------------------I am going to decrypt here........................................................
        // NativeInteger n{ctKS->GetLWEParams->Getn()};
    NativeInteger Q3{ctKS->GetModulus()};
    NativeVector LweSec(n,0);
    ifstream inputFilel("LWE_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFilel.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line1;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int tttt=0;
    int temp455=0;
    while (getline(inputFilel, line1)) {
        // std::cout<<"line1"<<line1<<std::endl;
        temp455=stoll(line1);
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        LweSec[tttt]=temp455 % Q3.ConvertToInt();
        tttt++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFilel.close();
    NativeInteger temp2=0;
     for (uint32_t i=0; i<n; i++){
        temp2 +=ctKS->GetA()[i]*(LweSec[i]%Q3);
        // std::cout<<"skNPoly[i]"<<zeroPoly[i]<<std::endl;
    }
    if(qp>1){
        temp2=NativeInteger(static_cast<BasicInteger>( std::floor(0.5 + temp2.ConvertToDouble() * (Q3/qp).ConvertToDouble() / Q3.ConvertToDouble())));
        // RoundqQ( temp2, Q3/qp, Q3);
    }
    else{
        temp2=(temp2+Q3)%Q3;
    }
    
    NativeInteger bb2=ctKS->GetB();
    bb2=((bb2-temp2)+(Q3/qp))%(Q3/qp);
    // std::cout<<"Before ModSwitch Decrypted result of LWE cipher is:"<<Q3<<"\t\t"<<bb2<<"\t\t"<<((bb2+Q3/(qp*8))*4)/Q3<<std::endl;
    //-----------------------------------------------------------------------Decryption done I am going to decrypt here........................................................

    //Modswitch 
    if(qp>1){
        return LWEscheme->ModSwitchLwr(qp, ct1->GetModulus(), ctKS);
    }
    else{
        return LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
    }

    





    // std::cout<<"(Q, LWEParams->GetqKS(), ct1->GetModulus(), qp) are:("<<Q<<" , "<<LWEParams->GetqKS()<<" , "<<ct1->GetModulus()<<" , "<<qp<< " ) "<<std::endl;
    //define secret key
    //..........................................This code is to decrypt NTRU ciphertext.......................................
    // NativeInteger Q_ks{params->GetLWEParams()->GetqKS()};
    
    /*
    uint32_t N{params->GetLWEParams()->GetN()};
    // int Q_ks = Q_ks.ConvertToInt();
    // LWEParams->GetqKS().ConvertToInt();

    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    int Q_ks = params->GetLWEParams()->GetqKS().ConvertToInt();
    
    int Array1[N];
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<N; i++)
    {
        Array1[i]=0;
    }
    ifstream inputFile("NTRU_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFile.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int ttt=0;
    int temp45=0;
    while (getline(inputFile, line)) {
        temp45=stoi(line) % Q_ks;
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        Array1[ttt]=temp45;
        ttt++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile.close();
    // for (uint32_t i=0; i<N; i++)
    // {
    //     std::cout<<Array1[i]<<std::endl;
    // }
    // std::cout<<"Here is the value o ciphertext:"<<std::endl;
    int result54 =0, temp67=0,temp76=0;
    for (uint32_t i=0; i<N-1; i++)
    {   
        temp67=Array1[i+1];
        temp76 =(Q_ks-ctMS ->GetA()[i].ConvertToInt()); 
        // std::cout<<temp76<<std::endl;
        result54 = (result54 + temp76 * temp67) % Q_ks;
    }
    result54 = (result54 + ctMS ->GetA()[N-1].ConvertToInt() * Array1[0] + Q_ks ) % Q_ks;
    // std::cout<<"Result is: "<<result54<<"\t"<<ctMS ->GetB().ConvertToInt()<<std::endl;
    result54= (result54 % Q_ks + Q_ks) % Q_ks ;
    temp45 =ctMS->GetB().ConvertToInt(); // <<" , "<<std::endl;
    // std::cout<<"Result of result % Q_p is:"<<result<<"\t"<<temp<<std::endl;
    result54 = ((temp45-result54) % Q_ks +  Q_ks) % Q_ks ;
    
    result54 = ((result54+ Q_ks/(2*4)) % Q_ks); // .ModAddFastEq((mod / (qp * p * 2)), mod/qp);
    
    double result_temp01= (double) result54;
    double Q_ks1= (double) Q_ks;
    double result129= ((4*result_temp01) /(Q_ks1));
    // result129= ((4*result129) /(Q_ks1));
    // *result = ((NativeInteger(p) * r) / (mod/qp)).ConvertToInt();

    // std::cout<<"Result is:"<<result<<std::endl;
    // test=(Q/(2*4)).ConvertToDouble() * result;
    std::cout<<"Decrypted result of NTRU ciphertext is     :"<<"\t"<<result129<<std::endl;
    */
    
    //..........................................Code to decrypt NTRU ciphertext ends here .......................................

    // std::cout<<"I am here1:"<<std::endl;
    // NativePoly skNPoly(polyParams);
    // std::cout<<"I am here:2"<<std::endl;
    // skNPoly.SetValues(NatVec, Format::COEFFICIENT);
    // std::cout<<"I am here:3"<<std::endl;
    // skNPoly.SetFormat(Format::EVALUATION);
    // std::cout<<"I am here:"<<std::endl;
    // NativePoly ct(polyParams);
    // ct.SetValues(ctMS, Format::COEFFICIENT);
    // ct.SetFormat(Format::EVALUATION);
    
    //'''''''''''''''''Below lines are commented out.




    // ..........................................This code is to decrypt LWE ciphertext starts here .......................................
    /*
    uint32_t dim{params->GetLWEParams()->GetN()};
    int Q_ks = Q_ks.ConvertToInt();
    int Q_pp = ctKS->GetModulus().ConvertToInt();
    
    int Arraylwe[dim];
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<dim; i++)
    {
        Arraylwe[i]=0;
    }

    ifstream inputFile1("LWE_Secret_key1.txt");
    // Check if the file is successfully opened
    if (!inputFile1.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int tt=0;
    int temp11=0;
    while (getline(inputFile1, line)) {
        temp11=stoi(line);
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        Arraylwe[tt]=temp11;
        tt++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile1.close();
    // for (uint32_t i=0; i<dim; i++)
    // {
    //     std::cout<<Arraylwe[i]<<std::endl;
    // }
    int result1 =0, temp31=0, temp32=0;
    for (uint32_t i=0; i<dim; i++)
    {   
        temp31=Arraylwe[i];
        temp32 =ctKS ->GetA()[i].ConvertToInt(); // <<" , "<<std::endl;
        result1 = (result1 + temp31 * temp32)% Q_pp;
    }
    // std::cout<<"Result is: "<<result<<"\t"<<ctExt ->GetB().ConvertToInt()<<std::endl;
    result1= (result1 % Q_pp + Q_pp) % Q_pp ;
    temp31 =ctKS ->GetB().ConvertToInt(); // <<" , "<<std::endl;
    // std::cout<<"Result of result % Q_p is:"<<result<<"\t"<<temp<<std::endl;
    result1=((temp31-result1) % Q_pp +  Q_pp) % Q_pp ;

    // result1 = ((result1+ Q_pp/(2*4)) % Q_pp); // .ModAddFastEq((mod / (qp * p * 2)), mod/qp);
    double result_temp= (double) result1;
    double Q_ppp1= (double) Q_pp;
    double result12= ((4*result_temp) /( (double) Q_ppp1));
    

    std::cout<<"Decrypted result is:"<<result12<<std::endl;
    // test=(Q/(2*4)).ConvertToDouble() (4*result12);
    // std::cout<<"Decrypted result of LWE ciphertext is      :"<<"\t"<<result12<<std::endl;
    */
    //.......................................... Code to decrypt NTRU ciphertext ends here .......................................
    // auto ctksks=LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);

    // int Q_ppp = ctksks->GetModulus().ConvertToInt();
    //  auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    //.......................................... Code is to decrypt LWE ciphertext starts here .......................................
    /*
    int Arraylwelwe[dim];
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<dim; i++)
    {
        Arraylwelwe[i]=0;
    }

    ifstream inputFile2("LWE_Secret_key2.txt");
    // Check if the file is successfully opened
    if (!inputFile2.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    // string line;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int t=0;
    int temp1=0;
    while (getline(inputFile2, line)) {
        temp1=stoi(line);
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        Arraylwelwe[t]=temp1;
        t++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile2.close();
    // for (uint32_t i=0; i<dim; i++)
    // {
    //     std::cout<<Arraylwe[i]<<std::endl;
    // }
    // std::cout<<"Here is the result of a*s ffrom self dec:" <<std::endl;
    int result=0, temp=0;
    temp1=0;
    for (uint32_t i=0; i<dim; i++)
    {   
        temp1=Arraylwelwe[i];
        temp =ctksks ->GetA()[i].ConvertToInt(); // <<" , "<<std::endl;
        result = (result + temp * temp1)% Q_ppp;
        // std::cout<<result<<std::endl;
    }
    // std::cout<<"This is a from self dec:"<<ctksks ->GetA()[0]<<"\t"<<ctksks ->GetA()[10]<<"\t"<<ctksks ->GetA()[20]<<std::endl;
    // std::cout<<"This is s from self dec:"<<Arraylwelwe[0]<<"\t"<<Arraylwelwe[10]<<"\t"<<Arraylwelwe[20]<<std::endl;
    // std::cout<<"Q_ppp is is: "<<Q_ppp<<std::endl; //result<<"\t"<<ctExt ->GetB().ConvertToInt()<<std::endl;
    result= (result % Q_ppp+Q_ppp) % Q_ppp ;
    temp =ctksks ->GetB().ConvertToInt(); // <<" , "<<std::endl;
    
    int result15=((temp-result)% Q_ppp + Q_ppp ) % Q_ppp ;
    // std::cout<<"comp Result of (as,b,b-as) % Q_p is:"<<result<<"\t"<<temp<<"\t"<<result15<<std::endl;
    // int result1234 = ((result15+ Q_ppp/(2*4)) % Q_ppp); // .ModAddFastEq((mod / (qp * p * 2)), mod/qp);
    // std::cout<<"Self Impl computed b is::"<<result1234<< std::endl;

    double result_temp1= (double) result15; //result1234;
    double Q_pppp= (double) Q_ppp;
    double result123= ((4*result_temp1) /( (double) Q_pppp));
    

    // std::cout<<"Result is:"<<result<<std::endl;
    // test=(Q/(2*4)).ConvertToDouble() * result;
    std::cout<<"Decrypted result of LWE ciphertext is      :"<<"\t"<<result123<<std::endl;
    */
    //..........................................This code is to decrypt LWE ciphertext ends here .......................................

    // Modulus switching
    //LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
    // return LWEscheme->ModSwitchLwr(ct1->GetModulus(), qp, qpmodKS, ctKS);
    
}

LWECiphertext BinFHEScheme::EvalBinGateS(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const VectorNTRUBTKey& EK, ConstLWECiphertext& ct1,
                                        ConstLWECiphertext& ct2) const {

    // std::cout<<"Reached to EvalBinGate"<<std::endl;
    if (ct1 == ct2)
        OPENFHE_THROW(config_error, "Input ciphertexts should be independant");

    // By default, we compute XOR/XNOR using a combination of AND, OR, and NOT gates 迭代计算
    if ((gate == XOR) || (gate == XNOR)) {
        const auto& ctAND1 = EvalBinGateS(params, AND, EK, ct1, EvalNOT(params, ct2));
        const auto& ctAND2 = EvalBinGateS(params, AND, EK, EvalNOT(params, ct1), ct2);
        const auto& ctOR   = EvalBinGateS(params, OR, EK, ctAND1, ctAND2);

        // NOT is free so there is not cost to do it an extra time for XNOR
        return (gate == XOR) ? ctOR : EvalNOT(params, ctOR);
    }

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ct1);
    //构造(5q/8,0)
    auto n = params->GetLWEParams()->Getn();
    uint32_t N{params->GetLWEParams()->GetN()};

    NativeVector zero(n,0);
    uint32_t q  = params->GetLWEParams()->Getq().ConvertToInt<uint32_t>();
    uint32_t qp = params->GetLWEParams()->Getqp(); //.ConvertToInt<uint32_t>();
    // uint32_t qpmodKS = params->GetLWEParams()->GetqpKS(); //.ConvertToInt<uint32_t>();
    zero.SetModulus(q);

    //For LWE or lwr 
    // NativeInteger temp_b= 5*q/(8*qp);  // 5*p/(8) This had been changed due to new cipher mod p
    // std::cout<<"temp_b is="<<temp_b.Mod(q/qp)<<std::endl;
    // LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q/qp)));
    // LWECiphertext ct_temp = std::make_shared<LWECiphertextImpl>(LWECiphertextImpl(std::move(zero), temp_b.Mod(q/qp)));


    // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
    // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
    if ((gate == XOR_FAST) || (gate == XNOR_FAST)) {
        LWEscheme->EvalSubEq(qp, ctprep, ct2);
        LWEscheme->EvalAddEq(qp, ctprep, ctprep);
    }
    else {
        // LWEscheme->EvalAddEq(qp, ctprep, ct2);
        // LWEscheme->EvalSubEq(qp, ct_temp, ctprep);
        LWEscheme->EvalAddEq(qp, ctprep, ct2);

    }
    //This is for normal Bootstrapping using automorphisms
    // auto acc{BootstrapGateCore(params, gate, EK.BSkey, ct_temp)};  //一个NativePoly，多项式
    
    //This is for the window style bootstrapping need to have special property to satisy the isomorphic group
    auto acc{BootstrapGateCoreS(params, gate, EK.BSkey, ctprep)};  //一个NativePoly

    // auto acc{BootstrapGateCoreS(params, gate, EK.BSkey, ct_temp)};  //一个NativePoly，多项式
    // This is to check the computed output ..................................................................................................................
    /*
    NativePoly ct(acc->GetElements());
    ct.SetFormat(Format::COEFFICIENT);
    
    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    NativePoly sum(polyParams, Format::EVALUATION, true);

    uint32_t Q1{VNTRUParams->GetQ().ConvertToInt<uint32_t>()};
    NativeVector NtruSecO(N,0);
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<N; i++)
    {
        NtruSecO[i]=0;
    }
    ifstream inputFile1("Self_created_NTRU_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFile1.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line_1;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int tt_t=0;
    int temp045=0;
    while (getline(inputFile1, line_1)) {
        temp045=stoi(line_1) % Q1;
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        NtruSecO[tt_t]=temp045;
        tt_t++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile1.close();
    

    // NativeVector NatVec(N, Q);
    // NativeVector NatVec_inv(N, Q);
    // Get_invertible_NativeVector(NatVec, NatVec_inv, Q, N);
    // LWEPrivateKey LWEskN = std::make_shared<LWEPrivateKeyImpl>(LWEPrivateKeyImpl(NatVec));

    // ek.KSkey = LWEscheme->KeySwitchGen(LWEParams, LWEsk, LWEskN);

    NativePoly skNPoly(polyParams, Format::COEFFICIENT);
    skNPoly.SetValuesToZero();
    for (uint32_t i=0; i<N; i++){
        skNPoly[i]= NtruSecO[i];
        // std::cout<<"skNPoly[i]"<<zeroPoly[i]<<std::endl;
    }
    // NativePoly skNPoly(polyParams);
    // std::cout<<"working fine"<<std::endl;
    // skNPoly.SetValues(zeroPoly, Format::COEFFICIENT);
    // std::cout<<"working fine"<<std::endl;
    // NativePoly invskNPoly(polyParams);
    // invskNPoly.SetValues(NatVec_inv, Format::COEFFICIENT);
    skNPoly.SetFormat(Format::EVALUATION);
    ct.SetFormat(Format::EVALUATION);
    // invskNPoly.SetFormat(Format::EVALUATION);
    // std::cout<<"working fine"<<std::endl;
    sum=skNPoly*ct;
    sum.SetFormat(Format::COEFFICIENT);

    // std::cout<<"Content of the accumulator polynomial after NTRU decryption is:..........................."<<((sum[0]+(Q1 >> 3) + 1)*4)/Q1<<std::endl;
    // std::cout<<"Content of the accumulator polynomial after NTRU decryption is:..........................."<<(sum[0]*4)/Q1<<std::endl;
    std::cout<<"\nHere is the content of the NTRU decrypted ciphertext:"<<std::endl;
    for (uint32_t i=0; i<N; i++){
        std::cout<<sum[i]<<" , ";
    }
    std::cout<<" \n ";
    */
    // This is to check the computed output ends here ..................................................................................................................





    NativePoly& accVec{acc->GetElements()};
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    // std::cout<<"Content of accumulator vector is:"<<std::endl;
    // for (uint32_t i=0;i<N;i++)
    // {
    //     std::cout<<accVec[i]<<",";
    // }
    // std::cout<<"Before transpose accmulator content is:"<<std::endl;
    // for (u_int32_t i=0;i<N;i++)
    //     std::cout<<(acc->GetElements()[i]*4)/Q<<" , ";
    // std::cout<<"\n";
    accVec = accVec.Transpose();

    // std::cout<<"After transpose accmulator content is:"<<std::endl;
    // for (u_int32_t i=0;i<N;i++)
    //     std::cout<<(accVec[i]*4)/Q<<" , ";
    // std::cout<<"\n";


    // std::cout<<"After transpose"<<std::endl;
    // for (u_int32_t i=0;i<N;i++)
    //     std::cout<<"accVec["<<i<<"]"<<acc->GetElements()[i]<<std::endl;
    // std::cout<<"\n\n\nContent of accumulator vector transpose is:"<<std::endl;
    // for (uint32_t i=0;i<N;i++)
    // {
    //     std::cout<<accVec[i]<<",";
    // }
    accVec.SetFormat(Format::COEFFICIENT);
    // std::cout<<"\n\n\nContent of accumulator vector coefficient is:"<<std::endl;
    // for (uint32_t i=0;i<N;i++)
    // {
    //     std::cout<<accVec[i]<<",";
    // }
    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    // const auto& LWEParams = params->GetLWEParams();
    // NativeInteger Q{LWEParams->GetQ()};

    // uint32_t dim{LWEParams->Getn()};
    NativeInteger b{(Q >> 3) + 1};
    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec.GetValues()), std::move(b));
    
    //Check the returned ciphertext is correct ...........................................................................................
    // NativePoly ct(acc->GetElements());
    // ct.SetFormat(Format::COEFFICIENT);
    
    // const auto& VNTRUParams = params->GetVectorNTRUParams();
    // const auto& polyParams  = VNTRUParams->GetPolyParams();  
    // NativePoly sum(polyParams, Format::EVALUATION, true);

    // uint32_t Q1{VNTRUParams->GetQ().ConvertToInt<uint32_t>()};

    NativeVector NtruSec(N,0);
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<N; i++)
    {
        NtruSec[i]=0;
    }
    ifstream inputFile("Self_created_NTRU_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFile.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    // Read each line of the file and print it 
    NativeInteger Qhalf{Q>>1};
    string line;
    int tttt=0;
    int temp451=0;
    while (getline(inputFile, line)) {
        temp451=stoll(line);
        // std::cout<<"temp45 :"<<temp45<<std::endl;
        if (temp451 > Qhalf){
            NtruSec[tttt]=temp451-Q;
        }
        else{
            NtruSec[tttt]=temp451;
        }
        tttt++;
    }
    inputFile.close();
    NativeInteger temp{0};
     for (uint32_t i=0; i<N; i++){
        temp +=ctExt->GetA()[i]*NtruSec[i];   
    }
    temp=(temp+Q)%Q;
    NativeInteger bb=ctExt->GetB();
    bb=((bb-temp)+Q)%Q;
    std::cout<<"Decrypted result of NTRU/LWE ciphertext  is:"<<Q<<"\t\t"<<bb<<"\t\t"<<(((bb+Q/8)*4)/Q) % 4<<std::endl;
    //Checking completed ..................... ...........................................................................................

    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    












    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    //-----------------------------------------------------------------------I am going to decrypt here........................................................
        // NativeInteger n{ctKS->GetLWEParams->Getn()};
    NativeInteger Q3{ctKS->GetModulus()};
    NativeVector LweSec(n,0);
    ifstream inputFilel("LWE_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFilel.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line1;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int ttt=0;
    int temp45=0;
    while (getline(inputFilel, line)) {
        temp45=stoll(line);
        // std::cout<<"temp45 :"<<temp45<<std::endl;
        if (temp45 > Qhalf){
            NtruSec[ttt]=temp45-Q;
        }
        else{
            NtruSec[ttt]=temp45;
        }
        ttt++;
    }
    // Close the file
    inputFilel.close();
    NativeInteger temp2=0;
     for (uint32_t i=0; i<n; i++){
        temp2 +=ctKS->GetA()[i]*(LweSec[i]%Q3);
        // std::cout<<"skNPoly[i]"<<zeroPoly[i]<<std::endl;
    }
    if(qp>1){
        temp2=NativeInteger(static_cast<BasicInteger>( std::floor(0.5 + temp2.ConvertToDouble() * (Q3/qp).ConvertToDouble() / Q3.ConvertToDouble())));
        // RoundqQ( temp2, Q3/qp, Q3);
    }
    else{
        temp2=(temp2+Q3)%Q3;
    }
    
    NativeInteger bb2=ctKS->GetB();
    bb2=((bb2-temp2)+(Q3/qp))%(Q3/qp);
    // std::cout<<"Before ModSwitch Decrypted result of LWE cipher is:"<<Q3<<"\t\t"<<bb2<<"\t\t"<<((bb2+Q3/(qp*8))*4)/Q3<<std::endl;
    //-----------------------------------------------------------------------Decryption done I am going to decrypt here........................................................

    //Modswitch 
    if(qp>1){
        return LWEscheme->ModSwitchLwr(qp, ct1->GetModulus(), ctKS);
    }
    else{
        return LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
    }

    





    // std::cout<<"(Q, LWEParams->GetqKS(), ct1->GetModulus(), qp) are:("<<Q<<" , "<<LWEParams->GetqKS()<<" , "<<ct1->GetModulus()<<" , "<<qp<< " ) "<<std::endl;
    //define secret key
    //..........................................This code is to decrypt NTRU ciphertext.......................................
    // NativeInteger Q_ks{params->GetLWEParams()->GetqKS()};
    
    /*
    uint32_t N{params->GetLWEParams()->GetN()};
    // int Q_ks = Q_ks.ConvertToInt();
    // LWEParams->GetqKS().ConvertToInt();

    const auto& VNTRUParams = params->GetVectorNTRUParams();
    const auto& polyParams  = VNTRUParams->GetPolyParams();  
    int Q_ks = params->GetLWEParams()->GetqKS().ConvertToInt();
    
    int Array1[N];
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<N; i++)
    {
        Array1[i]=0;
    }
    ifstream inputFile("NTRU_Secret_key.txt");
    // Check if the file is successfully opened
    if (!inputFile.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int ttt=0;
    int temp45=0;
    while (getline(inputFile, line)) {
        temp45=stoi(line) % Q_ks;
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        Array1[ttt]=temp45;
        ttt++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile.close();
    // for (uint32_t i=0; i<N; i++)
    // {
    //     std::cout<<Array1[i]<<std::endl;
    // }
    // std::cout<<"Here is the value o ciphertext:"<<std::endl;
    int result54 =0, temp67=0,temp76=0;
    for (uint32_t i=0; i<N-1; i++)
    {   
        temp67=Array1[i+1];
        temp76 =(Q_ks-ctMS ->GetA()[i].ConvertToInt()); 
        // std::cout<<temp76<<std::endl;
        result54 = (result54 + temp76 * temp67) % Q_ks;
    }
    result54 = (result54 + ctMS ->GetA()[N-1].ConvertToInt() * Array1[0] + Q_ks ) % Q_ks;
    // std::cout<<"Result is: "<<result54<<"\t"<<ctMS ->GetB().ConvertToInt()<<std::endl;
    result54= (result54 % Q_ks + Q_ks) % Q_ks ;
    temp45 =ctMS->GetB().ConvertToInt(); // <<" , "<<std::endl;
    // std::cout<<"Result of result % Q_p is:"<<result<<"\t"<<temp<<std::endl;
    result54 = ((temp45-result54) % Q_ks +  Q_ks) % Q_ks ;
    
    result54 = ((result54+ Q_ks/(2*4)) % Q_ks); // .ModAddFastEq((mod / (qp * p * 2)), mod/qp);
    
    double result_temp01= (double) result54;
    double Q_ks1= (double) Q_ks;
    double result129= ((4*result_temp01) /(Q_ks1));
    // result129= ((4*result129) /(Q_ks1));
    // *result = ((NativeInteger(p) * r) / (mod/qp)).ConvertToInt();

    // std::cout<<"Result is:"<<result<<std::endl;
    // test=(Q/(2*4)).ConvertToDouble() * result;
    std::cout<<"Decrypted result of NTRU ciphertext is     :"<<"\t"<<result129<<std::endl;
    */
    
    //..........................................Code to decrypt NTRU ciphertext ends here .......................................

    // std::cout<<"I am here1:"<<std::endl;
    // NativePoly skNPoly(polyParams);
    // std::cout<<"I am here:2"<<std::endl;
    // skNPoly.SetValues(NatVec, Format::COEFFICIENT);
    // std::cout<<"I am here:3"<<std::endl;
    // skNPoly.SetFormat(Format::EVALUATION);
    // std::cout<<"I am here:"<<std::endl;
    // NativePoly ct(polyParams);
    // ct.SetValues(ctMS, Format::COEFFICIENT);
    // ct.SetFormat(Format::EVALUATION);
    
    //'''''''''''''''''Below lines are commented out.




    // ..........................................This code is to decrypt LWE ciphertext starts here .......................................
    /*
    uint32_t dim{params->GetLWEParams()->GetN()};
    int Q_ks = Q_ks.ConvertToInt();
    int Q_pp = ctKS->GetModulus().ConvertToInt();
    
    int Arraylwe[dim];
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<dim; i++)
    {
        Arraylwe[i]=0;
    }

    ifstream inputFile1("LWE_Secret_key1.txt");
    // Check if the file is successfully opened
    if (!inputFile1.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    string line;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int tt=0;
    int temp11=0;
    while (getline(inputFile1, line)) {
        temp11=stoi(line);
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        Arraylwe[tt]=temp11;
        tt++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile1.close();
    // for (uint32_t i=0; i<dim; i++)
    // {
    //     std::cout<<Arraylwe[i]<<std::endl;
    // }
    int result1 =0, temp31=0, temp32=0;
    for (uint32_t i=0; i<dim; i++)
    {   
        temp31=Arraylwe[i];
        temp32 =ctKS ->GetA()[i].ConvertToInt(); // <<" , "<<std::endl;
        result1 = (result1 + temp31 * temp32)% Q_pp;
    }
    // std::cout<<"Result is: "<<result<<"\t"<<ctExt ->GetB().ConvertToInt()<<std::endl;
    result1= (result1 % Q_pp + Q_pp) % Q_pp ;
    temp31 =ctKS ->GetB().ConvertToInt(); // <<" , "<<std::endl;
    // std::cout<<"Result of result % Q_p is:"<<result<<"\t"<<temp<<std::endl;
    result1=((temp31-result1) % Q_pp +  Q_pp) % Q_pp ;

    // result1 = ((result1+ Q_pp/(2*4)) % Q_pp); // .ModAddFastEq((mod / (qp * p * 2)), mod/qp);
    double result_temp= (double) result1;
    double Q_ppp1= (double) Q_pp;
    double result12= ((4*result_temp) /( (double) Q_ppp1));
    

    std::cout<<"Decrypted result is:"<<result12<<std::endl;
    // test=(Q/(2*4)).ConvertToDouble() (4*result12);
    // std::cout<<"Decrypted result of LWE ciphertext is      :"<<"\t"<<result12<<std::endl;
    */
    //.......................................... Code to decrypt NTRU ciphertext ends here .......................................
    // auto ctksks=LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);

    // int Q_ppp = ctksks->GetModulus().ConvertToInt();
    //  auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    //.......................................... Code is to decrypt LWE ciphertext starts here .......................................
    /*
    int Arraylwelwe[dim];
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<dim; i++)
    {
        Arraylwelwe[i]=0;
    }

    ifstream inputFile2("LWE_Secret_key2.txt");
    // Check if the file is successfully opened
    if (!inputFile2.is_open()) {
        cerr << "Error opening the file!" << endl;
    }
    // string line;
    // Read each line of the file and print it to the
    // standard output stream
    // cout << "File Content: " << endl;
    int t=0;
    int temp1=0;
    while (getline(inputFile2, line)) {
        temp1=stoi(line);
        // std::cout<<line<<"\t"<<temp1<<std::endl;
        Arraylwelwe[t]=temp1;
        t++;
        // NatVec[ii].ConvertToInt<uint32_t>();
        // std::cout << NatVec[ii] << endl; // Print the current line
    }
    // Close the file
    inputFile2.close();
    // for (uint32_t i=0; i<dim; i++)
    // {
    //     std::cout<<Arraylwe[i]<<std::endl;
    // }
    // std::cout<<"Here is the result of a*s ffrom self dec:" <<std::endl;
    int result=0, temp=0;
    temp1=0;
    for (uint32_t i=0; i<dim; i++)
    {   
        temp1=Arraylwelwe[i];
        temp =ctksks ->GetA()[i].ConvertToInt(); // <<" , "<<std::endl;
        result = (result + temp * temp1)% Q_ppp;
        // std::cout<<result<<std::endl;
    }
    // std::cout<<"This is a from self dec:"<<ctksks ->GetA()[0]<<"\t"<<ctksks ->GetA()[10]<<"\t"<<ctksks ->GetA()[20]<<std::endl;
    // std::cout<<"This is s from self dec:"<<Arraylwelwe[0]<<"\t"<<Arraylwelwe[10]<<"\t"<<Arraylwelwe[20]<<std::endl;
    // std::cout<<"Q_ppp is is: "<<Q_ppp<<std::endl; //result<<"\t"<<ctExt ->GetB().ConvertToInt()<<std::endl;
    result= (result % Q_ppp+Q_ppp) % Q_ppp ;
    temp =ctksks ->GetB().ConvertToInt(); // <<" , "<<std::endl;
    
    int result15=((temp-result)% Q_ppp + Q_ppp ) % Q_ppp ;
    // std::cout<<"comp Result of (as,b,b-as) % Q_p is:"<<result<<"\t"<<temp<<"\t"<<result15<<std::endl;
    // int result1234 = ((result15+ Q_ppp/(2*4)) % Q_ppp); // .ModAddFastEq((mod / (qp * p * 2)), mod/qp);
    // std::cout<<"Self Impl computed b is::"<<result1234<< std::endl;

    double result_temp1= (double) result15; //result1234;
    double Q_pppp= (double) Q_ppp;
    double result123= ((4*result_temp1) /( (double) Q_pppp));
    

    // std::cout<<"Result is:"<<result<<std::endl;
    // test=(Q/(2*4)).ConvertToDouble() * result;
    std::cout<<"Decrypted result of LWE ciphertext is      :"<<"\t"<<result123<<std::endl;
    */
    //..........................................This code is to decrypt LWE ciphertext ends here .......................................

    // Modulus switching
    //LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
    // return LWEscheme->ModSwitchLwr(ct1->GetModulus(), qp, qpmodKS, ctKS);
    
}

// Full evaluation as described in https://eprint.iacr.org/2020/086
LWECiphertext BinFHEScheme::EvalBinGate(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const RingGSWBTKey& EK, ConstLWECiphertext& ct1,
                                        ConstLWECiphertext& ct2) const {
    if (ct1 == ct2)
        OPENFHE_THROW(config_error, "Input ciphertexts should be independant");

    uint32_t qp = params->GetLWEParams()->Getqp();
    // uint32_t qpmodKS = params->GetLWEParams()->GetqpKS(); //.ConvertToInt<uint32_t>();
    // By default, we compute XOR/XNOR using a combination of AND, OR, and NOT gates 迭代计算
    if ((gate == XOR) || (gate == XNOR)) {
        const auto& ctAND1 = EvalBinGate(params, AND, EK, ct1, EvalNOT(params, ct2));
        const auto& ctAND2 = EvalBinGate(params, AND, EK, EvalNOT(params, ct1), ct2);
        const auto& ctOR   = EvalBinGate(params, OR, EK, ctAND1, ctAND2);

        // NOT is free so there is not cost to do it an extra time for XNOR
        return (gate == XOR) ? ctOR : EvalNOT(params, ctOR);
    }

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ct1);
    // the additive homomorphic operation for XOR/NXOR is different from the other gates we compute
    // 2*(ct1 - ct2) mod 4 for XOR, me map 1,2 -> 1 and 3,0 -> 0
    if ((gate == XOR_FAST) || (gate == XNOR_FAST)) {
        LWEscheme->EvalSubEq(qp, ctprep, ct2);
        LWEscheme->EvalAddEq(qp, ctprep, ctprep);
    }
    else {
        // for all other gates, we simply compute (ct1 + ct2) mod 4
        // for AND: 0,1 -> 0 and 2,3 -> 1
        // for OR: 1,2 -> 1 and 3,0 -> 0
        LWEscheme->EvalAddEq(qp, ctprep, ct2);
    }

    auto acc{BootstrapGateCore(params, gate, EK.BSkey, ctprep)};

    // the accumulator result is encrypted w.r.t. the transposed secret key
    // we can transpose "a" to get an encryption under the original secret key
    std::vector<NativePoly>& accVec{acc->GetElements()};
    accVec[0] = accVec[0].Transpose();  //对a重排列
    accVec[0].SetFormat(Format::COEFFICIENT);
    accVec[1].SetFormat(Format::COEFFICIENT);

    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    NativeInteger b{(Q >> 3) + 1};
    b.ModAddFastEq(accVec[1][0], Q);
    
    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(b));
    
    //Check the returned ciphertext is correct ...........................................................................................

    uint32_t n{params->GetLWEParams()->Getn()};
    uint32_t N{params->GetLWEParams()->GetN()};
    int QI=Q.ConvertToInt();
    int Qhalf{QI>>1};
    NativeInteger qKS{LWEParams->GetqKS()};
    int NtruSec[N]; 
    int Lwe_base_sec[n];
    // //to decrypt NTRU cipher 
    for (uint32_t i=0; i<N; i++)
    {
        NtruSec[i]=0;
    }
    for (uint32_t i=0; i<n; i++)
    {
        Lwe_base_sec[i]=0;
    }
    ifstream inputFileN("LWE_Secret_key1.txt");
    ifstream inputFileL("LWE_Secret_key_base.txt");
    // Check if the file is successfully opened
    if (!inputFileN.is_open()) {
        cerr << "Error opening the file NTRU!" << endl;
    }
    if (!inputFileL.is_open()) {
        cerr << "Error opening the file lwe!" << endl;
    }
    
    // Read each line of the file NtruSec and print it 
    string line;
    int ttt=0;
    int temp45=0;
    
    while (getline(inputFileN, line)) {
        temp45=stoll(line);
        if (temp45 > Qhalf){
            NtruSec[ttt]=temp45-QI;
        }
        else{
            NtruSec[ttt]=temp45;    
        }
        ttt++;
    }
    inputFileN.close();

    // Read each line of the file Lwe_base_sec and print it  
    string line1;
    int tttt=0;
    int temp445=0;
    int qKSI=qKS.ConvertToInt();
    int qKShalf=(qKSI>>1);
    while (getline(inputFileL, line1)) {
        temp445=stoll(line1);
        if (temp445 > qKShalf){
            int tiks=temp445-qKSI;
            // std::cout<<"tiks: "<<tiks<<"\t";
            Lwe_base_sec[tttt]=tiks;
            // std::cout<<"Lwe_base_sec[tttt]: "<<Lwe_base_sec[tttt]<<"\t";
        }
        else{
            // std::cout<<temp445<<std::endl;            
            Lwe_base_sec[tttt]=temp445;
        }
        tttt++;
    }
    inputFileL.close();

    // std::cout<<"Ntru secret key is:"<<"\t";
    // for (uint32_t i=0; i<N; i++){
    //     std::cout<<NtruSec[i]<<"\t";
    // }
    // std::cout<<"\n";

    // std::cout<<"Lwe_base_sec is:"<<"\t";
    // for (uint32_t i=0; i<n; i++){
    //     std::cout<<Lwe_base_sec[i]<<"\t";
    // }
    // std::cout<<"\n";

    NativeInteger temp{0};
    for (uint32_t i=0; i<N; i++){
        temp +=ctExt->GetA()[i]*NtruSec[i];   
    }
    temp=(temp+Q)%Q;
    NativeInteger bb=ctExt->GetB();
    std::cout<<"\ntemp0 and bb0 are:"<<temp<<"\t"<<bb<<std::endl;
    bb=((bb-temp)+Q)%Q;
    std::cout<<"Decrypted result of bootstrapped ciphertext  is:"<<Q<<"\t\t"<<bb<<"\t\t"<<((((bb+Q/8)*4)/Q) % 4)<<std::endl;
    //Checking completed ..................... ..........................................................................................................
   
    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    // auto ctMS = LWEscheme->ModSwitchLwr(qp, LWEParams->GetqKS(), ctExt);

    //Checking started ..................... ............................................................................................................
    NativeInteger temp1{0};
    for (uint32_t i=0; i<N; i++){
        temp1 += ctMS->GetA()[i] * NtruSec[i];   
    }
    temp1=(temp1+qKS) % qKS;
    NativeInteger bb1=ctMS->GetB();
    bb1=((bb1-temp1)+qKS) % qKS;
    std::cout<<"Decrypted result of ModSwitched ciphertext  is:"<<qKS<<"\t\t\t"<<bb1<<"\t\t\t"<<(((((bb1+qKS/8)*4)/qKS)) % 4)<<std::endl;
    //Checking completed ..................... ............................................................................................................


    // Key switching
    auto ctKS = LWEscheme->KeySwitchM(LWEParams, EK.KSkey, ctMS);
    
    NativeInteger temp2{0};
    for (uint32_t i=0; i<n; i++){
        temp2 += ctKS->GetA()[i] * Lwe_base_sec[i];   
    }
    temp2 = (((temp2 + qKS) % qKS)/qp) % (qKS/qp); 
    NativeInteger bb2 = ctKS->GetB();
    bb2=((bb2-temp2) + qKS/qp) % (qKS/qp);
    std::cout<<"Decrypted result of KeySwitched ciphertext is:"<<(qKS/qp)<<"\t\t\t"<<bb2<<"\t\t\t"<<((bb2+(qKS/(qp*8)))*4)/(qKS/qp) % 4<<std::endl;

    NativeInteger qpKS{params->GetLWEParams()->GetqpKS()};

    // Modulus switching
    if(qp>1){                                                           //for LWR     
        auto ctMSLwr = LWEscheme->ModSwitchLwr(qp, ct1->GetModulus(), ctKS);
        return ctMSLwr;
    }
    else{                                                               //for LWE     
        auto ctMSLwe = LWEscheme->ModSwitch(ct1->GetModulus(), ctKS);
        return ctMSLwe;

    }
}

// Full evaluation as described in https://eprint.iacr.org/2020/086
LWECiphertext BinFHEScheme::EvalBinGate(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                        const RingGSWBTKey& EK, const std::vector<LWECiphertext>& ctvector) const {
    // check if the ciphertexts are all independent
    for (size_t i = 0; i < ctvector.size(); i++) {
        for (size_t j = i + 1; j < ctvector.size(); j++) {
            if (ctvector[j] == ctvector[i]) {
                OPENFHE_THROW(config_error, "Input ciphertexts should be independent");
            }
        }
    }

    NativeInteger p = ctvector[0]->GetptModulus();
    uint32_t qp = params->GetLWEParams()->Getqp();

    LWECiphertext ctprep = std::make_shared<LWECiphertextImpl>(*ctvector[0]);
    ctprep->SetptModulus(p);
    if ((gate == MAJORITY) || (gate == AND3) || (gate == OR3) || (gate == AND4) || (gate == OR4)) {
        // we simply compute sum(ctvector[i]) mod p
        for (size_t i = 1; i < ctvector.size(); i++) {
            LWEscheme->EvalAddEq(qp, ctprep, ctvector[i]);
        }
        auto acc = BootstrapGateCore(params, gate, EK.BSkey, ctprep);

        std::vector<NativePoly>& accVec = acc->GetElements();
        // the accumulator result is encrypted w.r.t. the transposed secret key
        // we can transpose "a" to get an encryption under the original secret key
        accVec[0] = accVec[0].Transpose();
        accVec[0].SetFormat(Format::COEFFICIENT);
        accVec[1].SetFormat(Format::COEFFICIENT);

        // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
        auto& LWEParams = params->GetLWEParams();
        NativeInteger Q = LWEParams->GetQ();
        NativeInteger b = Q / NativeInteger(2 * p) + 1;
        b.ModAddFastEq(accVec[1][0], Q);

        auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(b));
        // Modulus switching to a middle step Q'
        auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
        // Key switching
        auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
        // Modulus switching
        return LWEscheme->ModSwitch(ctvector[0]->GetModulus(), ctKS);
    }
    else if (gate == CMUX) {
        if (ctvector.size() != 3)
            OPENFHE_THROW(not_implemented_error, "CMUX gate implemented for ciphertext vectors of size 3");

        auto ccNOT   = EvalNOT(params, ctvector[2]);
        auto ctNAND1 = EvalBinGate(params, NAND, EK, ctvector[0], ccNOT);
        auto ctNAND2 = EvalBinGate(params, NAND, EK, ctvector[1], ctvector[2]);
        auto ctCMUX  = EvalBinGate(params, NAND, EK, ctNAND1, ctNAND2);
        return ctCMUX;
    }
    else {
        OPENFHE_THROW(not_implemented_error, "This gate is not implemented for vector of ciphertexts at this time");
    }
}

// Full evaluation as described in https://eprint.iacr.org/2020/086
LWECiphertext BinFHEScheme::Bootstrap(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                      ConstLWECiphertext& ct) const {
    NativeInteger p = ct->GetptModulus();
    uint32_t qp = params->GetLWEParams()->Getqp();
    LWECiphertext ctprep{std::make_shared<LWECiphertextImpl>(*ct)};
    // ctprep = ct + q/4
    LWEscheme->EvalAddConstEq(qp, ctprep, (ct->GetModulus() >> 2));

    auto acc{BootstrapGateCore(params, AND, EK.BSkey, ctprep)};

    // the accumulator result is encrypted w.r.t. the transposed secret key
    // we can transpose "a" to get an encryption under the original secret key
    std::vector<NativePoly>& accVec{acc->GetElements()};
    accVec[0] = accVec[0].Transpose();
    accVec[0].SetFormat(Format::COEFFICIENT);
    accVec[1].SetFormat(Format::COEFFICIENT);

    // we add Q/8 to "b" to to map back to Q/4 (i.e., mod 2) arithmetic.
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger Q{LWEParams->GetQ()};
    NativeInteger b = Q / NativeInteger(2 * p) + 1;
    b.ModAddFastEq(accVec[1][0], Q);

    auto ctExt = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(b));
    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    // Modulus switching
    return LWEscheme->ModSwitch(ct->GetModulus(), ctKS);
}

// Evaluation of the NOT operation; no key material is needed
LWECiphertext BinFHEScheme::EvalNOT(const std::shared_ptr<BinFHECryptoParams>& params, ConstLWECiphertext& ct) const {
    NativeInteger q{ct->GetModulus()};
    uint32_t n{ct->GetLength()};
    uint32_t qp=params->GetLWEParams()->Getqp();

    NativeVector a(n, q);
    for (uint32_t i = 0; i < n; ++i)
        a[i] = ct->GetA(i) == 0 ? 0 : q - ct->GetA(i);
    return std::make_shared<LWECiphertextImpl>(std::move(a), (q/qp >> 2).ModSubFast(ct->GetB(), q/qp));
}

// Evaluate Arbitrary Function homomorphically
// Modulus of ct is q | 2N
LWECiphertext BinFHEScheme::EvalFunc(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                     ConstLWECiphertext& ct, const std::vector<NativeInteger>& LUT,
                                     const NativeInteger& beta) const {
    auto ct1 = std::make_shared<LWECiphertextImpl>(*ct);
    NativeInteger q{ct->GetModulus()};
    uint32_t functionProperty{this->checkInputFunction(LUT, q)};
    uint32_t qp = params->GetLWEParams()->Getqp();

    if (functionProperty == 0) {  // negacyclic function only needs one bootstrap
        auto fLUT = [LUT](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            return LUT[x.ConvertToInt()];
        };
        LWEscheme->EvalAddConstEq(qp, ct1, beta);
        return BootstrapFunc(params, EK, ct1, fLUT, q);
    }

    if (functionProperty == 2) {  // arbitary funciton
        const auto& LWEParams = params->GetLWEParams();
        uint32_t N{LWEParams->GetN()};
        if (q.ConvertToInt() > N) {  // need q to be at most = N for arbitary function
            std::string errMsg =
                "ERROR: ciphertext modulus q needs to be <= ring dimension for arbitrary function evaluation";
            OPENFHE_THROW(not_implemented_error, errMsg);
        }

        // TODO: figure out a way to not do this :(

        // repeat the LUT to make it periodic
        std::vector<NativeInteger> LUT2 = LUT;
        LUT2.insert(LUT2.end(), LUT.begin(), LUT.end());

        NativeInteger dq{q << 1};
        // raise the modulus of ct1 : q -> 2q
        ct1->GetA().SetModulus(dq);

        auto ct2 = std::make_shared<LWECiphertextImpl>(*ct1);
        LWEscheme->EvalAddConstEq(qp, ct2, beta);
        // this is 1/4q_small or -1/4q_small mod q
        auto f0 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            if (x < (q >> 1))
                return Q - (q >> 2);
            else
                return (q >> 2);
        };
        auto ct3 = BootstrapFunc(params, EK, ct2, f0, dq);
        LWEscheme->EvalSubEq2(qp, ct1, ct3);
        LWEscheme->EvalAddConstEq(qp, ct3, beta);
        LWEscheme->EvalSubConstEq(qp, ct3, q >> 1);

        // Now the input is within the range [0, q/2).
        // Note that for non-periodic function, the input q is boosted up to 2q
        auto fLUT2 = [LUT2](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            if (x < (q >> 1))
                return LUT2[x.ConvertToInt()];
            else
                return Q - LUT2[x.ConvertToInt() - q.ConvertToInt() / 2];
        };
        auto ct4 = BootstrapFunc(params, EK, ct3, fLUT2, dq);
        ct4->SetModulus(q);
        return ct4;
    }

    // Else it's periodic function so we evaluate directly
    LWEscheme->EvalAddConstEq(qp, ct1, beta);
    // this is 1/4q_small or -1/4q_small mod q
    auto f0 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 1))
            return Q - (q >> 2);
        else
            return (q >> 2);
    };
    auto ct2 = BootstrapFunc(params, EK, ct1, f0, q);
    LWEscheme->EvalSubEq2(qp, ct, ct2);
    LWEscheme->EvalAddConstEq(qp, ct2, beta);
    LWEscheme->EvalSubConstEq(qp, ct2, q >> 2);

    // Now the input is within the range [0, q/2).
    // Note that for non-periodic function, the input q is boosted up to 2q
    auto fLUT1 = [LUT](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 1))
            return LUT[x.ConvertToInt()];
        else
            return Q - LUT[x.ConvertToInt() - q.ConvertToInt() / 2];
    };
    return BootstrapFunc(params, EK, ct2, fLUT1, q);
}

// Evaluate Homomorphic Flooring
LWECiphertext BinFHEScheme::EvalFloor(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                      ConstLWECiphertext& ct, const NativeInteger& beta, uint32_t roundbits) const {
    const auto& LWEParams = params->GetLWEParams();
    NativeInteger q{roundbits == 0 ? LWEParams->Getq() : beta * (1 << roundbits + 1)};
    NativeInteger mod{ct->GetModulus()};
    uint32_t qp = params->GetLWEParams()->Getqp();
    auto ct1 = std::make_shared<LWECiphertextImpl>(*ct);
    LWEscheme->EvalAddConstEq(qp, ct1, beta);

    auto ct1Modq = std::make_shared<LWECiphertextImpl>(*ct1);
    ct1Modq->SetModulus(q);
    // this is 1/4q_small or -1/4q_small mod q
    auto f1 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 1))
            return Q - (q >> 2);
        else
            return (q >> 2);
    };
    auto ct2 = BootstrapFunc(params, EK, ct1Modq, f1, mod);
    LWEscheme->EvalSubEq(qp, ct1, ct2);

    auto ct2Modq = std::make_shared<LWECiphertextImpl>(*ct1);
    ct2Modq->SetModulus(q);

    // now the input is only within the range [0, q/2)
    auto f2 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
        if (x < (q >> 2))
            return Q - (q >> 1) - x;
        else if (((q >> 2) <= x) && (x < 3 * (q >> 2)))
            return x;
        else
            return Q + (q >> 1) - x;
    };
    auto ct3 = BootstrapFunc(params, EK, ct2Modq, f2, mod);
    LWEscheme->EvalSubEq(qp, ct1, ct3);

    return ct1;
}

// Evaluate large-precision sign
LWECiphertext BinFHEScheme::EvalSign(const std::shared_ptr<BinFHECryptoParams>& params,
                                     const std::map<uint32_t, RingGSWBTKey>& EKs, ConstLWECiphertext& ct,
                                     const NativeInteger& beta, bool schemeSwitch) const {
    auto mod{ct->GetModulus()};
    const auto& LWEParams = params->GetLWEParams();
    auto q{LWEParams->Getq()};
    uint32_t qp = params->GetLWEParams()->Getqp();

    if (mod <= q) {
        std::string errMsg =
            "ERROR: EvalSign is only for large precision. For small precision, please use bootstrapping directly";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    const auto& RGSWParams = params->GetRingGSWParams();
    const auto curBase     = RGSWParams->GetBaseG();
    auto search            = EKs.find(curBase);
    if (search == EKs.end()) {
        std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
        OPENFHE_THROW(openfhe_error, errMsg);
    }
    RingGSWBTKey curEK(search->second);

    auto cttmp = std::make_shared<LWECiphertextImpl>(*ct);
    while (mod > q) {
        cttmp = EvalFloor(params, curEK, cttmp, beta);
        // round Q to 2betaQ/q
        //  mod   = mod / q * 2 * beta;
        mod   = (mod << 1) * beta / q;
        cttmp = LWEscheme->ModSwitch(mod, cttmp);

        // if dynamic
        if (EKs.size() == 3) {
            // TODO: use GetMSB()?
            uint32_t binLog = static_cast<uint32_t>(ceil(GetMSB(mod.ConvertToInt()) - 1));
            uint32_t base{0};
            if (binLog <= static_cast<uint32_t>(17))
                base = static_cast<uint32_t>(1) << 27;
            else if (binLog <= static_cast<uint32_t>(26))
                base = static_cast<uint32_t>(1) << 18;

            if (0 != base) {  // if base is to change ...
                RGSWParams->Change_BaseG(base);

                auto search = EKs.find(base);
                if (search == EKs.end()) {
                    std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
                    OPENFHE_THROW(openfhe_error, errMsg);
                }
                curEK = search->second;
            }
        }
    }
    LWEscheme->EvalAddConstEq(qp, cttmp, beta);

    if (!schemeSwitch) {
        // if the ended q is smaller than q, we need to change the param for the final boostrapping
        auto f3 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            return (x < q / 2) ? (Q / 4) : (Q - Q / 4);
        };
        cttmp = BootstrapFunc(params, curEK, cttmp, f3, q);  // this is 1/4q_small or -1/4q_small mod q
        LWEscheme->EvalSubConstEq(qp, cttmp, q >> 2);
    }
    else {  // return the negated f3 and do not subtract q/4 for a more natural encoding in scheme switching
        // if the ended q is smaller than q, we need to change the param for the final boostrapping
        auto f3 = [](NativeInteger x, NativeInteger q, NativeInteger Q) -> NativeInteger {
            return (x < q / 2) ? (Q - Q / 4) : (Q / 4);
        };
        cttmp = BootstrapFunc(params, curEK, cttmp, f3, q);  // this is 1/4q_small or -1/4q_small mod q
    }
    RGSWParams->Change_BaseG(curBase);
    return cttmp;
}

// Evaluate Ciphertext Decomposition
std::vector<LWECiphertext> BinFHEScheme::EvalDecomp(const std::shared_ptr<BinFHECryptoParams>& params,
                                                    const std::map<uint32_t, RingGSWBTKey>& EKs, ConstLWECiphertext& ct,
                                                    const NativeInteger& beta) const {
    auto mod         = ct->GetModulus();
    auto& LWEParams  = params->GetLWEParams();
    auto& RGSWParams = params->GetRingGSWParams();

    NativeInteger q = LWEParams->Getq();
    if (mod <= q) {
        std::string errMsg =
            "ERROR: EvalDecomp is only for large precision. For small precision, please use bootstrapping directly";
        OPENFHE_THROW(not_implemented_error, errMsg);
    }

    const auto curBase = RGSWParams->GetBaseG();
    auto search        = EKs.find(curBase);
    if (search == EKs.end()) {
        std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
        OPENFHE_THROW(openfhe_error, errMsg);
    }
    RingGSWBTKey curEK(search->second);

    auto cttmp = std::make_shared<LWECiphertextImpl>(*ct);
    std::vector<LWECiphertext> ret;
    while (mod > q) {
        auto ctq = std::make_shared<LWECiphertextImpl>(*cttmp);
        ctq->SetModulus(q);
        ret.push_back(std::move(ctq));

        // Floor the input sequentially to obtain the most significant bit
        cttmp = EvalFloor(params, curEK, cttmp, beta);
        mod   = mod / q * 2 * beta;
        // round Q to 2betaQ/q
        cttmp = LWEscheme->ModSwitch(mod, cttmp);

        if (EKs.size() == 3) {  // if dynamic
            uint32_t binLog = static_cast<uint32_t>(ceil(log2(mod.ConvertToInt())));
            uint32_t base   = 0;
            if (binLog <= static_cast<uint32_t>(17))
                base = static_cast<uint32_t>(1) << 27;
            else if (binLog <= static_cast<uint32_t>(26))
                base = static_cast<uint32_t>(1) << 18;

            if (0 != base) {  // if base is to change ...
                RGSWParams->Change_BaseG(base);

                auto search = EKs.find(base);
                if (search == EKs.end()) {
                    std::string errMsg("ERROR: No key [" + std::to_string(curBase) + "] found in the map");
                    OPENFHE_THROW(openfhe_error, errMsg);
                }
                curEK = search->second;
            }
        }
    }
    RGSWParams->Change_BaseG(curBase);
    ret.push_back(std::move(cttmp));
    return ret;
}

// private:
NTRUCiphertext BinFHEScheme::BootstrapGateCore(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               ConstVectorNTRUACCKey& ek, ConstLWECiphertext& ct) const {
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();
    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger p  = ct->GetptModulus();                                     //4 明文模数
    NativeInteger q  = ct->GetModulus(); 
    NativeInteger qp     = LWEParams->Getqp();            //

    // depending on whether the value is the range, it will be set
    // to either Q/8 or -Q/8 to match binary arithmetic                                      //1024
    NativeInteger Q      = LWEParams->GetQ();             //
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1
    NativeInteger Q2pNeg = Q - Q2p;                       //7/8Q-1

    uint32_t N = LWEParams->GetN();
                   
    
    // Since q | (2*N), we deal with a sparse embedding of Z_Q[x]/(X^{q/2}+1) to
    // Z_Q[x]/(X^N+1)
    uint32_t factor = (2 * N)/(q.ConvertToInt()); 
    const NativeInteger& b = (qp*ct->GetB())*factor; // % (2*N); // 0~2N //making it domain q
    
    NativeVector m(N, Q);
    NativeVector new_m(N, Q);
    for (size_t j = 0; j < N; ++j) {
        m[j] = j<N/2 ?  Q2p:Q2pNeg;
    }
    for (size_t j = 0; j < N; ++j) {
        auto k = b.ConvertToInt()+j;
        if (k>=N && k<2*N )
        {
            new_m[k%N]=Q- m[j];
        }
        else
        {
             new_m[k%N]= m[j];
        }
    }

    // std::cout<<"Here is the value of accumulator function :"<<std::endl;
    // for (size_t i=0;i<N;i++){
    //     // std::cout<<(new_m[i]+Q/4)/Q<<" , ";
    //     std::cout<<new_m[i]<<" , ";
    // }
    // std::cout<<"\n";


    NativeInteger azero = ct->GetA()[0];
    uint32_t wzero = factor * azero.ConvertToInt() + 1;
    uint32_t invw = ModInverse(wzero, 2 * N) % (2 * N);

    NativePoly polym(polyParams);
    polym.SetValues(new_m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);
    auto polym2{polym.AutomorphismTransform(invw)};  
    auto acc = std::make_shared<NTRUCiphertextImpl>(std::move(polym2));
    NACCscheme->EvalAcc(NTRUParams, ek, acc, ct->GetA());
    return acc;
}

NTRUCiphertext BinFHEScheme::BootstrapGateCoreS(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               ConstVectorNTRUACCKey& ek, ConstLWECiphertext& ct) const {
    
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }
    
    auto& LWEParams  = params->GetLWEParams();
    auto& NTRUParams = params->GetVectorNTRUParams();
    auto polyParams  = NTRUParams->GetPolyParams();

    
    NativeInteger p  = ct->GetptModulus();                                     //4 明文模数
    NativeInteger q  = ct->GetModulus(); 
    NativeInteger qp = LWEParams->Getqp();            //

//    std::cout<<"vector a is:"<<ct->GetA()<<endl;
//    std::cout<<"vector b is:"<<ct->GetB()<<endl;

    // Specifies the range [q1,q2) that will be used for mapping
    uint32_t qHalf   = q.ConvertToInt() >> 1;
    // NativeInteger q11 = (3*q.ConvertToInt())/8; //RGSWParams->GetGateConst()[static_cast<size_t>(gate)];  //3/8q
    NativeInteger q1 = NTRUParams->GetGateConst()[static_cast<size_t>(gate)];
    // std::cout<<"q11 and q1 are:"<<q11<<"\t"<<q1<<std::endl;
    NativeInteger q2 = q1.ModAddFast(NativeInteger(qHalf), q);                 //7/8
    
    std::cout<<"Value of qp, q, q1 and q2 is:"<<qp<<","<<q<<","<<q1<<","<<q2<<std::endl;
    // depending on whether the value is the range, it will be set
    // to either Q/8 or -Q/8 to match binary arithmetic                                      //1024
    NativeInteger Q      = LWEParams->GetQ();             //
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1
    NativeInteger Q2pNeg = Q - Q2p;                       //7Q/8-1
    uint32_t N = LWEParams->GetN();
    NativeVector m(N, Q);

    uint32_t factor =(2 * N)/(q.ConvertToInt()); 
    const NativeInteger& b = qp*ct->GetB();

    // Since q | (2*N), we deal with a sparse embedding of Z_Q[x]/(X^{q/2}+1) to
    // Z_Q[x]/(X^N+1)
    // std::cout<<"Value of factor is:"<<factor<<std::endl;
    // std::cout<<"vector b is:"<<b<<endl;
    for (size_t j = 0; j < qHalf; ++j) {
        NativeInteger temp = b.ModSub(j, q);
        if (q1 < q2)
            m[j * factor] = ((temp >= q1) && (temp < q2)) ?  Q2p : Q2pNeg; //Q2pNeg : Q2p;
        else
            m[j * factor] = ((temp >= q2) && (temp < q1)) ? Q2pNeg : Q2p ; //Q2p : Q2pNeg;
    }

    // for (size_t j = 0; j < N; ++j) {
    //     std::cout<<m[j]<<" , ";
    // }        
    // std::cout<<" \n ";

    NativePoly polym(polyParams);
    // polym.SetValues(new_m, Format::COEFFICIENT);
    polym.SetValues(m, Format::COEFFICIENT);
    polym.SetFormat(EVALUATION);
    // auto polym2{polym.AutomorphismTransform(invw)};  
    auto acc = std::make_shared<NTRUCiphertextImpl>(std::move(polym));    
    // NACCscheme->EvalAccS(NTRUParams, ek, acc, ct->GetA());
    // NACCscheme->EvalAccTS(NTRUParams, ek, acc, ct->GetA());    
    NACCscheme->EvalAccTSC(NTRUParams, ek, acc, ct->GetA());    
    return acc;
}

RLWECiphertext BinFHEScheme::BootstrapGateCore(const std::shared_ptr<BinFHECryptoParams>& params, BINGATE gate,
                                               ConstRingGSWACCKey& ek, ConstLWECiphertext& ct) const {
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen "
            "before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& RGSWParams = params->GetRingGSWParams();
    auto polyParams  = RGSWParams->GetPolyParams();

    // Specifies the range [q1,q2) that will be used for mapping
    NativeInteger p  = ct->GetptModulus();  //4 明文模数
    NativeInteger q  = ct->GetModulus();
    NativeInteger qp = LWEParams->Getqp();
    uint32_t qHalf   = q.ConvertToInt() >> 1;
    NativeInteger q1 = RGSWParams->GetGateConst()[static_cast<size_t>(gate)];  //3/8q
    NativeInteger q2 = q1.ModAddFast(NativeInteger(qHalf), q);                 //7/8

    // depending on whether the value is the range, it will be set
    // to either Q/8 or -Q/8 to match binary arithmetic
    NativeInteger Q      = LWEParams->GetQ();
    NativeInteger Q2p    = Q / NativeInteger(2 * p) + 1;  //Q/8+1
    NativeInteger Q2pNeg = Q - Q2p;                       //7/8Q-1

    // std::cout<<"Value of qp, q, q1 and q2 is:"<<qp<<","<<q<<","<<q1<<","<<q2<<std::endl;
    // std::cout<<"Value of Q is:"<<Q<<std::endl;
    // std::cout<<"Value of Q2p is:"<<Q2p<<std::endl;
    // std::cout<<"Value of Q2pNeg is:"<<Q2pNeg<<std::endl;


    uint32_t N = LWEParams->GetN();
    NativeVector m(N, Q);
    // Since q | (2*N), we deal with a sparse embedding of Z_Q[x]/(X^{q/2}+1) to
    // Z_Q[x]/(X^N+1)
    uint32_t factor = (2 * N / q.ConvertToInt());
    
    const NativeInteger& b = (qp * ct->GetB()) % q;
    for (size_t j = 0; j < qHalf; ++j) {
        NativeInteger temp = b.ModSub(j, q);
        // std::cout<<"temp is:"<<temp<<endl;
        if (q1 < q2)
            m[j * factor] = ((temp >= q1) && (temp < q2)) ? Q2pNeg : Q2p; //Q2p : Q2pNeg; //
        else
            m[j * factor] = ((temp >= q2) && (temp < q1)) ? Q2p : Q2pNeg; //Q2pNeg : Q2p; 
    }
    // for (size_t j = 0; j < N; ++j) {
    //     std::cout<<m[j]<<" , ";
    // }        
    // std::cout<<" \n ";

    //m(x)-m(x^w)
    std::vector<NativePoly> res(2);
    // no need to do NTT as all coefficients of this poly are zero
    res[0] = NativePoly(polyParams, Format::EVALUATION, true);
    res[1] = NativePoly(polyParams, Format::COEFFICIENT, false);
    res[1].SetValues(std::move(m), Format::COEFFICIENT);
    res[1].SetFormat(Format::EVALUATION);

    // main accumulation computation
    // the following loop is the bottleneck of bootstrapping/binary gate evaluation
    //主累加计算
    //下面的循环是引导/二进制门的瓶颈
    //评估
    auto acc = std::make_shared<RLWECiphertextImpl>(std::move(res));
    // ACCscheme->EvalAcc(RGSWParams, ek, acc, ct->GetA());
    // ACCscheme->EvalAccTS(RGSWParams, ek, acc, ct->GetA());
    // ACCscheme->EvalAccTSS(RGSWParams, ek, acc, ct->GetA());
    return acc;
}

// Functions below are for large-precision sign evaluation, 大精度符号评估
// flooring, homomorphic digit decomposition, and arbitrary 向下取整、同态数字分解以及任意函数评估相关
// funciton evaluation, from https://eprint.iacr.org/2021/1337
template <typename Func>
RLWECiphertext BinFHEScheme::BootstrapFuncCore(const std::shared_ptr<BinFHECryptoParams>& params,
                                               ConstRingGSWACCKey& ek, ConstLWECiphertext& ct, const Func f,
                                               const NativeInteger& fmod) const {
    if (ek == nullptr) {
        std::string errMsg =
            "Bootstrapping keys have not been generated. Please call BTKeyGen before calling bootstrapping.";
        OPENFHE_THROW(config_error, errMsg);
    }

    auto& LWEParams  = params->GetLWEParams();
    auto& RGSWParams = params->GetRingGSWParams();
    auto polyParams  = RGSWParams->GetPolyParams();

    NativeInteger Q = LWEParams->GetQ();
    uint32_t N      = LWEParams->GetN();
    NativeVector m(N, Q);
    // For specific function evaluation instead of general bootstrapping
    NativeInteger ctMod    = ct->GetModulus();
    uint32_t factor        = (2 * N / ctMod.ConvertToInt());
    const NativeInteger& b = ct->GetB();
    for (size_t j = 0; j < (ctMod >> 1); ++j) {
        NativeInteger temp = b.ModSub(j, ctMod);
        m[j * factor]      = Q.ConvertToInt() / fmod.ConvertToInt() * f(temp, ctMod, fmod);
    }
    std::vector<NativePoly> res(2);
    // no need to do NTT as all coefficients of this poly are zero
    res[0] = NativePoly(polyParams, Format::EVALUATION, true);
    res[1] = NativePoly(polyParams, Format::COEFFICIENT, false);
    res[1].SetValues(std::move(m), Format::COEFFICIENT);
    res[1].SetFormat(Format::EVALUATION);

    // main accumulation computation
    // the following loop is the bottleneck of bootstrapping/binary gate
    // evaluation
    auto acc = std::make_shared<RLWECiphertextImpl>(std::move(res));
    ACCscheme->EvalAcc(RGSWParams, ek, acc, ct->GetA());
    return acc;
}

// Full evaluation as described in https://eprint.iacr.org/2020/086
template <typename Func>
LWECiphertext BinFHEScheme::BootstrapFunc(const std::shared_ptr<BinFHECryptoParams>& params, const RingGSWBTKey& EK,
                                          ConstLWECiphertext& ct, const Func f, const NativeInteger& fmod) const {
    auto acc = BootstrapFuncCore(params, EK.BSkey, ct, f, fmod);

    std::vector<NativePoly>& accVec = acc->GetElements();
    // the accumulator result is encrypted w.r.t. the transposed secret key
    // we can transpose "a" to get an encryption under the original secret key
    accVec[0] = accVec[0].Transpose();
    accVec[0].SetFormat(Format::COEFFICIENT);
    accVec[1].SetFormat(Format::COEFFICIENT);

    auto ctExt      = std::make_shared<LWECiphertextImpl>(std::move(accVec[0].GetValues()), std::move(accVec[1][0]));
    auto& LWEParams = params->GetLWEParams();
    // Modulus switching to a middle step Q'
    auto ctMS = LWEscheme->ModSwitch(LWEParams->GetqKS(), ctExt);
    // Key switching
    auto ctKS = LWEscheme->KeySwitch(LWEParams, EK.KSkey, ctMS);
    // Modulus switching
    return LWEscheme->ModSwitch(fmod, ctKS);
}
}
;  // namespace lbcrypto
