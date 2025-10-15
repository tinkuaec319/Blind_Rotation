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
  Example for the FHEW scheme using the AP bootstrapping
 */

#include "binfhecontext.h"
using namespace std::chrono;  
using namespace lbcrypto;

int main() {
    // Sample Program: Step 1: Set CryptoContext

    auto cc = BinFHEContext();

    // STD128 is the security level of 128 bits of security based on LWE Estimator
    // and HE standard. Other common options are TOY, MEDIUM, STD192, and STD256.
    // MEDIUM corresponds to the level of more than 100 bits for both quantum and
    // classical computer attacks. The second argument is the bootstrapping method
    // (AP or GINX). The default method is GINX. Here we explicitly set AP. GINX
    // typically provides better performance: the bootstrapping key is much
    // smaller in GINX (by 20x) while the runtime is roughly the same.
    cc.GenerateBinFHEContext(STD128, AP);
    // cc.GenerateBinFHEContext(STD128_LMKCDEY, AP);
       

    // Sample Program: Step 2: Key Generation

    // Generate the secret key
    auto sk = cc.KeyGen();

    std::cout << "Generating the bootstrapping keys..." << std::endl;

    // Generate the bootstrapping keys (refresh and switching keys)
    cc.BTKeyGen(sk);

    std::cout << "Completed the key generation." << std::endl;

    // Sample Program: Step 3: Encryption

    // Encrypt two ciphertexts representing Boolean True (1)
    // By default, freshly encrypted ciphertexts are bootstrapped.
    // If you wish to get a fresh encryption without bootstrapping, write
    // auto   ct1 = cc.Encrypt(sk, 1, FRESH);
    auto ct1 = cc.Encrypt(sk, 1);
    auto ct2 = cc.Encrypt(sk, 1);

    // Sample Program: Step 4: Evaluation

    // Compute (1 AND 1) = 1; Other binary gate options are OR, NAND, and NOR
    auto ctAND1 = cc.EvalBinGate(AND, ct1, ct2);

    // Compute (NOT 1) = 0
    auto ct2Not = cc.EvalNOT(ct2);

    // Compute (1 AND (NOT 1)) = 0
    auto ctAND2 = cc.EvalBinGate(AND, ct2Not, ct1);

    // Computes OR of the results in ctAND1 and ctAND2 = 1
    auto ctResult = cc.EvalBinGate(OR, ctAND1, ctAND2);

    // Sample Program: Step 5: Decryption

    LWEPlaintext result;

    cc.Decrypt(sk, ctResult, &result);
    std::cout << "Result of encrypted computation of (1 AND 1) OR (1 AND (NOT 1)) = " << result << std::endl;


    int m1=1;
    int m2=1;
    auto total_time=0;
    int Enc_total=1000;

    for (int i=0; i<Enc_total; i++){

        // Sample Program: Step 3: Encryption
        auto ct1 = cc.Encrypt(sk, m1);
        auto ct2 = cc.Encrypt(sk, m2);
        
        // LWEPlaintext result1;
        // LWEPlaintext result2;
        // cc.Decrypt(sk, ct1, &result1);
        // cc.Decrypt(sk, ct2, &result2);
        // std::cout << "Result of enc/dec of.............................................................. ("<<m1<<" , "<<m2<<" )= "<<"( "<<result1<<" , "<<result2<<" )"<<std::endl;

    
        LWEPlaintext result;
        LWECiphertext ctAND1;
        // // Sample Program: Step 4: Evaluation
        // std::cout << "Start  the  gate bootstrapping " << std::endl;
        //Notice: We have only made specific modifications for NAND gates, and will add other gates in the future.
        auto start = high_resolution_clock::now();        
        ctAND1 = cc.EvalBinGate(NAND, ct1, ct2);   
        auto end = high_resolution_clock::now();
        auto comp_time=duration_cast<microseconds>(end-start);
        // std::cout<<"Time taken is:"<<comp_time.count()<<std::endl;
        total_time=total_time + comp_time.count();

        cc.Decrypt(sk, ctAND1, &result);

        // std::cout << "Result of enc/dec of 1/0 is:"<<result1<<"/"<<result2<<std::endl;
        //std::cout << "Result of enc/dec of........................................................... ("<<m1<<" NAND "<<m2<<")=(0)--> "<<result<<std::endl;        
        if (result)
        {
            std::cout<<"Decryption incorrect....please have a look on the code"<<std::endl;
            break;
        }
        
    }    
    std::cout<<"Total computation time is:"<<(total_time/(1000*Enc_total))<<std::endl;    

    return 0;
}
