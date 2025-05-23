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

#include "binfhecontext.h"

using namespace lbcrypto;

int main() {
    // Sample Program: Step 1: Set CryptoContext
    auto cc = BinFHEContext();

    cc.GenerateBinFHEContext(P192G, XZDDF);

    // Sample Program: Step 2: Key Generation
    auto sk = cc.KeyGen();
    int m0=0;
    int m1=1;
    LWEPlaintext result1; //, result1, result;  
    LWECiphertext ctAND1;

    // LWECiphertext ctnot1;
    // LWECiphertext ctnot2;
    // Generate the bootstrapping keys (refresh and switching keys)
    // std::cout << "Generating the bootstrapping keys..." << std::endl;
    cc.NBTKeyGen(sk);
    NativeVector s   = sk->GetElement();
    // NativeInteger q  = cc.GetParams()->GetLWEParams()->Getq();
    // NativeInteger n  = cc.GetParams()->GetLWEParams()->Getn();
    // NativeInteger q_ks  = cc.GetParams()->GetLWEParams()->GetqKS();
    // // std::cout<<"q is:"<<q<<std::endl;
    // write secret key into a file

    // //write NTRU secret key into a file
    // std::ofstream outputFile1("LWE_Secret_key1.txt");  // Open/create a file named "test.txt" for writing
    // if (outputFile1.is_open()) {  // Check if the file was successfully opened
    //     // Write some text into the file
    //     for (uint32_t i=0; i<n; i++)
    //     {
    //         outputFile1 << (s[i]+q_ks)%q_ks <<std::endl;
    //     }
    // }
    // std::ofstream outputFile2("LWE_Secret_key2.txt");  // Open/create a file named "test.txt" for writing
    // if (outputFile2.is_open()) {  // Check if the file was successfully opened
    //     // Write some text into the file
    //     for (uint32_t i=0; i<n; i++)
    //     {
    //         outputFile2 << (s[i]+q) % q <<std::endl;
    //     }
    // }

    // std::cout << "Completed the key generation." << std::endl;

    for (int i=0; i<10; i++){

        // Sample Program: Step 3: Encryption
        auto ct1 = cc.Encrypt(sk, m0);
        auto ct2 = cc.Encrypt(sk, m1);

        // ctnot1 =cc.EvalNOT(ct1);
        // ctnot2 =cc.EvalNOT(ct2);

        // cc.Decrypt(sk, ct1, &result1);
        // cc.Decrypt(sk, ctnot1, &result1);
        // cc.Decrypt(sk, ctnot2, &result2);

        // cc.Decrypt(sk, ctnot2, &result2);

        // std::cout << "Result of not computation of encrypted computation of (0,1) is :("<<result1<<" , " <<result2<<" ) "<<std::endl;
        // m0<<" NAND "<<m1<<" ) = " << result << std::endl;

        // Sample Program: Step 4: Evaluation
        // std::cout << "...................................Start  the  gate bootstrapping................................... " << std::endl;

        //Notice: We have only made specific modifications for NAND gates, and will add other gates in the future.
        ctAND1 = cc.EvalBinGate(NAND, ct1, ct2);
        cc.Decrypt(sk, ctAND1, &result1);
        // cc.m_NBTKey
        // NativeInteger res = cc.EvalBinGate(NAND, ct1, ct2);
        
        
        std::cout << "Result of encrypted computation of ( "<<m0<<" NAND "<<m1<<" ) = " << result1 << std::endl;
        // std::cout << "...................................end  the  gate bootstrapping................................... " << std::endl;
    }    
    return 0;
}
