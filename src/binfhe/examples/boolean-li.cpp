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
#include <chrono>
using namespace std::chrono;  
using namespace lbcrypto;

int main() {
    // Sample Program: Step 1: Set CryptoContext
    auto cc = BinFHEContext();

    // cc.GenerateBinFHEContext(P192G_LWE, XZDDF);
    // cc.GenerateBinFHEContext(P192G_LWR, XZDDF);
    // cc.GenerateBinFHEContext(P128G_LWR, XZDDF);    
    cc.GenerateBinFHEContext(P128G, XZDDF);        

    // Sample Program: Step 2: Key Generation
    auto sk = cc.KeyGen();

    // int n=cc.GetParams()->GetLWEParams()->Getn();
    // NativeInteger q=cc.GetParams()->GetLWEParams()->Getq();
    // std::cout<<"n,q is:"<<n<<"\t"<<q<<std::endl;
    
    //Save LWE secret key....
    //write NTRU secret key into a file
    // std::ofstream outputFile("LWE_Secret_key.txt");  // Open/create a file named "test.txt" for writing
    // if (outputFile.is_open()) {  // Check if the file was successfully opened
    //     // Write some text into the file
    //     for (int i=0; i<n ; i++)
    //     {
    //         outputFile << sk->GetElement()[i] <<std::endl;
    //     }
    // }
    // Close the file
    // outputFile.close();  // Close the file after writing

    int m1=1;
    int m2=1;
    auto total_time=0;
    int Enc_total=1000;

    // Generate the bootstrapping keys (refresh and switching keys)
    std::cout << "Generating the bootstrapping keys..." << std::endl;
    
    //This is for norml bootstrapping .... this do not requires any special property/NTRU dimension
    // cc.NBTKeyGen(sk);
    
    //This is for winow based bootstrapping .... this requires power 2 NTRU dimension
    cc.NBTKeyGenS(sk);

    std::cout << "Completed the key generation." << std::endl;

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
        // std::cout << "Result of enc/dec of........................................................... ("<<m1<<" NAND "<<m2<<")=(0)--> "<<result<<std::endl;        
        // if (result)
        // {
        //     std::cout<<"Decryption incorrect....please have a look on the code"<<std::endl;
        //     break;
        // }
        
    }    
    std::cout<<"Total computation time is:"<<(total_time/(1000*Enc_total))<<std::endl;    
    return 0;
}
