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

#ifndef _RGSW_ACC_LMKCDEY_H_
#define _RGSW_ACC_LMKCDEY_H_

#include "rgsw-acc.h"

#include <memory>

namespace lbcrypto {

/**
 * @brief Ring GSW accumulator schemes described in
 * https://eprint.iacr.org/2022/198
 */
class RingGSWAccumulatorLMKCDEY final : public RingGSWAccumulator {
public:
    RingGSWAccumulatorLMKCDEY() = default;

    /**
   * Key generation for internal Ring GSW as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skNTT secret key polynomial in the EVALUATION representation
   * @param LWEsk the secret key
   * @return a shared pointer to the resulting keys
   */
    RingGSWACCKey KeyGenAcc(const std::shared_ptr<RingGSWCryptoParams>& params, const NativePoly& skNTT,
                            ConstLWEPrivateKey& LWEsk) const override;

    /**
   * Main accumulator function used in bootstrapping - LMKCDEY variant
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek the accumulator key
   * @param acc previous value of the accumulator
   * @param a value to update the accumulator with
   */
    void EvalAcc(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek, RLWECiphertext& acc,
                 const NativeVector& a) const override;
                
    /**
   * Main accumulator function used in bootstrapping - LMKCDEY variant with reduced acc key switching
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek the accumulator key
   * @param acc previous value of the accumulator
   * @param a value to update the accumulator with
   */
    void EvalAccTS(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek, RLWECiphertext& acc,
                 const NativeVector& a) const override;

   /*
   * Main accumulator function used in bootstrapping - LMKCDEY variant with reduced acc key switching
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek the accumulator key
   * @param acc previous value of the accumulator
   * @param a value to update the accumulator with
   */
    void EvalAccTSW(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek, RLWECiphertext& acc,
                 const NativeVector& a) const override;

   /* 
   * Main accumulator function used in bootstrapping - LMKCDEY variant with reduced acc key switching
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek the accumulator key
   * @param acc previous value of the accumulator
   * @param a value to update the accumulator with
   */
    void EvalAccTSS(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWACCKey& ek, RLWECiphertext& acc,
                 const NativeVector& a) const override;

    /**
   * Key combination for each set 
   * @param params a shared pointer to RingGSW scheme parameters
   * @param SkeySq secret key square polynomial in the EVALUATION representation
   * @param key encryption of the secret key
   * @return a shared pointer pointing to the multiplication of keys
   * Create RGSW ciphertext of the given RLWE' ciphertext
   */
    RingGSWACCKey RLWEtoRGSW(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                     ConstRingGSWEvalKey& SkeySq, ConstRingGSWEvalKey& key) const override;

    /**
   * This function is used to multiply RLWE and RGSW ciphertexts
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek encryption of secret key element 
   * @param rgswkey key is the encryption of sk.m and m in two components of the key 
   * @return a shared pointer pointing to the multiplication of two keys resulting one single key
   */

    RingGSWEvalKey MulRLWEtoRGSW(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                ConstRingGSWEvalKey& ek, ConstRingGSWEvalKey& rgswkey0, ConstRingGSWEvalKey& rgswkey1) const override;

private:
    /**
   * LMKCDEY Key generation for internal Ring GSW as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skNTT secret key polynomial in the EVALUATION representation
   * @param m a plaintext
   * @return a shared pointer to the resulting keys
   */
    RingGSWEvalKey KeyGenLMKCDEY(const std::shared_ptr<RingGSWCryptoParams>& params, const NativePoly& skNTT,
                                 LWEPlaintext m) const;
    /**
   * LMKCDEYC Key generation for internal Ring GSW as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skNTT secret key polynomial in the EVALUATION representation
   * @param m a plaintext
   * @param t a representing total rotation for automorphism computation
   * @return a shared pointer to the resulting keys
   */
    RingGSWEvalKey KeyGenLMKCDEYC(const std::shared_ptr<RingGSWCryptoParams>& params, const NativePoly& skNTT, 
                                LWEPlaintext m, uint32_t j, LWEPlaintext t) const;


    /**
   * LMKCDEY Key generation for internal Ring GSW as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skNTT secret key polynomial in the EVALUATION representation
   * @return a shared pointer to the resulting keys
   */

    RingGSWEvalKey KeyGenSecSqr(const std::shared_ptr<RingGSWCryptoParams>& params,
                                                     const NativePoly& skNTT) const;
    /**
   * Automorphism keys generation for internal Ring GSW as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param skNTT secret key polynomial in the EVALUATION representation
   * @param k a plaintext
   * @return a shared pointer to the resulting keys
   */
    RingGSWEvalKey KeyGenAuto(const std::shared_ptr<RingGSWCryptoParams>& params, const NativePoly& skNTT,
                              LWEPlaintext k) const;

    /**
   * LMKCDEY Accumulation as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek evaluation key for Ring GSW
   * @param acc previous value of the accumulator
   * @return
   */
    void AddToAccLMKCDEY(const std::shared_ptr<RingGSWCryptoParams>& params, ConstRingGSWEvalKey& ek,
                         RLWECiphertext& acc) const;
    /**
   * LMKCDEY Accumulation accompaning key switching key as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param ek evaluation key for Ring GSW
   * @param acc previous value of the accumulator
   * @return
   */
    void AddToAccLMKCDEY_A(const std::shared_ptr<RingGSWCryptoParams>& params, const NativeInteger& a,
                        ConstRingGSWEvalKey& ek, RLWECiphertext& acc) const;


    /**
   * LMKCDEY Accumulation automorphism evaluation as described in https://eprint.iacr.org/2022/198
   *
   * @param params a shared pointer to RingGSW scheme parameters
   * @param a index
   * @param ak evaluation key for Ring GSW
   * @param acc previous value of the accumulator
   * @return
   */
    void Automorphism(const std::shared_ptr<RingGSWCryptoParams>& params, const NativeInteger& a,
                      ConstRingGSWEvalKey& ak, RLWECiphertext& acc) const;
};

}  // namespace lbcrypto

#endif  // _RGSW_ACC_LMKCDEY_H_
