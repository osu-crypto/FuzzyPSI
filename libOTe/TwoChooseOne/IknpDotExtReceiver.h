#pragma once
// This file and the associated implementation has been placed in the public domain, waiving all copyright. No restrictions are placed on its use. 
#include "libOTe/config.h"
#ifdef ENABLE_DELTA_IKNP

#include "libOTe/TwoChooseOne/OTExtInterface.h"
#include <array>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Common/Timer.h>
#include "libOTe/Tools/LinearCode.h"
#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"

namespace osuCrypto
{

    class IknpDotExtReceiver :
        public IknpOtExtReceiver
    {
    public:

        IknpDotExtReceiver()
        {
            mDeltaOT = true;
        }
        IknpDotExtReceiver(const IknpDotExtReceiver&) = delete;
        IknpDotExtReceiver(IknpDotExtReceiver&&) = default;

        IknpDotExtReceiver(span<std::array<block, 2>> baseSendOts)
        {
            mDeltaOT = true;
            setUniformBaseOts(baseSendOts);
        }

        IknpDotExtReceiver &operator=(IknpDotExtReceiver&& v) = default;


        // returns an independent instance of this extender which can securely be
        // used concurrently to this current one. The base OTs for the new instance 
        // are derived from the orginial base OTs.
        IknpDotExtReceiver splitBase();

        // returns an independent (type eased) instance of this extender which can securely be
        // used concurrently to this current one. The base OTs for the new instance 
        // are derived from the orginial base OTs.
        std::unique_ptr<OtExtReceiver> split() override;


    };

}
#endif