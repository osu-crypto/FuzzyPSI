#pragma once
#include "libOTe/config.h"
#ifdef ENABLE_MR

#include "libOTe/TwoChooseOne/OTExtInterface.h"
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Crypto/PRNG.h>

#include "coproto/Proto.h"
#include "libOTe/Tools/CoprotoSock.h"

namespace osuCrypto
{

    class MasnyRindal : public OtReceiver, public OtSender
    {
    public:

        coproto::Proto receive(
            const BitVector& choices,
            span<block> messages,
            PRNG& prng)override;



        coproto::Proto send(
            span<std::array<block, 2>> messages,
            PRNG& prng)override;

        using OtReceiver::receive;
        using OtSender::send;


        void receive(
            const BitVector& choices,
            span<block> messages,
            PRNG& prng,
            Channel& chl,
            u64 numThreads)
        {
            receive(choices, messages, prng, chl);
        }

        void send(
            span<std::array<block, 2>> messages,
            PRNG& prng,
            Channel& chl,
            u64 numThreads)
        {
            send(messages, prng, chl);
        }

        //void receive(
        //    const BitVector& choices,
        //    span<block> messages,
        //    PRNG& prng,
        //    Channel& chl) override
        //{
        //    auto sock = CoprotoSock(chl);
        //    receive(choices, messages, prng).evaluate(sock);
        //}

        //void send(
        //    span<std::array<block, 2>> messages,
        //    PRNG& prng,
        //    Channel& chl) override
        //{
        //    auto sock = CoprotoSock(chl);
        //    send(messages, prng).evaluate(sock);
        //}
    };

}
#endif