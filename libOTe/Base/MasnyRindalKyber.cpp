#include "MasnyRindalKyber.h"
#ifdef ENABLE_MR_KYBER

#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Network/Channel.h>
#include "coproto/NativeProto.h"
#include "coproto/Macros.h"
#include "libOTe/Tools/CoprotoSock.h"

namespace osuCrypto
{

    void MasnyRindalKyber::receive(
        const BitVector & choices, 
        span<block> messages, 
        PRNG & prng, 
        Channel & chl)
    {
        CoprotoSock s(chl);
        receive(choices, messages, prng).evaluate(s);
    }

    coproto::Proto MasnyRindalKyber::receive(const BitVector& choices, span<block> messages, PRNG& prng)
    {
        struct MPProto : public coproto::NativeProto
        {
            const BitVector& choices;
            span<block> messages;
            PRNG& prng;
            MPProto(const BitVector& c, span<block> m, PRNG& p)
                :choices(c)
                ,messages(m)
                ,prng(p)
            {}

            u64 n;
            std::vector<KyberOTRecver> ot;
            std::vector<KyberOtRecvPKs> pkBuff;
            static_assert(std::is_trivial<KyberOtRecvPKs>::value, "");
            std::vector<KyberOTCtxt> ctxts;
            KyberOtRecvPKs* iter;
            coproto::error_code resume() override
            {
                CP_BEGIN();

                n = choices.size();

                ot.resize(n);

                pkBuff.resize(n);

                iter = pkBuff.data();

                for (u64 i = 0; i < n; ++i)
                {
                    ot[i].b = choices[i];

                    //get receivers message and secret coins
                    KyberReceiverMessage(&ot[i], iter++);
                }

                CP_SEND(std::move(pkBuff));


                static_assert(std::is_pod<KyberOTCtxt>::value, "");
                ctxts.resize(n);

                CP_RECV(ctxts);

                for (u64 i = 0; i < n; ++i)
                {
                    KyberReceiverStrings(&ot[i], &ctxts[i]);
                    memcpy(&messages[i], ot[i].rot, sizeof(block));
                }

                CP_END();
                return {};
            }
        };

        return coproto::makeProto<MPProto>(choices, messages, prng);
    }

    void MasnyRindalKyber::send(
        span<std::array<block, 2>> messages, 
        PRNG & prng, 
        Channel & chl)
    {
        CoprotoSock s(chl);
        send(messages, prng).evaluate(s);
    }
    coproto::Proto MasnyRindalKyber::send(span<std::array<block, 2>> messages, PRNG& prng)
    {
        struct MPProto : public coproto::NativeProto
        {
            span<std::array<block, 2>> messages;
            PRNG& prng;
            MPProto(span<std::array<block, 2>> m, PRNG& p)
                :messages(m)
                ,prng(p)
            {}

            u64 n;
            std::vector<KyberOtRecvPKs> pkBuff;
            std::vector<KyberOTCtxt> ctxts;

            coproto::error_code resume() override
            {
                CP_BEGIN();
                n = messages.size();
                pkBuff.resize(n);
                ctxts.resize(n);

                prng.get(messages.data(), messages.size());


                CP_RECV(pkBuff);
                KyberOTPtxt ptxt;

                for (u64 i = 0; i < n; ++i)
                {
                    memcpy(ptxt.sot[0], &messages[i][0], sizeof(block));
                    memset(ptxt.sot[0] + sizeof(block), 0, sizeof(ptxt.sot[0]) - sizeof(block));

                    memcpy(ptxt.sot[1], &messages[i][1], sizeof(block));
                    memset(ptxt.sot[1] + sizeof(block), 0, sizeof(ptxt.sot[1]) - sizeof(block));

                    //get senders message, secret coins and ot strings
                    KyberSenderMessage(&ctxts[i], &ptxt, &pkBuff[i]);
                }

                CP_SEND(std::move(ctxts));

                CP_END();
                return {};
            }
        };

        return coproto::makeProto<MPProto>(messages, prng);
    }
}
#endif