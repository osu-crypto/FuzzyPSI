#include "OTExtInterface.h"
#include "libOTe/Base/BaseOT.h"
#include <cryptoTools/Common/BitVector.h>
#include <vector>
#include <cryptoTools/Network/Channel.h>
#include "coproto/NativeProto.h"
#include "coproto/Macros.h"
#include "libOTe/Tools/CoprotoSock.h"

coproto::Proto osuCrypto::OtExtReceiver::setBaseOts(span<std::array<block, 2>> baseSendOts, PRNG& prng)
{
    struct Proto : public coproto::NativeProto
    {
        OtExtReceiver& ot;
        span<std::array<block, 2>> baseSendOts;
        Proto(OtExtReceiver& o, span<std::array<block, 2>> b)
            :ot(o), baseSendOts(b) {}

        coproto::error_code resume() override
        {
            ot.setUniformBaseOts(baseSendOts);
            return {};
        }
    };

    return coproto::makeProto<Proto>(*this, baseSendOts);
}

void osuCrypto::OtExtReceiver::genBaseOts(PRNG & prng, Channel & chl)
{
    CoprotoSock s(chl);
    auto ec = genBaseOts(prng).evaluate(s);
    if (ec)
        throw std::runtime_error("genBaseOts(), " + ec.message());
}

coproto::Proto osuCrypto::OtExtReceiver::genBaseOts(PRNG& prng)
{
#ifdef LIBOTE_HAS_BASE_OT
    struct Base : coproto::NativeProto
    {
        OtExtReceiver& ot;
        PRNG& prng;
        DefaultBaseOT base;
        std::vector<std::array<block, 2>> msgs;

        Base(OtExtReceiver& o, PRNG& p)
            :ot(o)
            ,prng(p)
        {}

        coproto::error_code resume() override
        {
            CP_BEGIN();
            msgs.resize(ot.baseOtCount());
            CP_AWAIT(base.send(msgs, prng));
            CP_AWAIT(ot.setBaseOts(msgs, prng));
            CP_END();
            return {};
        }
    }; 
     

    return coproto::makeProto<Base>(*this, prng);

#else
    throw std::runtime_error("The libOTe library does not have base OTs. Enable them to call this. " LOCATION);
#endif
}

coproto::Proto osuCrypto::OtExtSender::setBaseOts(span<block> baseRecvOts, const BitVector& choices)
{
    struct Proto : public coproto::NativeProto
    {
        OtExtSender& ot;
        span<block> baseRecvOts;
        const BitVector& choices;
        Proto(OtExtSender& o, span<block> b, const BitVector& c)
            :ot(o), baseRecvOts(b), choices(c) {}

        coproto::error_code resume() override
        {
            ot.setUniformBaseOts(baseRecvOts, choices);
            return {};
        }
    };

    return coproto::makeProto<Proto>(*this, baseRecvOts, choices);

}

void osuCrypto::OtExtSender::genBaseOts(PRNG & prng, Channel & chl)
{
#ifdef LIBOTE_HAS_BASE_OT
    DefaultBaseOT base;
    auto count = baseOtCount();
    std::vector<block> msgs(count);
    BitVector bv(count);
    bv.randomize(prng);

    base.receive(bv, msgs, prng, chl);
    setBaseOts(msgs, bv, chl);
#else
    throw std::runtime_error("The libOTe library does not have base OTs. Enable them to call this. " LOCATION);
#endif

}

coproto::Proto osuCrypto::OtExtSender::genBaseOts(PRNG& prng)
{
    struct Base : coproto::NativeProto
    {
        OtExtSender& ot;
        PRNG& prng;
        DefaultBaseOT base;
        std::vector<block> msgs;
        BitVector bv;

        Base(OtExtSender& o, PRNG& p)
            :ot(o)
            , prng(p)
        {}

        coproto::error_code resume() override
        {
            CP_BEGIN();
            msgs.resize(ot.baseOtCount());
            bv.resize(msgs.size());
            bv.randomize(prng);

            CP_AWAIT(base.receive(bv, msgs, prng));
            CP_AWAIT(ot.setBaseOts(msgs, bv));
            CP_END();
            return {};
        }
    };


    return coproto::makeProto<Base>(*this, prng);
    //return coproto::Proto();
}

void osuCrypto::OtReceiver::receiveChosen(
    const BitVector & choices, 
    span<block> recvMessages,
    PRNG & prng, 
    Channel & chl)
{
    receive(choices, recvMessages, prng, chl);
    std::vector<std::array<block,2>> temp(recvMessages.size());
    chl.recv(temp.data(), temp.size());
    auto iter = choices.begin();
    for (u64 i = 0; i < temp.size(); ++i)
    {
        recvMessages[i] = recvMessages[i] ^ temp[i][*iter];
        ++iter;
    }
}

void osuCrypto::OtReceiver::receiveCorrelated(const BitVector& choices, span<block> recvMessages, PRNG& prng, Channel& chl)
{
    receive(choices, recvMessages, prng, chl);
    std::vector<block> temp(recvMessages.size());
    chl.recv(temp.data(), temp.size());
    auto iter = choices.begin();
    
    for (u64 i = 0; i < temp.size(); ++i)
    {
        recvMessages[i] = recvMessages[i] ^ (zeroAndAllOne[*iter] & temp[i]);
        ++iter;
    }

}

void osuCrypto::OtSender::sendChosen(
    span<std::array<block, 2>> messages, 
    PRNG & prng, 
    Channel & chl)
{
    std::vector<std::array<block, 2>> temp(messages.size());
    send(temp, prng, chl);

    for (u64 i = 0; i < static_cast<u64>(messages.size()); ++i)
    {
        temp[i][0] = temp[i][0] ^ messages[i][0];
        temp[i][1] = temp[i][1] ^ messages[i][1];
    }

    chl.asyncSend(std::move(temp));
}
