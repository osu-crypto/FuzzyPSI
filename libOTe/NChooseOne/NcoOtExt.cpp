#include "NcoOtExt.h"
#include "libOTe/Base/BaseOT.h"
#include "libOTe/TwoChooseOne/KosOtExtSender.h"
#include "libOTe/TwoChooseOne/KosOtExtReceiver.h"
#include "libOTe/TwoChooseOne/IknpOtExtSender.h"
#include "libOTe/TwoChooseOne/IknpOtExtReceiver.h"
#include <cryptoTools/Common/Matrix.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Network/Channel.h>

#include "coproto/NativeProto.h"
#include "coproto/Macros.h"
#include "libOTe/Tools/CoprotoSock.h"




void osuCrypto::NcoOtExtReceiver::genBaseOts(PRNG & prng, Channel & chl)
{
    CoprotoSock s(chl);
    auto ec = genBaseOts(prng).evaluate(s);
    if (ec)
        throw std::runtime_error("NcoOtExtReceiver::genBaseOts(), " + ec.message());
}

void osuCrypto::NcoOtExtSender::genBaseOts(PRNG& prng, Channel& chl)
{
    CoprotoSock s(chl);
    auto ec = genBaseOts(prng).evaluate(s);
    if (ec)
        throw std::runtime_error("NcoOtExtSender::genBaseOts(), " + ec.message());
}

void osuCrypto::NcoOtExtReceiver::setBaseOts(
    span<std::array<block, 2>> baseSendOts,
    PRNG& prng, Channel& chl)
{
    CoprotoSock s(chl);
    auto ec = setBaseOts(baseSendOts, prng).evaluate(s);
    if (ec)
        throw std::runtime_error("NcoOtExtReceiver::setBaseOts(), " + ec.message());
}

void osuCrypto::NcoOtExtReceiver::init(u64 numOtExt, PRNG& prng, Channel& chl)
{
    CoprotoSock s(chl);
    auto ec = init(numOtExt, prng).evaluate(s);
    if (ec)
        throw std::runtime_error("NcoOtExtReceiver::init(), " + ec.message());
}

coproto::Proto osuCrypto::NcoOtExtSender::setBaseOts(span<block> baseRecvOts, const BitVector& choices)
{
    struct Proto : public coproto::NativeProto
    {
        NcoOtExtSender& ot;
        span<block> baseRecvOts;
        const BitVector& choices;
        Proto(NcoOtExtSender&o,span<block> b, const BitVector& c)
            : ot(o)
            , baseRecvOts(b)
            , choices(c)
        {}

        block seed;
        PRNG prng;
        coproto::error_code resume() override
        {
            CP_BEGIN();

            if (ot.isMalicious())
            {
                CP_RECV(seed);
                prng.SetSeed(seed);
                BitVector bv(choices.size());
                bv.randomize(prng);
                bv ^= choices;

                ot.setUniformBaseOts(baseRecvOts, bv);

            }
            else
            {
                ot.setUniformBaseOts(baseRecvOts, choices);
            }

            CP_END();
            return {};
        }
    };

    return coproto::makeProto<Proto>(*this, baseRecvOts, choices);
}



coproto::Proto osuCrypto::NcoOtExtReceiver::setBaseOts(span<std::array<block, 2>> baseSendOts, PRNG& prng)
{

    struct Proto : public coproto::NativeProto
    {
        NcoOtExtReceiver& ot;
        span<std::array<block, 2>> baseSendOts;
        PRNG& prng;

        Proto(NcoOtExtReceiver& o, span<std::array<block, 2>> b,PRNG& p)
            : ot(o)
            , baseSendOts(b)
            , prng(p)
        {}

        block seed;
        PRNG seedPrng;
        coproto::error_code resume() override
        {
            CP_BEGIN();

            if (ot.isMalicious())
            {
                seed = prng.get();
                CP_SEND(seed);
                seedPrng.SetSeed(seed);
                BitVector bv(baseSendOts.size());
                bv.randomize(seedPrng);

                auto iter = bv.begin();
                std::vector< std::array<block, 2>> mm(baseSendOts.size());
                for (u64 i = 0; i < baseSendOts.size(); ++i)
                {
                    mm[i][0] = baseSendOts[i][*iter];
                    mm[i][1] = baseSendOts[i][*iter ^ 1];
                    ++iter;
                }

                ot.setUniformBaseOts(mm);

            }
            else
            {
                ot.setUniformBaseOts(baseSendOts);
            }

            CP_END();
            return {};
        }
    };

    return coproto::makeProto<Proto>(*this, baseSendOts, prng);
}


void osuCrypto::NcoOtExtSender::setBaseOts(span<block> baseRecvOts, const BitVector& choices, Channel& chl)
{
    CoprotoSock s(chl);
    auto ec = setBaseOts(baseRecvOts, choices).evaluate(s);
    if (ec)
        throw std::runtime_error("NcoOtExtSender::setBaseOts(), " + ec.message());
}

void osuCrypto::NcoOtExtSender::init(u64 numOtExt, PRNG& prng, Channel& chl)
{
    CoprotoSock s(chl);
    auto ec = init(numOtExt, prng).evaluate(s);
    if (ec)
        throw std::runtime_error("NcoOtExtSender::init(), " + ec.message());
}



coproto::Proto osuCrypto::NcoOtExtReceiver::genBaseOts(PRNG& prng)
{
    struct Proto : public coproto::NativeProto
    {
        NcoOtExtReceiver& ot;
        PRNG& prng;
        Proto(NcoOtExtReceiver& o, PRNG& p)
            :ot(o)
            ,prng(p)
        {}

        std::vector<std::array<block, 2>> msgs;
        bool mal;
#ifdef ENABLE_IKNP
        IknpOtExtSender iknp;
#endif
#ifdef ENABLE_KOS
        KosOtExtSender kos;
#elif defined LIBOTE_HAS_BASE_OT
        DefaultBaseOT base;
#endif

        coproto::error_code resume() override
        {
            CP_BEGIN();

            msgs.resize(ot.getBaseOTCount());
            mal = ot.isMalicious();
#ifdef ENABLE_IKNP
            if (!mal)
            {
                CP_AWAIT(iknp.send(msgs, prng));
                ot.setUniformBaseOts(msgs);
                return {};
            }
#endif

#ifdef ENABLE_KOS
            CP_AWAIT(kos.send(msgs, prng));
#elif defined LIBOTE_HAS_BASE_OT
            CP_AWAIT(base.send(msgs, prng));
#else
            throw std::runtime_error("The libOTe library does not have base OTs. Enable them to call this. " LOCATION);
#endif
            if(mal)
                CP_AWAIT(ot.setBaseOts(msgs, prng));
            else
                ot.setUniformBaseOts(msgs);


            CP_END();
            return {};
        }
    };

    return coproto::makeProto<Proto>(*this, prng);
}

coproto::Proto osuCrypto::NcoOtExtSender::genBaseOts(PRNG& prng)
{
    struct Proto : public coproto::NativeProto
    {
        NcoOtExtSender& ot;
        PRNG& prng;

        Proto(NcoOtExtSender& o,
            PRNG& p)
            :ot(o)
            ,prng(p)
        {}

        bool mal;
        std::vector<block> msgs;
        BitVector bv;
#ifdef ENABLE_IKNP
        IknpOtExtReceiver iknp;
#endif
#ifdef ENABLE_KOS
        KosOtExtReceiver kos;
#elif defined(LIBOTE_HAS_BASE_OT)
        DefaultBaseOT base;
#endif
        coproto::error_code resume() override
        {
            CP_BEGIN();
            
            msgs.resize(ot.getBaseOTCount());
            bv.resize(msgs.size());
            bv.randomize(prng);
            mal = ot.isMalicious();
#ifdef ENABLE_IKNP
            if (!mal)
            {
                CP_AWAIT(iknp.receive(bv, msgs, prng));
                ot.setUniformBaseOts(msgs, bv);
                return {};
            }
#endif

#ifdef ENABLE_KOS
            CP_AWAIT(kos.receive(bv, msgs, prng));
#elif defined LIBOTE_HAS_BASE_OT
            CP_AWAIT(base.receive(bv, msgs, prng));
#else 
            throw std::runtime_error("The libOTe library does not have base OTs. Enable them to call this. " LOCATION);
#endif
            if(mal)
                CP_AWAIT(ot.setBaseOts(msgs, bv));
            else
                ot.setUniformBaseOts(msgs, bv);

            CP_END();
            return{};
        }
    };

    return coproto::makeProto<Proto>(*this, prng);
}


void osuCrypto::NcoOtExtSender::check(Channel& chl, block seed)
{
    CoprotoSock s(chl);
    auto ec = check(seed).evaluate(s);
    if (ec)
        throw std::runtime_error("NcoOtExtSender::check(), " + ec.message());
}

coproto::Proto osuCrypto::NcoOtExtSender::sendChosen(MatrixView<block> messages, PRNG& prng)
{

    struct Proto : public coproto::NativeProto
    {
        NcoOtExtSender& ot;
        MatrixView<block> messages;
        PRNG& prng;

        Proto(NcoOtExtSender& o, MatrixView<block> m, PRNG& p)
            : ot(o)
            , messages(m)
            , prng(p)
        {}
        // must be at least 128 bits.
        Matrix<block> temp;

        coproto::error_code resume() override
        {

            CP_BEGIN();


            if (ot.hasBaseOts() == false)
                CP_AWAIT(ot.genBaseOts(prng));

            CP_AWAIT(ot.init(messages.rows(), prng));
            CP_AWAIT(ot.recvCorrection(messages.rows()));

            if (ot.isMalicious())
            {
                CP_AWAIT(ot.check(prng.get<block>()));
            }

            {
                std::array<u64, 2> choice{ 0,0 };
                u64& j = choice[0];

                temp.resize(messages.rows(), messages.cols());
                for (u64 i = 0; i < messages.rows(); ++i)
                {
                    for (j = 0; j < messages.cols(); ++j)
                    {
                        ot.encode(i, choice.data(), &temp(i, j));
                        temp(i, j) = temp(i, j) ^ messages(i, j);
                    }
                }
            }

            CP_SEND(std::move(temp));

            CP_END();
            return {};
        }
    };

    return coproto::makeProto<Proto>(*this, messages, prng);
}

void osuCrypto::NcoOtExtSender::sendChosen(MatrixView<block> messages, PRNG & prng, Channel & chl)
{

    auto ec = evaluate(sendChosen(messages, prng), chl);
    if (ec)
        throw std::runtime_error("NcoOtExtSender::sendChosen(), " + ec.message());
}

void osuCrypto::NcoOtExtReceiver::check(Channel& chl, block seed)
{
    auto ec = evaluate(check(seed),chl);
    if (ec)
        throw std::runtime_error("NcoOtExtReceiver::check(), " + ec.message());
}

coproto::Proto osuCrypto::NcoOtExtReceiver::receiveChosen(u64 numMsgsPerOT, span<block> messages, span<u64> choices, PRNG& prng)
{
    struct Proto : public coproto::NativeProto
    {
        NcoOtExtReceiver& ot;
        u64 numMsgsPerOT;
        span<block> messages; 
        span<u64> choices;
        PRNG& prng;

        Proto(NcoOtExtReceiver& o, u64 n, span<block> m, span<u64> c, PRNG& p)
            :ot(o)
            , numMsgsPerOT(n)
            , messages(m)
            , choices(c)
            , prng(p)
        {}
        // must be at least 128 bits.
        std::array<u64, 2> choiceVec{ 0,0 };
        Matrix<block> temp;

        coproto::error_code resume() override
        {
            u64& choice = choiceVec[0];
                
            CP_BEGIN();

            if (ot.hasBaseOts() == false)
                CP_AWAIT(ot.genBaseOts(prng));


            CP_AWAIT(ot.init(messages.size(), prng));

            for (i64 i = 0; i < messages.size(); ++i)
            {
                choice = choices[i];
                ot.encode(i, &choice, &messages[i]);
            }

            CP_AWAIT(ot.sendCorrection(messages.size()));


            if (ot.isMalicious())
            {
                CP_AWAIT(ot.check(prng.get<block>()));
            }

            temp.resize(messages.size(), numMsgsPerOT);
            CP_RECV(temp);
            //chl.recv(temp.data(), temp.size());

            for (i64 i = 0; i < messages.size(); ++i)
            {
                messages[i] = messages[i] ^ temp(i, choices[i]);
            }

            CP_END();
            return {};
        }
    };

    return coproto::makeProto<Proto>(*this, numMsgsPerOT, messages, choices, prng);
}

void osuCrypto::NcoOtExtReceiver::receiveChosen(
    u64 numMsgsPerOT, 
    span<block> messages, 
    span<u64> choices, PRNG & prng, Channel & chl)
{
    auto ec = evaluate(receiveChosen(numMsgsPerOT, messages, choices, prng), chl);
    if (ec)
        throw std::runtime_error("NcoOtExtReceiver::receiveChosen(), " + ec.message());
}
