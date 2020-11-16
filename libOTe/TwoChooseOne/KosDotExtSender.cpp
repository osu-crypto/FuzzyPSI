#include "KosDotExtSender.h"
#ifdef ENABLE_DELTA_KOS

#include "libOTe/Tools/Tools.h"
#include <cryptoTools/Common/Matrix.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Crypto/Commit.h>
#include <cryptoTools/Network/Channel.h>
#include "TcoOtDefines.h"
#include "coproto/NativeProto.h"
#include <coproto/Macros.h>

namespace osuCrypto
{
    //#define KOS_DEBUG

    using namespace std;

    KosDotExtSender KosDotExtSender::splitBase()
    {
        std::vector<block> baseRecvOts(mGens.size());
        for (u64 i = 0; i < mGens.size(); ++i)
            baseRecvOts[i] = mGens[i].get<block>();

        return KosDotExtSender(baseRecvOts, mBaseChoiceBits);
    }

    std::unique_ptr<OtExtSender> KosDotExtSender::split()
    {
        std::vector<block> baseRecvOts(mGens.size());
        for (u64 i = 0; i < mGens.size(); ++i)
            baseRecvOts[i] = mGens[i].get<block>();

        return std::make_unique<KosDotExtSender>(baseRecvOts, mBaseChoiceBits);
    }

    void KosDotExtSender::setUniformBaseOts(span<block> baseRecvOts, const BitVector& choices)
    {
        mBaseChoiceBits = choices;
        mGens.resize(choices.size());
        mBaseChoiceBits.resize(roundUpTo(mBaseChoiceBits.size(), 8));
        for (u64 i = mBaseChoiceBits.size() - 1; i >= choices.size(); --i)
            mBaseChoiceBits[i] = 0;

        mBaseChoiceBits.resize(choices.size());
        for (u64 i = 0; i < mGens.size(); i++)
            mGens[i].SetSeed(baseRecvOts[i]);
    }

    void KosDotExtSender::setDelta(const block& delta)
    {
        mDelta = delta;
    }

    coproto::Proto KosDotExtSender::send(
        span<std::array<block, 2>> messages,
        PRNG& prng)
    {
        struct Proto : public coproto::NativeProto
        {
            KosDotExtSender& ot;
            span<std::array<block, 2>> messages;
            PRNG& prng;
            Proto(KosDotExtSender& o, span<std::array<block, 2>> m, PRNG& p)
                : ot(o)
                , messages(m)
                , prng(p)
            {}


            // round up
            u64 numOtExt, numSuperBlocks, superBlkIdx, step, doneIdx = 0,
                xx, end;

            // a temp that will be used to transpose the sender's matrix
            Matrix<u8> t;
            std::vector<std::array<block, superBlkSize>> u;

            std::vector<block> choiceMask;
            std::array<block, 2> delta{ ZeroBlock, ZeroBlock };

            std::array<std::array<block, 2>, 128> extraBlocks;
            std::array<block, 2>* xIter;
            Commit theirSeedComm;
            span<std::array<block, 2>>::iterator mIter, mIterPartial;
            block* uIter, * uEnd, * cIter, * tIter;
            block seed, curDelta, offset, theirSeed;
            PRNG codePrng, commonPrng;
            LinearCode code;
            std::array<block, 4> q, qi, tt;
            span<u8> data;
            std::array<block, 128> challenges, challenges2;
            std::vector<block> proof;
            coproto::error_code resume() override
            {

                CP_BEGIN();

                if (ot.hasBaseOts() == false)
                    CP_AWAIT(ot.genBaseOts(prng));

                ot.setTimePoint("KosDot.send.start");

                // round up
                numOtExt = roundUpTo(messages.size(), 128);
                numSuperBlocks = (numOtExt / 128 + superBlkSize) / superBlkSize;

                // a temp that will be used to transpose the sender's matrix
                t.resize(ot.mGens.size(), superBlkSize * sizeof(block));
                u.resize(ot.mGens.size() * commStepSize);

                choiceMask.resize(ot.mBaseChoiceBits.size());

                memcpy(delta.data(), ot.mBaseChoiceBits.data(), ot.mBaseChoiceBits.sizeBytes());


                for (u64 i = 0; i < choiceMask.size(); ++i)
                {
                    if (ot.mBaseChoiceBits[i]) choiceMask[i] = AllOneBlock;
                    else choiceMask[i] = ZeroBlock;
                }


                xIter = extraBlocks.data();

                CP_RECV(theirSeedComm);

                mIter = messages.begin();
                mIterPartial = messages.end() - std::min<u64>(128 * superBlkSize, messages.size());


                // set uIter = to the end so that it gets loaded on the first loop.
                uIter = (block*)u.data() + superBlkSize * ot.mGens.size() * commStepSize;
                uEnd = uIter;

                for (superBlkIdx = 0; superBlkIdx < numSuperBlocks; ++superBlkIdx)
                {

                    if (uIter == uEnd)
                    {
                        step = std::min<u64>(numSuperBlocks - superBlkIdx, (u64)commStepSize);
                        data = span<u8>((u8*)u.data(), step * superBlkSize * ot.mGens.size() * sizeof(block));

                        CP_RECV(data);
                        uIter = (block*)u.data();
                    }

                    cIter = choiceMask.data();
                    tIter = (block*)t.data();

                    // transpose 128 columns at at time. Each column will be 128 * superBlkSize = 1024 bits long.
                    for (u64 colIdx = 0; colIdx < ot.mGens.size(); ++colIdx)
                    {
                        // generate the columns using AES-NI in counter mode.
                        ot.mGens[colIdx].mAes.ecbEncCounterMode(ot.mGens[colIdx].mBlockIdx, superBlkSize, tIter);
                        ot.mGens[colIdx].mBlockIdx += superBlkSize;

                        uIter[0] = uIter[0] & *cIter;
                        uIter[1] = uIter[1] & *cIter;
                        uIter[2] = uIter[2] & *cIter;
                        uIter[3] = uIter[3] & *cIter;
                        uIter[4] = uIter[4] & *cIter;
                        uIter[5] = uIter[5] & *cIter;
                        uIter[6] = uIter[6] & *cIter;
                        uIter[7] = uIter[7] & *cIter;

                        tIter[0] = tIter[0] ^ uIter[0];
                        tIter[1] = tIter[1] ^ uIter[1];
                        tIter[2] = tIter[2] ^ uIter[2];
                        tIter[3] = tIter[3] ^ uIter[3];
                        tIter[4] = tIter[4] ^ uIter[4];
                        tIter[5] = tIter[5] ^ uIter[5];
                        tIter[6] = tIter[6] ^ uIter[6];
                        tIter[7] = tIter[7] ^ uIter[7];

                        ++cIter;
                        uIter += 8;
                        tIter += 8;
                    }



                    if (mIter >= mIterPartial)
                    {
                        Matrix<u8> tOut(128 * superBlkSize, sizeof(block) * 2);

                        // transpose our 128 columns of 1024 bits. We will have 1024 rows,
                        // each 128 bits wide.
                        transpose(t, tOut);

                        auto mCount = std::min<u64>(128 * superBlkSize, messages.end() - mIter);
                        auto xCount = std::min<u64>(128 * superBlkSize - mCount, extraBlocks.data() + extraBlocks.size() - xIter);


                        //std::copy(mIter, mIter + mCount, tOut.begin());
                        if (mCount) memcpy(&*mIter, tOut.data(), mCount * sizeof(block) * 2);
                        mIter += mCount;


                        memcpy(xIter, tOut.data() + mCount * sizeof(block) * 2, xCount * sizeof(block) * 2);
                        xIter += xCount;
                    }
                    else
                    {
                        MatrixView<u8> tOut(
                            (u8*)&*mIter,
                            128 * superBlkSize,
                            sizeof(block) * 2);

                        mIter += std::min<u64>(128 * superBlkSize, messages.end() - mIter);

                        // transpose our 128 columns of 1024 bits. We will have 1024 rows,
                        // each 128 bits wide.
                        transpose(t, tOut);
                    }

                }

                ot.setTimePoint("KosDot.send.transposeDone");

                seed = prng.get<block>();

                CP_SEND(seed);
                //chl.asyncSend((u8*)&seed, sizeof(block));

                codePrng.SetSeed(seed);
                code.random(codePrng, ot.mBaseChoiceBits.size(), 128);
                code.encode((u8*)delta.data(), (u8*)&curDelta);

                if (eq(ot.mDelta, ZeroBlock))
                    ot.mDelta = prng.get<block>();

                offset = curDelta ^ ot.mDelta;
                CP_SEND(offset);
                //chl.asyncSend(offset);
                CP_RECV(theirSeed);
                //chl.recv((u8*)&theirSeed, sizeof(block));
                ot.setTimePoint("KosDot.send.cncSeed");

                if (Commit(theirSeed) != theirSeedComm)
                    throw std::runtime_error("bad commit " LOCATION);

                commonPrng.SetSeed(seed ^ theirSeed);
                memset(q.data(), 0, sizeof(q));
                ot.setTimePoint("KosDot.send.checkStart");

                xx = 0;
                end = (messages.size() + 127 + 128) / 128;
                for (u64 blockIdx = 0; blockIdx < end; ++blockIdx)
                {
                    commonPrng.mAes.ecbEncCounterMode(doneIdx, 128, challenges.data());
                    commonPrng.mAes.ecbEncCounterMode(doneIdx, 128, challenges2.data());

                    u64 stop0 = std::min<u64>(messages.size(), doneIdx + 128);
                    u64 stop1 = std::min<u64>(messages.size() + 128, doneIdx + 128);

                    u64 i = 0, dd = doneIdx;
                    for (; dd < stop0; ++dd, ++i)
                    {
                        mul256(messages[dd][0], messages[dd][1],
                            challenges[i], challenges2[i],
                            qi[0], qi[1], qi[2], qi[3]);

                        q[0] = q[0] ^ qi[0];
                        q[1] = q[1] ^ qi[1];
                        q[2] = q[2] ^ qi[2];
                        q[3] = q[3] ^ qi[3];


                        code.encode((u8*)messages[dd].data(), (u8*)&messages[dd][0]);
                        messages[dd][1] = messages[dd][0] ^ ot.mDelta;
                        //code.encode((u8*)messages1.data(), (u8*)&messages[dd][1]);
                    }


                    for (; dd < stop1; ++dd, ++i, ++xx)
                    {
                        mul256(extraBlocks[xx][0], extraBlocks[xx][1], challenges[i], challenges2[i],
                            qi[0], qi[1], qi[2], qi[3]);

                        q[0] = q[0] ^ qi[0];
                        q[1] = q[1] ^ qi[1];
                        q[2] = q[2] ^ qi[2];
                        q[3] = q[3] ^ qi[3];
                    }

                    doneIdx = stop1;
                }

                ot.setTimePoint("KosDot.send.checkSummed");
                proof.resize(8);
                CP_RECV(proof);
                ot.setTimePoint("KosDot.send.proofReceived");
                {

                    auto& received_x = ((std::array<block, 4>*)proof.data())[0];
                    auto& received_t = ((std::array<block, 4>*)proof.data())[1];

                    // check t = x * Delta + q
                    mul256(received_x[0], received_x[1],
                        delta[0], delta[1],
                        tt[0], tt[1], tt[2], tt[3]);

                    tt[0] = tt[0] ^ q[0];
                    tt[1] = tt[1] ^ q[1];
                    tt[2] = tt[2] ^ q[2];
                    tt[3] = tt[3] ^ q[3];

                    if (eq(tt[0], received_t[0]) && eq(tt[1], received_t[1]) &&
                        eq(tt[2], received_t[2]) && eq(tt[3], received_t[3]))
                    {
                        //std::cout << "\tCheck passed\n";
                    }
                    else
                    {
                        std::cout << "OT Ext Failed Correlation check failed" << std::endl;
                        std::cout << "rec t[0] = " << received_t[0] << std::endl;
                        std::cout << "rec t[1] = " << received_t[1] << std::endl;
                        std::cout << "rec t[2] = " << received_t[2] << std::endl;
                        std::cout << "rec t[3] = " << received_t[3] << std::endl << std::endl;
                        std::cout << "exp t[0] = " << tt[0] << std::endl;
                        std::cout << "exp t[1] = " << tt[1] << std::endl;
                        std::cout << "exp t[2] = " << tt[2] << std::endl;
                        std::cout << "exp t[3] = " << tt[3] << std::endl << std::endl;
                        std::cout << "q  = " << q[0] << std::endl;
                        throw std::runtime_error("Exit");;
                    }

                    ot.setTimePoint("KosDot.send.done");

                    static_assert(gOtExtBaseOtCount == 128, "expecting 128");

                }
                CP_END();
                return {};
            }
        };

        return coproto::makeProto<Proto>(*this, messages, prng);

    }


}
#endif