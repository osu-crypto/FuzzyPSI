#include "naor-pinkas.h"

#include <cryptoTools/Common/Log.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Crypto/RandomOracle.h>
#include <cryptoTools/Network/Channel.h>

#include <cryptoTools/Crypto/RCurve.h>
#include <cryptoTools/Crypto/Curve.h>

#include "libOTe/Tools/CoprotoSock.h"
#include "coproto/NativeProto.h"
#include "coproto/Macros.h"

#include <iomanip>
#ifdef ENABLE_NP

#include <memory>

namespace osuCrypto
{

#ifdef ENABLE_RELIC
    using Curve = REllipticCurve;
    using Point = REccPoint;
    using Brick = REccPoint;
    using Number = REccNumber;
#else    
    using Curve = EllipticCurve;
    using Point = EccPoint;
    using Brick = EccBrick;
    using Number = EccNumber;
#endif
    //static const  u64 minMsgPerThread(16);


    //std::string hexStr(span<u8> data)
    //{
    //    std::stringstream ss;
    //    auto min = std::min<u64>(data.size(), 32);
    //    for (u64 i = 0; i < min; ++i)
    //    {
    //        ss << std::hex << std::setw(2) << std::setfill('0') << int(data[i]);
    //    }

    //    if (min != data.size())
    //        ss << "...";

    //    return ss.str();
    //}

    //template<typename T>
    //typename std::enable_if<std::is_trivial<T>::value,
    //    std::string>::type hexStr(const T& data)
    //{
    //    span<u8> dd((u8*)&data, sizeof(data));
    //    return hexStr(dd);
    //}



    coproto::Proto NaorPinkas::receive(
        const BitVector& choices,
        span<block> messages,
        PRNG& prng)
    {

        struct NPProto : public coproto::NativeProto
        {
            const BitVector& choices;
            span<block> messages;
            PRNG& prng;
            NPProto(const BitVector& c, span<block> m, PRNG& p)
                : choices(c)
                , messages(m)
                , prng(p)
                , g(_curve)
                , PK0(_curve)
                , bg(_curve)
            {}

            Curve _curve;
            u64 nSndVals = 2;
            Point g;
            u64 fieldElementSize;
            //std::vector<std::thread> thrds;// (numThreads);
            std::vector<u8> buff;// (messages.size()* fieldElementSize);
            std::vector<u8> cBuff;//(nSndVals* fieldElementSize);
            std::array<u8, RandomOracle::HashSize> comm, comm2;
            block R;

            Point PK0;
            Brick bg;

            std::vector<Number> pK;
            std::vector<Point>
                PK_sigma,
                pC;

            std::vector<u8>::iterator pBufIdx;
            u8* iter;
            RandomOracle ro;


            coproto::error_code resume() override
            {
                Curve curve;
                //auto mStart = 0;// t* messages.size() / numThreads;
                //auto mEnd = messages.size();// (t + 1)* messages.size() / numThreads;

                CP_BEGIN();
                g = curve.getGenerator();
                fieldElementSize = g.sizeBytes();

                //thrds.resize(numThreads);
                buff.resize(messages.size() * fieldElementSize);

                pK.reserve(messages.size());
                PK_sigma.reserve(messages.size());
                pC.reserve(nSndVals);
                bg = g;

                for (u64 i =  0; i < (u64)messages.size(); ++i)
                {
                    // get a random value from Z_p
                    pK.emplace_back(curve);
                    pK[i].randomize(prng);

                    // using brickexp which has the base of g, compute
                    //
                    //      PK_sigma[i] = g ^ pK[i]
                    //
                    // where pK[i] is just a random number in Z_p
                    PK_sigma.emplace_back(curve);
                    PK_sigma[i] = bg * pK[i];

                    //if (i < 4)
                    //    std::cout << "pK " << i << " " << pK[j] << std::endl;
                }

                cBuff.resize(nSndVals * fieldElementSize);

                CP_RECV(cBuff);
                //std::cout << "p0 recv " << cBuff.size() << " " << hexStr(cBuff) << std::endl;
                //cRecvFuture.get();
                pBufIdx = cBuff.begin();
                for (u64 u = 0; u < nSndVals; u++)
                {
                    pC.emplace_back(curve);

                    pC[u].fromBytes(&*pBufIdx);
                    pBufIdx += fieldElementSize;
                }

                iter = buff.data();

                for (u64 i = 0; i < (u64)messages.size(); ++i)
                {
                    u8 choice = choices[i];
                    if (choice != 0) {
                        PK0 = pC[choice] - PK_sigma[i];
                    }
                    else {
                        PK0 = PK_sigma[i];
                    }

                    //if(i < 4)
                    //    std::cout << "PK0 " << i << " " << PK0 << std::endl;
                    PK0.toBytes(iter);
                    iter += fieldElementSize;
                }

                //if (--remainingPK0s == 0)
                //    PK0Prom.set_value();
                CP_RECV(comm);
                //std::cout << "p0 recv comm " << hexStr(comm) << std::endl;
                //std::cout << "p0 send " << buff.size() << " " << hexStr(buff) << std::endl;
                CP_SEND(std::move(buff));

                // resuse this space, not the data of PK0...

                CP_END_OF_ROUND();

                buff.resize(fieldElementSize);
                CP_RECV(R);
                //std::cout << "p0 recv R " << hexStr(R) << std::endl;
                ro.Update(R);
                ro.Final(comm2);
                if (comm != comm2)
                    throw std::runtime_error("bad commitment " LOCATION);

                ro = RandomOracle(sizeof(block));
                for (u64 i = 0; i < (u64)messages.size(); ++i)
                {
                    // now compute g ^(a * k) = (g^a)^k
                    PK0 = pC[0] * pK[i];
                    PK0.toBytes(buff.data());

                    ro.Reset();
                    ro.Update((u8*)&i, sizeof(i));
                    ro.Update(buff.data(), buff.size());
                    ro.Update(R);
                    ro.Final(messages[i]);
                }

                CP_END();
                return {};
            }
        };
        return coproto::makeProto<NPProto>(choices, messages, prng);
    }

    coproto::Proto NaorPinkas::send(
        span<std::array<block, 2>> messages,
        PRNG& prng)
    {

        struct NPProto : public coproto::NativeProto
        {

            span<std::array<block, 2>> messages;
            PRNG& prng;

            NPProto(
                span<std::array<block, 2>> m,
                PRNG& p)
                : messages(m)
                , prng(p)
                , alpha(_curve)
                , tmp(_curve)
                , pPK0(_curve)
                , PK0a(_curve)
                , fetmp(_curve)
            {}

            block R;
            u64 nSndVals = 2;
            Curve _curve;
            Number alpha;//(curve, prng), tmp(curve);
            Number tmp;
            Point g;
            u64 fieldElementSize;
            std::vector<u8> buff;
            std::vector<Point> pC;
            RandomOracle ro;
            std::vector<u8> hashBuff;// (RandomOracle::HashSize);
            Point pPK0, PK0a, fetmp;
            std::vector<Point> c;

            coproto::error_code resume() override
            {
                Curve curve;
                //auto mStart = 0;
                //auto mEnd = messages.size();

                CP_BEGIN();
                R = prng.get<block>();

                alpha.randomize(prng);
                g = curve.getGenerator();
                fieldElementSize = g.sizeBytes();
                buff.resize(nSndVals * fieldElementSize);
                pC.reserve(nSndVals);


                pC.emplace_back(curve);
                pC[0] = g * alpha;
                pC[0].toBytes(buff.data());

                for (u64 u = 1; u < nSndVals; u++)
                {
                    pC.emplace_back(curve);
                    tmp.randomize(prng);

                    pC[u] = g * tmp;
                    pC[u].toBytes(buff.data() + u * fieldElementSize);
                }

                //std::cout << "p1 send " << buff.size() << " " << hexStr(buff) << std::endl;
                CP_SEND(std::move(buff));

                // sends a commitment to R. This strengthens the security of NP01 to 
                // make the protocol output uniform strings no matter what.
                hashBuff.resize(RandomOracle::HashSize);
                ro.Update(R);
                ro.Final(hashBuff.data());
                //std::cout << "p1 send comm " << hexStr(hashBuff) << std::endl;
                CP_SEND(std::move(hashBuff));


                for (u64 u = 1; u < nSndVals; u++)
                    pC[u] = pC[u] * alpha;

                c.reserve(nSndVals);
                for (u64 i = 0; i < nSndVals; ++i)
                    c.emplace_back(curve, pC[i]);

                hashBuff.resize(fieldElementSize);


                buff.resize(fieldElementSize * messages.size());
                CP_RECV(buff);
                //std::cout << "p1 recv " << buff.size() << " " << hexStr(buff) << std::endl;

                //std::cout << "p1 send R " << hexStr(R) << std::endl;
                CP_SEND(R);


                ro = RandomOracle(sizeof(block));

                for (u64 i = 0; i < (u64)messages.size(); i++)
                {

                    pPK0.fromBytes(buff.data() + i * fieldElementSize);
                    PK0a = pPK0 * alpha;
                    PK0a.toBytes(hashBuff.data());

                    ro.Reset();
                    ro.Update((u8*)&i, sizeof(i));
                    ro.Update(hashBuff.data(), hashBuff.size());
                    ro.Update(R);
                    ro.Final(messages[i][0]);

                    for (u64 u = 1; u < nSndVals; u++)
                    {
                        fetmp = c[u] - PK0a;
                        fetmp.toBytes(hashBuff.data());

                        ro.Reset();
                        ro.Update((u8*)&i, sizeof(i));
                        ro.Update(hashBuff.data(), hashBuff.size());
                        ro.Update(R);
                        ro.Final(messages[i][u]);
                    }
                }

                CP_END();

                return {};
            }
        };

        return coproto::makeProto<NPProto>(messages, prng);
    }


    void NaorPinkas::receive(
        const BitVector& choices,
        span<block> messages,
        PRNG& prng,
        Channel& socket,
        u64 numThreads)
    {
        CoprotoSock s(socket);
        receive(choices, messages, prng).evaluate(s);
    }
    void NaorPinkas::send(span<std::array<block, 2>> messages, PRNG& prng, Channel& sock, u64 numThreads)
    {

        CoprotoSock s(sock);
        send(messages, prng).evaluate(s);
    }



    //void NaorPinkas::receive(
    //    const BitVector& choices,
    //    span<block> messages,
    //    PRNG& prng,
    //    Channel& socket,
    //    u64 numThreads)
    //{
    //    // should generalize to 1 out of N by changing this. But isn't tested...
    //    auto nSndVals(2);
    //    Curve curve;
    //    auto g = curve.getGenerator();
    //    u64 fieldElementSize = g.sizeBytes();

    //    std::vector<std::thread> thrds(numThreads);
    //    std::vector<u8> sendBuff(messages.size() * fieldElementSize);
    //    //std::atomic<u32> remainingPK0s((u32)numThreads);
    //    //std::promise<void> PK0Prom;
    //    //std::future<void> PK0Furture(PK0Prom.get_future());
    //    std::vector<u8> cBuff(nSndVals * fieldElementSize);
    //    //auto cRecvFuture = socket.asyncRecv(cBuff.data(), cBuff.size()).share();
    //    block R;

    //    std::array<u8, RandomOracle::HashSize> comm, comm2;
    //    //socket.asyncRecv(comm);
    //    //auto RFuture = socket.asyncRecv(R).share();

    //    //for (u64 t = 0; t < numThreads; ++t)
    //    //{
    //    //    auto seed = prng.get<block>();

    //    //    thrds[t] = std::thread(
    //    //        [t, numThreads, &messages, seed,
    //    //        &sendBuff, &choices, cRecvFuture, &cBuff,
    //    //        &remainingPK0s, &PK0Prom, nSndVals, &RFuture, &R, &comm, &prng]()
    //    //        {
    //    auto t = 0;
    //    auto mStart = t * messages.size() / numThreads;
    //    auto mEnd = (t + 1) * messages.size() / numThreads;

    //    //PRNG prng(seed);
    //    //Curve curve;

    //    //auto g = curve.getGenerator();
    //    //u64 fieldElementSize = g.sizeBytes();

    //    Point PK0(curve);
    //    Brick bg(g);

    //    std::vector<Number> pK;
    //    std::vector<Point>
    //        PK_sigma,
    //        pC;

    //    pK.reserve(mEnd - mStart);
    //    PK_sigma.reserve(mEnd - mStart);
    //    pC.reserve(nSndVals);

    //    for (u64 i = mStart, j = 0; i < mEnd; ++i, ++j)
    //    {
    //        // get a random value from Z_p
    //        pK.emplace_back(curve);
    //        pK[j].randomize(prng);

    //        // using brickexp which has the base of g, compute
    //        //
    //        //      PK_sigma[i] = g ^ pK[i]
    //        //
    //        // where pK[i] is just a random number in Z_p
    //        PK_sigma.emplace_back(curve);
    //        PK_sigma[j] = bg * pK[j];


    //        if (i < 4)
    //            std::cout << "pK " << i << " " << pK[j] << std::endl;

    //    }

    //    socket.recv(cBuff.data(), cBuff.size());

    //    std::cout << "p0 recv " << cBuff.size() << " " << hexStr(cBuff) << std::endl;

    //    auto pBufIdx = cBuff.begin();
    //    for (auto u = 0; u < nSndVals; u++)
    //    {
    //        pC.emplace_back(curve);

    //        pC[u].fromBytes(&*pBufIdx);
    //        pBufIdx += fieldElementSize;
    //    }

    //    auto iter = sendBuff.data() + mStart * fieldElementSize;

    //    for (u64 i = mStart, j = 0; i < mEnd; ++i, ++j)
    //    {
    //        u8 choice = choices[i];
    //        if (choice != 0) {
    //            PK0 = pC[choice] - PK_sigma[j];
    //        }
    //        else {
    //            PK0 = PK_sigma[j];
    //        }

    //        if (i < 4)
    //            std::cout << "PK0 " << i << " " << PK0 << std::endl;

    //        PK0.toBytes(iter);
    //        iter += fieldElementSize;
    //    }

    //    //if (--remainingPK0s == 0)
    //    {

    //        socket.recv(comm);
    //        std::cout << "p0 send " << sendBuff.size() << " " << hexStr(sendBuff) << std::endl;
    //        socket.asyncSend(std::move(sendBuff));
    //        //PK0Prom.set_value();
    //    }

    //    // resuse this space, not the data of PK0...
    //    auto& gka = PK0;
    //    RandomOracle sha(sizeof(block));

    //    std::vector<u8>buff(fieldElementSize);
    //    Brick bc(pC[0]);

    //    //RFuture.get();
    //    socket.recv(R);
    //    std::cout << "p0 recv comm " << hexStr(comm) << std::endl;
    //    std::cout << "p0 recv R " << hexStr(R) << std::endl;


    //    for (u64 i = mStart, j = 0; i < mEnd; ++i, ++j)
    //    {
    //        // now compute g ^(a * k) = (g^a)^k
    //        gka = bc * pK[j];
    //        gka.toBytes(buff.data());

    //        sha.Reset();
    //        sha.Update((u8*)&i, sizeof(i));
    //        sha.Update(buff.data(), buff.size());
    //        sha.Update(R);
    //        sha.Final(messages[i]);
    //    }
    //    //        });
    //    //}

    //    //PK0Furture.get();



    //    //for (auto& thrd : thrds)
    //    //    thrd.join();

    //    //block comm = *(block*)(cBuff.data() + nSndVals * fieldElementSize);
    //    RandomOracle ro;
    //    ro.Update(R);
    //    ro.Final(comm2);
    //    if (comm != comm2)
    //        throw std::runtime_error("bad commitment " LOCATION);

    //}


    //void NaorPinkas::send(
    //    span<std::array<block, 2>> messages,
    //    PRNG& prng,
    //    Channel& socket,
    //    u64 numThreads)
    //{

    //    block R = prng.get<block>();
    //    // one out of nSndVals OT.
    //    u64 nSndVals(2);
    //    std::vector<std::thread> thrds(numThreads);
    //    //auto seed = prng.get<block>();
    //    Curve curve;
    //    Number alpha(curve, prng), tmp(curve);
    //    const Point g = curve.getGenerator();
    //    u64 fieldElementSize = g.sizeBytes();
    //    std::vector<u8> sendBuff(nSndVals * fieldElementSize);
    //    std::vector<Point> pC;
    //    pC.reserve(nSndVals);


    //    pC.emplace_back(curve);
    //    pC[0] = g * alpha;
    //    pC[0].toBytes(sendBuff.data());

    //    for (u64 u = 1; u < nSndVals; u++)
    //    {
    //        pC.emplace_back(curve);
    //        tmp.randomize(prng);

    //        pC[u] = g * tmp;
    //        pC[u].toBytes(sendBuff.data() + u * fieldElementSize);
    //    }

    //    std::cout << "p1 send " << sendBuff.size() << " " << hexStr(sendBuff) << std::endl;
    //    socket.asyncSend(std::move(sendBuff));

    //    // sends a commitment to R. This strengthens the security of NP01 to 
    //    // make the protocol output uniform strings no matter what.
    //    RandomOracle ro;
    //    std::vector<u8> comm(RandomOracle::HashSize);
    //    ro.Update(R);
    //    ro.Final(comm.data());
    //    std::cout << "p1 send comm " << hexStr(comm) << std::endl;
    //    socket.asyncSend(std::move(comm));


    //    for (u64 u = 1; u < nSndVals; u++)
    //        pC[u] = pC[u] * alpha;

    //    std::vector<u8> buff(fieldElementSize * messages.size());
    //    //auto recvFuture = socket.asyncRecv(buff.data(), buff.size()).share();

    //    //for (u64 t = 0; t < numThreads; ++t)
    //    //{

    //    //    thrds[t] = std::thread([
    //    //        t, seed, fieldElementSize, &messages, recvFuture,
    //    //            numThreads, &buff, &alpha, nSndVals, &pC, &socket, &R]()
    //    //        {
    //    //            Curve curve;

    //    auto t = 0;
    //    Point pPK0(curve), PK0a(curve), fetmp(curve);
    //    Number alpha2(curve, alpha);

    //    std::vector<Point> c;
    //    c.reserve(nSndVals);
    //    for (u64 i = 0; i < nSndVals; ++i)
    //        c.emplace_back(curve, pC[i]);

    //    std::vector<u8> hashInBuff(fieldElementSize);
    //    RandomOracle sha(sizeof(block));
    //    //recvFuture.get();
    //    socket.recv(buff.data(), buff.size());
    //    std::cout << "p1 recv " << buff.size() << " " << hexStr(buff) << std::endl;

    //    if (t == 0)
    //    {
    //        std::cout << "p1 send R " << hexStr(R) << std::endl;
    //        socket.asyncSendCopy(R);
    //    }


    //    auto mStart = t * messages.size() / numThreads;
    //    auto mEnd = (t + 1) * messages.size() / numThreads;

    //    for (u64 i = mStart; i < mEnd; i++)
    //    {

    //        pPK0.fromBytes(buff.data() + i * fieldElementSize);
    //        PK0a = pPK0 * alpha2;
    //        PK0a.toBytes(hashInBuff.data());

    //        sha.Reset();
    //        sha.Update((u8*)&i, sizeof(i));
    //        sha.Update(hashInBuff.data(), hashInBuff.size());
    //        sha.Update(R);
    //        sha.Final(messages[i][0]);

    //        for (u64 u = 1; u < nSndVals; u++)
    //        {
    //            fetmp = c[u] - PK0a;
    //            fetmp.toBytes(hashInBuff.data());

    //            sha.Reset();
    //            sha.Update((u8*)&i, sizeof(i));
    //            sha.Update(hashInBuff.data(), hashInBuff.size());
    //            sha.Update(R);
    //            sha.Final(messages[i][u]);
    //        }
    //    }
    //    //        });
    //    //}

    //    //for (auto& thrd : thrds)
    //    //    thrd.join();

    //}

}

#endif
