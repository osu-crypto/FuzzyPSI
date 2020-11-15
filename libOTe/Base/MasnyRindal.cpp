#include "MasnyRindal.h"

#ifdef ENABLE_MR

#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Log.h>
#include <cryptoTools/Crypto/RandomOracle.h>
#include <cryptoTools/Network/Channel.h>

#include <cryptoTools/Crypto/RCurve.h>
#ifndef ENABLE_RELIC
static_assert(0, "ENABLE_RELIC must be defined to build MasnyRindal");
#endif

using Curve = oc::REllipticCurve;
using Point = oc::REccPoint;
using Brick = oc::REccPoint;
using Number = oc::REccNumber;

#include <libOTe/Base/SimplestOT.h>
#include "coproto/Buffers.h"
#include "coproto/NativeProto.h"
#include "coproto/AwaitProto.h"
#include "coproto/Macros.h"

namespace osuCrypto
{
	const u64 step = 16;

	coproto::Proto MasnyRindal::receive(const BitVector& choices, span<block> messages, PRNG& prng)
	{
		struct Frame : coproto::NativeProto
		{

			Frame(const BitVector& c, span<block>& m, PRNG& p)
				:choices(c)
				, messages(m)
				, prng(p)
			{}
			const BitVector& choices;
			span<block> messages;
			PRNG& prng;

			u64 n;
			Curve _curve;
			std::array<Point, 2> r;
			Point g;
			u64 pointSize;
			Point hPoint;
			std::vector<u8> hashBuff;
			std::vector<u8> ioBuff;
			std::vector<Number> sk;
			Point Mb, k;
			u64 i, curStep;
			u8* iter;

			coproto::error_code resume() override
			{
				Curve curve;

				CP_BEGIN();
				n = choices.size();
				g = curve.getGenerator();
				pointSize = g.sizeBytes();

				hashBuff.resize(roundUpTo(pointSize, 16));
				sk.reserve(n);


				for (i = 0; i < n;)
				{
					curStep = std::min<u64>(n - i, step);

					ioBuff.resize(pointSize * 2 * curStep);
					iter = ioBuff.data();


					for (u64 k = 0; k < curStep; ++k, ++i)
					{

						auto& rrNot = r[choices[i] ^ 1];
						auto& rr = r[choices[i]];

						rrNot.randomize();
						rrNot.toBytes(hashBuff.data());

						ep_map(hPoint, hashBuff.data(), int(pointSize));

						sk.emplace_back(prng);

						rr = g * sk[i];
						rr -= hPoint;

						r[0].toBytes(iter); iter += pointSize;
						r[1].toBytes(iter); iter += pointSize;
					}

					CP_SEND(std::move(ioBuff));
					//co_await coproto::send(std::move(ioBuff));
				}


				ioBuff.resize(pointSize);

				CP_RECV(ioBuff);
				Mb.fromBytes(ioBuff.data());

				for (u64 i = 0; i < n; ++i)
				{
					k = Mb;
					k *= sk[i];

					//lout << "g^ab  " << k << std::endl;

					k.toBytes(hashBuff.data());

					RandomOracle ro(sizeof(block));
					ro.Update(hashBuff.data(), pointSize);
					ro.Update(i);
					ro.Final(messages[i]);
				}

				CP_END();

				return {};
			}
		};

		return coproto::makeProto<Frame>(choices, messages, prng);

	}
	//struct Controller
	//{
	//	error_code send(span<u8>);
	//	error_code recv(span<u8>);

	//};


	coproto::Proto MasnyRindal::send(span<std::array<block, 2>> messages, PRNG& prng)
	{
		struct Frame : public coproto::NativeProto
		{
			Frame(span<std::array<block, 2>>& m, PRNG& p)
				:messages(m)
				, prng(p)
			{}

			span<std::array<block, 2>> messages;
			PRNG& prng;

			Curve _curve;
			u64 n;
			Point g;
			u64 pointSize;
			Number sk;
			Point Mb;
			std::vector<u8> buff, hashBuff;
			Point pHash, r;
			u64 i, curStep;

			coproto::error_code resume() override
			{
				// this initializes relic if its not already done.
				Curve curve;

				CP_BEGIN();

				n = static_cast<u64>(messages.size());

				g = curve.getGenerator();
				pointSize = g.sizeBytes();
				sk = Number(prng);
				Mb = g * sk;

				buff.resize(pointSize);
				hashBuff.resize(roundUpTo(pointSize, 16));

				Mb.toBytes(buff.data());


				CP_SEND(std::move(buff));

				for (i = 0; i < n; )
				{
					curStep = std::min<u64>(n - i, step);
					buff.resize(curStep * pointSize * 2);

					CP_RECV(buff);

					auto buffIter = buff.data();


					for (u64 k = 0; k < curStep; ++k, ++i)
					{
						std::array<u8*, 2> buffIters{
							buffIter,
							buffIter + pointSize
						};
						buffIter += pointSize * 2;

						for (u64 j = 0; j < 2; ++j)
						{
							r.fromBytes(buffIters[j]);
							ep_map(pHash, buffIters[j ^ 1], int(pointSize));

							r += pHash;
							r *= sk;

							r.toBytes(hashBuff.data());

							RandomOracle ro(sizeof(block));
							ro.Update(hashBuff.data(), pointSize);
							ro.Update(i);
							ro.Final(messages[i][j]);

						}
					}

				}

				CP_END();
				return {};
			}
		};


		return coproto::makeProto<Frame>(messages, prng);

		//for (u64 i = 0; i < n; )
		//{
		//    auto curStep = std::min<u64>(n - i, step);
		//    buff.resize(curStep * pointSize * 2);

		//    co_await coproto::recvFixedSize(buff);
		//    auto buffIter = buff.data();


		//    for (u64 k = 0; k < curStep; ++k, ++i)
		//    {
		//        std::array<u8*, 2> buffIters{
		//            buffIter,
		//            buffIter + pointSize
		//        };
		//        buffIter += pointSize * 2;

		//        for (u64 j = 0; j < 2; ++j)
		//        {
		//            r.fromBytes(buffIters[j]);
		//            ep_map(pHash, buffIters[j ^ 1], int(pointSize));

		//            r += pHash;
		//            r *= sk;

		//            r.toBytes(hashBuff.data());
		//            auto p = (block*)hashBuff.data();

		//            ro.Reset();
		//            ro.Update(hashBuff.data(), pointSize);
		//            ro.Update(i);
		//            ro.Final(messages[i][j]);
		//        }
		//    }
		//}
	}

	//void MasnyRindal::receive(
	//    const BitVector & choices,
	//    span<block> messages,
	//    PRNG & prng,
	//    Channel & chl)
	//{
	//    throw RTE_LOC;
	//}

	//void MasnyRindal::send(span<std::array<block, 2>> messages, PRNG & prng, Channel & chl)
	//{
	//    throw RTE_LOC;
	//}
}
#endif