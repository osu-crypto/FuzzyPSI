#include "IknpDotExtReceiver.h"
#ifdef ENABLE_DELTA_IKNP
#include "libOTe/Tools/Tools.h"

#include <cryptoTools/Common/BitVector.h>
#include <cryptoTools/Common/Matrix.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Crypto/PRNG.h>
#include <cryptoTools/Crypto/Commit.h>
#include <cryptoTools/Network/Channel.h>

#include "TcoOtDefines.h"
#include <queue>


namespace osuCrypto
{

    IknpDotExtReceiver IknpDotExtReceiver::splitBase()
    {
        std::vector<std::array<block, 2>>baseRecvOts(mGens.size());

        for (u64 i = 0; i < mGens.size(); ++i)
        {
            baseRecvOts[i][0] = mGens[i][0].get<block>();
            baseRecvOts[i][1] = mGens[i][1].get<block>();
        }

        return IknpDotExtReceiver(baseRecvOts);
    }

    std::unique_ptr<OtExtReceiver> IknpDotExtReceiver::split()
    {
        std::vector<std::array<block, 2>>baseRecvOts(mGens.size());

        for (u64 i = 0; i < mGens.size(); ++i)
        {
            baseRecvOts[i][0] = mGens[i][0].get<block>();
            baseRecvOts[i][1] = mGens[i][1].get<block>();
        }

        return std::make_unique<IknpDotExtReceiver>(baseRecvOts);
    }


}
#endif