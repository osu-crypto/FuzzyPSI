#include <iostream>
#include <cryptoTools/Common/Defines.h>
#include <cryptoTools/Network/Channel.h>
#include <cryptoTools/Network/Session.h>
#include <cryptoTools/Network/IOService.h>
#include <cryptoTools/Common/Timer.h>
#include <cryptoTools/Crypto/AES.h>
#include <cryptoTools/Crypto/PRNG.h>
#include "frontend/libPaXoS/ObliviousDictionary.h"
#include <cryptoTools/Common/BitVector.h>
#include <string>

//for OKVSDecode
#include <NTL/mat_GF2E.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>


void OKVSDecode(vector<osuCrypto::block> okvs, osuCrypto::block key);
void PaxosEncode(const std::vector<osuCrypto::block> setKeys, const std::vector<osuCrypto::block> setValues, std::vector<osuCrypto::block>& okvs, uint64_t fieldSize);
void ShareFss_farapart(uint64_t delta, int nSquares);