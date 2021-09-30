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
#include "cryptoTools/Common/BitIterator.h"
#include <cryptoTools/Crypto/RandomOracle.h>
#include <string>
#include <unordered_map>

//for OKVSDecode using blocks
#include <NTL/mat_GF2E.h>
#include <NTL/GF2E.h>
#include <NTL/GF2X.h>
#include <NTL/GF2XFactoring.h>

//void OKVSDecode(vector<osuCrypto::block> okvs, osuCrypto::block key);
//for arbitrary length val in terms of blocks
void OkvsEncode(std::vector<uint64_t> setKeys, const std::vector<array<osuCrypto::block, 2>> setValues0, std::vector<array<osuCrypto::block, 12>>& okvs0, uint64_t fieldSize);
//for a fixed length block value
void PaxosEncode(const std::vector<uint64_t> setKeys, const std::vector<osuCrypto::block> setValues0, const std::vector<osuCrypto::block> setValues1, std::vector<osuCrypto::block>& okvs0, std::vector<osuCrypto::block>& okvs1, uint64_t fieldSize);
void far_apart_FssShare(uint64_t delta, int nSquares, std::vector<osuCrypto::block>& okvs0, std::vector<osuCrypto::block>& okvs1);
void far_apart_FssEval(uint64_t x_coord, uint64_t y_coord, vector<osuCrypto::block> okvs, uint64_t delta, uint64_t hashSize);
//void FssShare_Enumerate(uint64_t delta, int nSquares, vector<vector<osuCrypto::block>> &okvs0, vector<vector<osuCrypto::block>> &okvs1);