#include "okvsfss.h"

using namespace osuCrypto;

/* How we write a square and identify the grid cells
    1. Grid cell label is bottom left co-coordinate
    2. We can identify a square also by bottom left coord

    Note: each square is processed independent of the other
    this allows us to sample the input squares systematically
    For now, we will assume that all the squares are in a straight line

    How to generate input squares
    1. Sample the bottom left square/center of the square
    2. Sample the side length = 2 * delta
    3. First square Bottomleftlabel - (delta, delta)
    4. To get next square add (2*delta) along the the x-axis
*/

void return_grid(uint64_t x, uint64_t y, uint64_t &grd_x, uint64_t &grd_y, uint64_t sqrlen){
        grd_x = x / sqrlen;
        grd_y = y / sqrlen;
}

void PaxosEncode(const std::vector<block> setKeys, const std::vector<block> setValues0, const std::vector<block> setValues1, std::vector<block>& okvs0, std::vector<block>& okvs1, uint64_t fieldSize)
{
    GF2E dhBitsVal;
    int hashSize=setKeys.size(), gamma = 60, v=20;
    double c1 = 2.4;
    vector<uint64_t> keys;
    keys.resize(hashSize);
    int fieldSizeBytes = fieldSize % 8 == 0 ? fieldSize/8 : fieldSize/8 + 1;
    int zeroBits = 8 - fieldSize % 8;

    //initialize 
    ObliviousDictionary * dict0 = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    ObliviousDictionary * dict1 = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict0->init();
    dict1->init();

    // hashing keys
    for (int i=0; i < setKeys.size(); i++){
        keys[i] = XXH64(&setKeys[i], 64, 0);
    }

    vector<byte> values0, values1;
    values0.resize(setValues0.size() * sizeof(block));
    values1.resize(setValues1.size() * sizeof(block));
    memcpy(values0.data(), setValues0.data(), setValues0.size() * sizeof(block));
    memcpy(values1.data(), setValues1.data(), setValues1.size() * sizeof(block));
    
    dict0->setKeysAndVals(keys, values0);
    dict1->setKeysAndVals(keys, values1);
    dict0->encode();
    dict1->encode();

    vector<byte> x0 = dict0->getVariables();
    vector<byte> x1 = dict1->getVariables();
    
    // vector of **bytes to block** , vector<block> Okvs;
    okvs0.resize(x0.size() / fieldSizeBytes);
    okvs1.resize(x1.size() / fieldSizeBytes);
    memcpy(okvs0.data(), x0.data(), x0.size());
    memcpy(okvs1.data(), x1.data(), x1.size());
    /*std::cout << "Okvs.size() " << Okvs.size() << std::endl;
    for (int i = 0; i < (Okvs.size() - 50); i++)
        std::cout << "okvs  " << i << "  " << Okvs[i] << std::endl;*/
    
    //below is function call to check output!!
    dict0->checkOutput();
    dict1->checkOutput();
    
    auto indices = dict0->dec(keys[0]);
    cout << "ENCODE " << setKeys[0] << " " << keys[0] << " " << indices[0] << " " << indices[1] << std::endl;
   
} 

/*void OKVSDecode(vector<block> okvs, block key){
    cout << "Let's write decode " << std::endl;
    // 1. convert block key -> to unint64_t keys using xxH64
    // nothing changes in setKeysandValues()
  // converting block to bytes
    vector<byte> okvs_bytes;
    okvs_bytes.resize(Okvs.size() * sizeof(block));
    memcpy(okvs_bytes.data(), Okvs.data(), Okvs.size() * sizeof(block));

    // from bytes to GF2E!! yahooo 
    GF2EVector okvs_gf2e(Okvs.size());
    GF2X temp;
    for (int i=0; i < okvs_gf2e.size(); i++){
        GF2XFromBytes(temp, okvs_bytes.data() + i*fieldSizeBytes ,fieldSizeBytes);
        okvs_gf2e[i] = to_GF2E(temp); 
        //std::cout << okvs_gf2e[i] << std::endl;
    }
    
    // decoding a single key
    std::cout << "setKey[4] " << setKeys[4] << std::endl;

    // decoding from checkoutput()
    GF2E dhBitsVals;
    auto indices = dic->dec(keys[5]);
    for (int j=0; j<indices.size(); j++){
        dhBitsVals += okvs_gf2e[indices[j]];
    }
    //std::cout << "dhBitsVals  " << dhBitsVals << std::endl;
    vector<byte> valBytes(fieldSizeBytes);
    BytesFromGF2X(valBytes.data(), rep(dhBitsVals), fieldSizeBytes);
    block returndecode(valBytes[15],valBytes[14],valBytes[13],valBytes[12],valBytes[11],valBytes[10],valBytes[9],valBytes[8],valBytes[7],valBytes[6],valBytes[5],valBytes[4],valBytes[3],valBytes[2],valBytes[1],valBytes[0]);
    std::cout << "setVals[4] " << setValues[5] << "return " << returndecode << std::endl;

    // 2. 
}*/


vector<block> share_trivialFSS(uint64_t delta, uint64_t grid_x, uint64_t grid_y, uint64_t point_x, bool x, uint64_t point_y, bool y, block rand){
    //std::cout << "Let's write the trivial FSS in 2 Dimensions " << std::endl;
    BitVector FSS_keyzero, FSS_keyone;
    FSS_keyzero.assign(rand);
    FSS_keyone.assign(rand);

    //(grid_x, grid_y) = grid label
    grid_x = grid_x * (2 * delta);
    grid_y = grid_y * (2 * delta);
    uint64_t int_start = 64 - (2 * delta);
 

    //Let's process the X-coord which is the second half of the block, access block.mData[0]
    // in a block(y, x)
    if (x) { // fill right
    uint64_t start_point_x = int_start + (point_x - grid_x); 
    for (int i = start_point_x; i < 64; i++)
        if(FSS_keyone[64 + i])
            FSS_keyone[64 + i] = 0;
        else 
            FSS_keyone[64 + i] = 1;
    }
    else { // fill left 
     uint64_t end_point_x = int_start + (point_x - grid_x); 
        for (int i = int_start; i < end_point_x; i++)
        if(FSS_keyone[64 + i])
            FSS_keyone[64 + i] = 0;
        else 
            FSS_keyone[64 + i] = 1;
        
    }

    if (y) { // fill top
    uint64_t start_point_y = int_start + (point_y - grid_y); 
    for (int i = start_point_y; i < 64; i++)
        if(FSS_keyone[i])
            FSS_keyone[i] = 0;
        else 
            FSS_keyone[i] = 1;
    }
    else { // fill bottom
     uint64_t end_point_y = int_start + (point_y - grid_y); 
        for (int i = int_start; i < end_point_y; i++)
        if(FSS_keyone[i])
            FSS_keyone[i] = 0;
        else 
            FSS_keyone[i] = 1;
        
    }
    //std::cout << "FSS interval  " << (FSS_keyzero ^ FSS_keyone) << std::endl;

    // converting the BitVector -> Block
    u8 * u8_key0 = FSS_keyzero.data();
    u8 * u8_key1 = FSS_keyone.data();
    block k0(u8_key0[15],u8_key0[14],u8_key0[13],u8_key0[12],u8_key0[11],u8_key0[10],u8_key0[9],u8_key0[8],u8_key0[7],u8_key0[6],u8_key0[5],u8_key0[4],u8_key0[3],u8_key0[2],u8_key0[1],u8_key0[0]);
    block k1(u8_key1[15],u8_key1[14],u8_key1[13],u8_key1[12],u8_key1[11],u8_key1[10],u8_key1[9],u8_key1[8],u8_key1[7],u8_key1[6],u8_key1[5],u8_key1[4],u8_key1[3],u8_key1[2],u8_key1[1],u8_key1[0]);
    vector<block> fss_shares;
    fss_shares.push_back(k0);
    fss_shares.push_back(k1);
    return fss_shares;
}

/*
    This implementation assumes 2-dimensions
*/
void far_apart_FssShare(uint64_t delta, int nSquares, vector<block> &okvs0, vector<block> &okvs1){

    std::cout << "balls are pairwise 3*delta apart " << std::endl;
    //initialize variables
    uint64_t len_sqr = 2 * delta;
    vector<block> okvsKeys, okvsVal0, okvsVal1, shares; // vals0 - fsskeys0, vals1 - fsskeys1
    uint64_t bl_x, bl_y, br_x, br_y, tl_x, tl_y, tr_x, tr_y, grd_bl_x, grd_bl_y, grd_br_x, grd_br_y, grd_tl_x, grd_tl_y, grd_tr_x, grd_tr_y;
    block rand;
    PRNG mprng(toBlock(31));

    // let's process the squares: identify the 2^dim = 4 grid cells -> compute fsskey0, fsskey1 
    // Input: square1 (delta, delta), square2 (delta + 4*delta, delta) ..... square_n(delta + 4(n -1)delta, delta) ...
    bl_x = delta;
    bl_y = delta; 
    for (int i = 0; i < nSquares; i++){
        //BL 
        bl_x = delta + (i*2*len_sqr);
        bl_y = delta; // make bl_y = delta + (i*2*len_sqr)
        return_grid(bl_x, bl_y, grd_bl_x, grd_bl_y, len_sqr);
        rand = mprng.get();
        shares = share_trivialFSS(delta, grd_bl_x, grd_bl_y, bl_x, true, bl_y, true, rand);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
        okvsKeys.push_back(toBlock(grd_bl_y, grd_bl_x));
        //std::cout << "keys " << okvsKeys.back() << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        
        //TL
        tl_x = bl_x;
        tl_y = bl_y + len_sqr;
        return_grid(tl_x, tl_y, grd_tl_x, grd_tl_y, len_sqr);
        rand = mprng.get();
        shares = share_trivialFSS(delta, grd_tl_x, grd_tl_y, tl_x, true, tl_y, false, rand);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
        okvsKeys.push_back(toBlock(grd_tl_y, grd_tl_x));
        //std::cout << "keys " << okvsKeys.back() << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        
        //BR
        br_x = bl_x + len_sqr;
        br_y = bl_y;
        return_grid(br_x, br_y, grd_br_x, grd_br_y, len_sqr);
        rand = mprng.get();
        shares = share_trivialFSS(delta, grd_br_x, grd_br_y, br_x, false, br_y, true, rand);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
        okvsKeys.push_back(toBlock(grd_br_y, grd_br_x));
        //std::cout << "keys " << okvsKeys.back() << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        //std::cout << "grd id " << grd_br_x << "  " << grd_br_y << std::endl;

        //TR
        tr_x = bl_x + len_sqr;
        tr_y = bl_y + len_sqr;
        return_grid(tr_x, tr_y, grd_tr_x, grd_tr_y, len_sqr);
        rand = mprng.get();
        shares = share_trivialFSS(delta, grd_tr_x, grd_tr_y, tr_x, false, tr_y, false, rand);
        okvsVal0.push_back(shares[0]);
        okvsVal1.push_back(shares[1]);
        okvsKeys.push_back(toBlock(grd_tr_y, grd_tr_x));
        //std::cout << "keys " << okvsKeys.back() << " " << okvsVal0.back()  << "  " << okvsVal1.back() << std::endl;
        //std::cout << "grd id " << grd_tr_x << "  " << grd_tr_y << std::endl;

    }
    //std::cout << "sizes " << okvsKeys.size() << "  " << okvsVal0.size()  << "  " << okvsVal1.size() << std::endl;

    // OKVS (keys, values = fsskeys0) || (keys, values = fsskeys1)
    PaxosEncode(okvsKeys, okvsVal0, okvsVal1, okvs0, okvs1, 128);
    //std::cout << "key okvsVals0[0] " << okvsKeys[0] << "  " << okvsVal0[0] << std::endl;
}
   
void far_apart_FssEval(uint64_t x_coord, uint64_t y_coord, vector<block> okvs, uint64_t delta, uint64_t hashSize){
    std::cout << "we are now evaluating the FSS, far apart " << std::endl;  
    uint64_t len_sqr = 2* delta;
    uint64_t grd_x, grd_y, x, y; 
    int gamma = 60, v=20, fieldSizeBytes = 16, fieldSize = 128; // for okvs
    double c1 = 2.4; // for okvs

    //handling the key
    x = x_coord;
    y = y_coord;
    return_grid(x, y, grd_x, grd_y, len_sqr);
    block okvskey(toBlock(grd_y, grd_x));
    std::cout << "okvskey" << okvskey << std::endl;
    uint64_t okvs_key = XXH64(&okvskey, 64, 0);

    //handling okvs object 
    // converting block to bytes
    vector<byte> okvs_bytes;
    okvs_bytes.resize(okvs.size() * sizeof(block));
    memcpy(okvs_bytes.data(), okvs.data(), okvs.size() * sizeof(block));
    // from bytes to GF2E!! yahooo 
    GF2EVector okvs_gf2e(okvs.size());
    GF2X temp;
    for (int i=0; i < okvs_gf2e.size(); i++){
        GF2XFromBytes(temp, okvs_bytes.data() + i*fieldSizeBytes ,fieldSizeBytes);
        okvs_gf2e[i] = to_GF2E(temp); 
        
    }
    std::cout << okvs[0] << std::endl;
    // decoding from checkoutput()

    //initialize just to be able to dec() -- below 3 calls only set up the hash functions
    ObliviousDictionary * dict = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dict->init();
    auto indices = dict->dec(okvs_key);
    cout << "ENCODE " << okvskey << " " << okvs_key << " " << indices[0] << " " << indices[1] << std::endl;
    GF2E dhBitsVals;
    for (int j=0; j<indices.size(); j++){
        //std::cout << indices[j] << std::endl;
        dhBitsVals += okvs_gf2e[indices[j]];
    }
    //std::cout << "dhBitsVals  " << dhBitsVals << std::endl;
    // can generalize this later
    vector<byte> valBytes(fieldSizeBytes);
    BytesFromGF2X(valBytes.data(), rep(dhBitsVals), fieldSizeBytes);
    block returndecode(valBytes[15],valBytes[14],valBytes[13],valBytes[12],valBytes[11],valBytes[10],valBytes[9],valBytes[8],valBytes[7],valBytes[6],valBytes[5],valBytes[4],valBytes[3],valBytes[2],valBytes[1],valBytes[0]);
    std::cout <<  "return block from okvs " << returndecode << std::endl;

}
   
   
   
    /*
    
    block test = toBlock(14, 1000);
    vector<uint64_t> coord;
    uint64_t a = 12345;
    uint64_t b = 0;
    coord.push_back(a);
    coord.push_back(b);
    vector<uint64_t> grid = return_grid(coord, delta);
    

    vector<block> shares = d_trivialFSS(10, 0, 0, 15, true, 5, true, test);
    BitVector share0, share1;
    share0.assign(shares[0]);
    share1.assign(shares[1]);
    std:: cout << "printing trivial FSS for square  " << (share0 ^ share1) << std::endl;
*/





/* useful notes
    PRNG mprng(toBlock(31));
    vector<uint64_t> coord;
    uint64_t a = 12345;
    uint64_t b = 0;
    coord.push_back(a);
    coord.push_back(b);
    vector<uint64_t> grid = return_grid(coord, delta);
    
    // trivial FSS
    block test = toBlock(a, b);
    std::cout << "a " << a << " b " << b << std::endl;
    std::cout << "test " << test << std::endl;
    std::cout << "a " << test.mData[1] << "  " << test.mData[0] << std::endl;

    // testing bitvector of size 10 where first half is 1 and second half is 0
    BitVector test1;
    
    for (int i = 0; i < 5; i++)
        test1.pushBack(1);
    for (int i = 0; i < 5; i++)
        test1.pushBack(0);

        
    test1.assign(test);
    std::cout << "bit vector " << test1 << std::endl;
    u8 * u8_p = test1.data();
    block g(u8_p[15],u8_p[14],u8_p[13],u8_p[12],u8_p[11],u8_p[10],u8_p[9],u8_p[8],u8_p[7],u8_p[6],u8_p[5],u8_p[4],u8_p[3],u8_p[2],u8_p[1],u8_p[0]);
    std::cout << "block " << g << "   actual block " << test << std::endl;

    vector<block> shares = d_trivialFSS(10, 0, 0, 15, true, 5, true, test);
    BitVector share0, share1;
    share0.assign(shares[0]);
    share1.assign(shares[1]);
    std:: cout << "printing trivial FSS for square  " << (share0 ^ share1) << std::endl;


*/