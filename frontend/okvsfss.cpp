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

void PaxosEncode(const std::vector<block> setKeys, const std::vector<block> setValues, std::vector<block>& okvs, uint64_t fieldSize)
{
    GF2E dhBitsVal;
    int hashSize=setKeys.size(), gamma = 60, v=20;
    double c1 = 2.4;
    vector<uint64_t> keys;
//    vector<unsigned char> values;
    keys.resize(hashSize);
    //okvs.resize(c1*hashSize);
    int fieldSizeBytes = fieldSize % 8 == 0 ? fieldSize/8 : fieldSize/8 + 1;
    int zeroBits = 8 - fieldSize % 8;
//    values.resize(hashSize*fieldSizeBytes);
    ObliviousDictionary * dic = new OBD3Tables(hashSize, c1, fieldSize, gamma, v);
    dic->init();
    for (int i=0; i < setKeys.size(); i++){
        keys[i] = XXH64(&setKeys[i], 64, 0);
    }
    cout << "Done with keys" << endl;
    vector<byte> values;
    values.resize(setValues.size() * sizeof(block));
    memcpy(values.data(), setValues.data(), setValues.size() * sizeof(block));
    cout << "setValue size " << setValues.size() << std::endl;
   // cout << "block size " << sizeof(block) << "  fieldSize " << fieldSizeBytes << std::endl;
    cout << "printting the keys, values size " << keys.size() << "    " << values.size() << std::endl;
    dic->setKeysAndVals(keys, values);
    dic->encode();

    vector<byte> x = dic->getVariables();
    //std::cout << "field size bytes " << fieldSizeBytes << std::endl;

    // vector of bytes to block 
    vector<block> Okvs;
    Okvs.resize(x.size() / fieldSizeBytes);
    memcpy(Okvs.data(), x.data(), x.size());
    std::cout << "Okvs.size() " << Okvs.size() << std::endl;
    for (int i = 0; i < (Okvs.size() - 50); i++)
        std::cout << "okvs  " << i << "  " << Okvs[i] << std::endl;
    dic->checkOutput();
    cout << "test decode below" << endl;

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

    //if ((variables[indices[0]] + variables[indices[1]] +  variables[indices[2]] + dhBitsVal) == val) {
    //    }
} 

void OKVSDecode(vector<block> okvs, block key){
    cout << "Let's write decode " << std::endl;
    // 1. convert block key -> to unint64_t keys using xxH64
    // nothing changes in setKeysandValues()

    // 2. 
}

vector<uint64_t> return_grid(vector<uint64_t> coords, uint64_t delta){
    vector<uint64_t> grid_id;
    uint64_t temp;
    for (int i = 0; i < coords.size(); i++){
        temp = coords[i] / delta;
        grid_id.push_back(temp);
    }
    
    return grid_id;
}

vector<block> d_trivialFSS(uint64_t delta, uint64_t grid_x, uint64_t grid_y, uint64_t point_x, bool x, uint64_t point_y, bool y, block rand){
    std::cout << "Let's write the trivial FSS in 2 Dimensions " << std::endl;
    BitVector FSS_keyzero, FSS_keyone;
    FSS_keyzero.assign(rand);
    FSS_keyone.assign(rand);

    //(grid_x, grid_y) = grid label
    grid_x = grid_x * (2 * delta);
    grid_y = grid_y * (2 * delta);
    uint64_t int_start = 64 - (2 * delta);
 

    //Let's process the X-coord which is the second half of the block, access block.mData[0]
   
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

    if (y) { // fill right
    uint64_t start_point_y = int_start + (point_y - grid_y); 
    for (int i = start_point_y; i < 64; i++)
        if(FSS_keyone[i])
            FSS_keyone[i] = 0;
        else 
            FSS_keyone[i] = 1;
    }
    else { // fill left 
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


void ShareFss_farapart(uint64_t delta, int nSquares){
    std::cout << "balls are pairwise 3*delta apart " << std::endl;
    PRNG mprng(toBlock(31));
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

}




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