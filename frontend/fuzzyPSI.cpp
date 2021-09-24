#include "fuzzyPSI.h"
//#include "paxos.h"

using namespace osuCrypto;
//using namespace cryptoTools;

/* experiment, issue is reusing te pointers, fix later
array<vector<block>, 2> s_psi_interval(Fss* f_sender, int n, uint64_t a, uint64_t b){

    array<vector<block>, 2> return_blks;
    vector<block> int_keys0, int_keys1;
    block * key_in_blks; 
    //pair of key
    ServerKeyLt l_lt_k0, l_lt_k1, r_lt_k0, r_lt_k1;
    for (int i = 0; i < n; i++){
        generateTreeLt(f_sender, &r_lt_k0, &r_lt_k1, b, 1);  
        key_in_blks = reinterpret_cast<block *>(&r_lt_k1); 
        for (int i = 0; i < 11; i++){
            int_keys1.push_back(key_in_blks[i]);
            std::cout << "one keys  " << key_in_blks[i] << std::endl;
        }
        key_in_blks = reinterpret_cast<block *>(&r_lt_k0); 
        for (int i = 0; i < 11; i++){
            int_keys0.push_back(key_in_blks[i]);
            std::cout << "zero keys  " << key_in_blks[i] << std::endl;
        }
    }
    return_blks[0] = int_keys0;
    return_blks[1] = int_keys1;
    return return_blks; 

    
}*/
void PaxosEncode(const std::vector<block> setKeys, const std::vector<block> setValues, std::vector<block>& okvs, uint64_t fieldSize)
{
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

    cout << "Copied vals" << endl;
    cout << "printting the keys, values size " << keys.size() << "    " << values.size() << std::endl;
    dic->setKeysAndVals(keys, values);
    dic->encode();

    //vector<byte> x = dic->getVariables();
    //memcpy(okvs, x.data(), x.size());
    //cout << "printting the okvs size " << okvs.size() << std::endl;
    dic->checkOutput();
    cout << "here" << endl;

}

void PaxosDecode(const std::vector<block> paxosMat, const std::vector<block> setKeys, std::vector<block>& setValues)
{
    cout << "Decode" << endl;
}
void fss_interval3(uint64_t a, uint64_t b){
    uint64_t lboundary = a;
    uint64_t rboundary = b;
    uint64_t l_lt_ans0, l_lt_ans1, l_lt_fin, r_lt_ans0, r_lt_ans1, r_lt_fin;
    Fss l_fClient, l_fServer, r_fClient, r_fServer;

    ServerKeyLt l_lt_k0;
    ServerKeyLt l_lt_k1;
    ServerKeyLt r_lt_k0;
    ServerKeyLt r_lt_k1;

    
    // Client does the key generation for the FSS

    // left interval (a, -1), x < a then - 1, x >= a then 0 
    common_test();
    initializeClient(&l_fClient, 64, 2);
    generateTreeLt(&l_fClient, &r_lt_k0, &r_lt_k1, b, 1);


    block* struct_in_blk = reinterpret_cast<block *>(&r_lt_k1); 
    vector<block> struct_in_blocks;
    for (int i = 0; i < 11; i++){
        std::cout << "blocks in struct " << struct_in_blk[i] << std::endl; 
        struct_in_blocks.push_back(struct_in_blk[i]);
    }
    generateTreeLt(&l_fClient, &r_lt_k0, &r_lt_k1, b, 1);
    block* struct_in_blk2 = reinterpret_cast<block *>(&r_lt_k1); 
    vector<block> struct_in_blocks2;
    for (int i = 0; i < 11; i++){
        std::cout << "blocks in struct2 " << struct_in_blk2[i] << std::endl; 
        struct_in_blocks2.push_back(struct_in_blk2[i]);
    }

    ServerKeyLt * original_struct = reinterpret_cast<ServerKeyLt*>(struct_in_blocks.data());
    l_lt_ans0 = evaluateLt(&l_fClient, &r_lt_k1, (a+1));
    r_lt_ans0 = evaluateLt(&l_fClient, &r_lt_k0, (a+1));
    l_lt_ans1 = evaluateLt(&l_fClient, original_struct, (a+1));

    std::cout << l_lt_ans0 << "   " << r_lt_ans0 << "   " << l_lt_ans1 << std::endl; 
}

void fuzzyPSI(u64 keysize, u64 y_input_size, u64 x_volume)
{
    //space to test out dependencies
    common_test();
    PRNG pprng(toBlock(123));
    vector<block> okvskeys(65536), okvsvals(65536), okvs;
    pprng.get(okvskeys.data(), okvskeys.size());
    pprng.get(okvsvals.data(), okvsvals.size());
    uint64_t fieldsize = 65;
    PaxosEncode(okvskeys, okvsvals, okvs, fieldsize);
    //test_paxos();
    //myfunc(toBlock(14), toBlock(15));

    // Setup networking. Setting up channels for PSI, sender and recver according to OT not PSI 
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // key length = # of blocks needed to represent the FSS keys
    // determine this according to the FSS function we calling PSI on
    // for a single interval the **keysize is 22**
    u64 fsskeySize =  2; 
    u64 x_input_volume = x_volume;

    // The number of OTs - we need \kappa OTs
    int baseCount = 2; // for  BFSS - baseCount = hamming distance
    
    // **PSI receiver** has 128 choices and OT received blocks 
    // **PSI receiver** has an unstructured set r_inputs as inputs
    // we just call the Hash * r_inputs 
    std::vector<osuCrypto::block> baseRecv(baseCount);
    BitVector choices(baseCount);
    u64 r_input_size = y_input_size;

    // **PSI sender** operates out of the receiver thread
    auto recverThread = std::thread([&]() {

        //initializing
        block ciphertxt, hash, theirHashseed, r_psiHash;
        

        // sampling the l - bit vector or choices as the OT receiver 
        PRNG r_prng(toBlock(14));
        choices.randomize(r_prng);
        MasnyRindal recver;
        
        // Receive the messages, establish **common hash key**
        recver.receiveChosen(choices, baseRecv, r_prng, recverChl);
        

        // basically receiver ciphertexts of form # base OT * [enc(k1, m1), enc(k2, m2)] here
        std::vector<block> recv_ciphertxt0, recv_ciphertxt1, fsskeysRecv;
        block recv_ciphertxt;
        recverChl.recv(recv_ciphertxt0);
        recverChl.recv(recv_ciphertxt1);

        //now we use the baseRecv as decryption key and choices 
        details::AESTypes::Portable;
        AESDec aeskey;
            for (int i = 0; i < baseCount; i++){
                aeskey.setKey(baseRecv[i]);
                for (int j = 0; j < fsskeySize; j++){
                        if (choices[i] == 0)
                            recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt0[(i*fsskeySize) + j]);  
                        else 
                            recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt1[(i*fsskeySize) + j]);
                        fsskeysRecv.push_back(recv_ciphertxt);

                }
            }

        
        // initializing variables, getting common hash seed
        std::vector<block> hashed_inputs, hashed_item(fsskeySize*baseCount+1);
        block r_HashSeed = r_prng.get<block>();
        details::AESTypes::Portable;
        recverChl.asyncSend(r_HashSeed);
        recverChl.recv(theirHashseed);
        r_psiHash = r_HashSeed ^ theirHashseed;
        AES hasher(r_psiHash);


        // hash each of the received keys with the appropriate inputs
        // ^^insert *FSS eval for each of the r_inputs happens here*
        // we Hash size (fsskeySize * baseCount )  ---> size (baseCount) for a single inputs!

        // suppose we have r_inputs_size, append input and hash (for now)
        
        for (u64 r = 0; r < r_input_size; r++){
            fsskeysRecv.push_back(toBlock(r));
            hasher.ecbEncBlocks(fsskeysRecv.data(), fsskeySize*baseCount+1, hashed_item.data());
            fsskeysRecv.pop_back();
            hash = hashed_item[0];
            for (int k = 1; k < hashed_item.size(); k++){ // here we xor all blocks associated with "one fsskey", there are baseCount such keys
                hash = hash ^ hashed_item[k];
            }
            hashed_inputs.push_back(hash);
        }
        //sender sending hash of FSS eval of his inputs
        recverChl.asyncSend(hashed_inputs); 

        });

    // *PSI receiver*
    block theirHash, s_HashSeed, s_psiHash, s_ciphertxt;
    std::vector<block> recvd_outputs;
    
    //initializing *OT sender*
    PRNG s_prng(toBlock(12));
    MasnyRindal sender;
    std::vector<std::array<osuCrypto::block, 2>> baseSend(baseCount);
    s_prng.get(baseSend.data(), baseSend.size());

    // Send the messages, that is, short keys (k0, k1) that will encrypt the actual FSS keys (FSSkey0, FSSkey1)
    sender.sendChosen(baseSend, s_prng, senderChl);
    
    //sampling uniformly random FSS keys ---> later this should come from some function call
    //Fss f_sender;
    vector<block> fsskeys0(baseCount*fsskeySize), fsskeys1(baseCount*fsskeySize);
    //= s_psi_interval(&f_sender, 1, 5, 10);
    std::vector<block> keys_ciphertxt0, keys_ciphertxt1;
    s_prng.get(fsskeys0.data(), fsskeys0.size());
    s_prng.get(fsskeys1.data(), fsskeys1.size());

   
    AES aeskey0, aeskey1;
    block ciphertxt0, ciphertxt1;  
    //encrypt the fsskeys using baseSend(encryption keys)
    for (int i = 0; i < baseCount; i++){
        aeskey0.setKey(baseSend[i][0]);
        aeskey1.setKey(baseSend[i][1]);
        for (int j = 0; j < fsskeySize; j++){
            aeskey0.ecbEncBlock(fsskeys0[(i*fsskeySize) + j], ciphertxt0);
            aeskey1.ecbEncBlock(fsskeys1[(i*fsskeySize) + j], ciphertxt1); 
            keys_ciphertxt0.push_back(ciphertxt0);
            keys_ciphertxt1.push_back(ciphertxt1);
        }
    }
    senderChl.asyncSend(keys_ciphertxt0);
    senderChl.asyncSend(keys_ciphertxt1);
    

    // establishing common hash seed for the last step 
    s_HashSeed = s_prng.get<block>();
    senderChl.recv(theirHash);
    senderChl.asyncSend(s_HashSeed);
    s_psiHash = s_HashSeed ^ theirHash; 
    AES s_hasher(s_psiHash);
    
    // now we recv the hash of the each of the PSI sender's inputs 
    vector<block> recvHashinputs(r_input_size), recver_hashes, s_hashed_item(fsskeySize*baseCount+1);
    block s_hash; 
    senderChl.recv(recvHashinputs);
    recverThread.join();
    // **PSI receiver learns the actual output here**
    // # idea 
    // use fsskeys0 and compute the hashes, |hashes| = x_volume (#structured input) and sort the vector hashinput

    for (u64 r = 0; r < x_input_volume; r++){
            fsskeys0.push_back(toBlock(r));
            s_hasher.ecbEncBlocks(fsskeys0.data(), fsskeySize*baseCount+1, s_hashed_item.data());
            fsskeys0.pop_back();
            s_hash = s_hashed_item[0];
            for (int k = 1; k < s_hashed_item.size(); k++) // here we xor all blocks associated with "one input"
                s_hash = s_hash ^ s_hashed_item[k];
            recver_hashes.push_back(s_hash);
        }

    
    std::sort(recver_hashes.begin(), recver_hashes.end());


    // search rechHashinputs against recver_hashes
    vector<block> psi_outputs; 
    int output_size = 0;
    for (u64 i = 0; i < recvHashinputs.size(); i++){
        if(binary_search(recver_hashes.begin(), recver_hashes.end(), recvHashinputs[i])){
            psi_outputs.push_back(recvHashinputs[i]);
            output_size = output_size + 1;
        }
        else 
            cout << "problematic point is " << i << std::endl;
    }

    if(psi_outputs.size() == recvHashinputs.size())
        cout << " PSI WORKS " << std::endl;
    else{
        cout << " PSI FAILS " << std::endl;
        cout << " output size " << output_size << std::endl;
    }

}


void fuzzyPSI_int(u64 keysize)
{
    //space to test out dependencies
    //fss_interval3(1, 10);

    // Setup networking. Setting up channels for PSI, sender and recver according to OT not PSI 
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // key length = # of blocks needed to represent the FSS keys
    // determine this according to the FSS function we calling PSI on
    // for a single interval the **keysize is 22**
    u64 fsskeySize =  11; 

    // The number of OTs - we need \kappa OTs
    int baseCount = 2; // for  BFSS - baseCount = hamming distance
    
    // **PSI receiver** has 128 choices and OT received blocks 
    std::vector<osuCrypto::block> baseRecv(baseCount);
    BitVector choices(baseCount);

    // **PSI sender** operates out of the receiver thread
    auto recverThread = std::thread([&]() {

        //initializing
        block ciphertxt, hash, theirHashseed, r_psiHash;
        std::vector<block> hashed_inputs, hashed_item(fsskeySize);

        // sampling the l - bit vector or choices as the OT receiver 
        PRNG r_prng(toBlock(14));
        choices.randomize(r_prng);
        MasnyRindal recver;
        
        // Receive the messages, establish **common hash key**
        recver.receiveChosen(choices, baseRecv, r_prng, recverChl);
        block r_HashSeed = r_prng.get<block>();
        details::AESTypes::Portable;
        recverChl.asyncSend(r_HashSeed);
        recverChl.recv(theirHashseed);
        r_psiHash = r_HashSeed ^ theirHashseed;
        AES hasher(r_psiHash);

        // basically receiver ciphertexts of form # base OT * [enc(k1, m1), enc(k2, m2)] here
        std::vector<block> recv_ciphertxt0, recv_ciphertxt1, fsskeysRecv;
        block recv_ciphertxt;
        recverChl.recv(recv_ciphertxt0);
        recverChl.recv(recv_ciphertxt1);

        //now we use the baseRecv as decryption key and choices 
        details::AESTypes::Portable;
        AESDec aeskey;
            for (int i = 0; i < baseCount; i++){
                aeskey.setKey(baseRecv[i]);
                for (int j = 0; j < fsskeySize; j++){
                        if (choices[i] == 0)
                            recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt0[(i*fsskeySize) + j]);  
                        else 
                            recv_ciphertxt = aeskey.ecbDecBlock(recv_ciphertxt1[(i*fsskeySize) + j]);
                        fsskeysRecv.push_back(recv_ciphertxt);
                }
            }

        // hash each of the received keys with the appropriate inputs
       
        });

    // *PSI receiver*
    block theirHash, s_HashSeed, s_psiHash, s_ciphertxt;
    std::vector<block> recvd_outputs;
    
    //initializing *OT sender*
    PRNG s_prng(toBlock(12));
    MasnyRindal sender;
    std::vector<std::array<osuCrypto::block, 2>> baseSend(baseCount);
    s_prng.get(baseSend.data(), baseSend.size());

    // Send the messages, that is, short keys (k0, k1) that will encrypt the actual FSS keys (FSSkey0, FSSkey1)
    sender.sendChosen(baseSend, s_prng, senderChl);
    s_HashSeed = s_prng.get<block>();
    senderChl.recv(theirHash);
    senderChl.asyncSend(s_HashSeed);
    s_psiHash = s_HashSeed ^ theirHash; 
    AES s_hasher(s_psiHash);
    
    //sampling uniformly random FSS keys ---> later this should come from some function call
    //Fss f_sender;
    vector<block> fsskeys0(baseCount*fsskeySize), fsskeys1(baseCount*fsskeySize);
    //= s_psi_interval(&f_sender, 1, 5, 10);
    std::vector<block> keys_ciphertxt0, keys_ciphertxt1;
    s_prng.get(fsskeys0.data(), fsskeys0.size());
    s_prng.get(fsskeys1.data(), fsskeys1.size());
   
    AES aeskey0, aeskey1;
    block ciphertxt0, ciphertxt1;  
    //encrypt the fsskeys using baseSend(encryption keys)
    for (int i = 0; i < baseCount; i++){
        aeskey0.setKey(baseSend[i][0]);
        aeskey1.setKey(baseSend[i][1]);
        for (int j = 0; j < fsskeySize; j++){
            aeskey0.ecbEncBlock(fsskeys0[(i*fsskeySize) + j], ciphertxt0);
            aeskey1.ecbEncBlock(fsskeys1[(i*fsskeySize) + j], ciphertxt1); 
            keys_ciphertxt0.push_back(ciphertxt0);
            keys_ciphertxt1.push_back(ciphertxt1);
        }
    }
    senderChl.asyncSend(keys_ciphertxt0);
    senderChl.asyncSend(keys_ciphertxt1);
    
    
    recverThread.join();


}


