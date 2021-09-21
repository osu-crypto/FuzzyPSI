#include "fuzzyPSI.h"

using namespace osuCrypto;
//using namespace cryptoTools;

/*
void test_paxos() {
	int v = 1; 
	int hashSize = 5;
    int fieldSize = 1;
    int gamma = 60;
    double c1 = 1.0;
    // Initialize
    OBD3Tables* dic = new OBD3tables(hashSize, fieldSize, c1, gamma, v);
    vector<uint64_t> keys{1,2,3,4,5};
    vector<byte> values{'a','b','c','d', 'e'};
    dic->encode();
    unsigned char ch = disc->decode(3);
    std::cout<<"3rd character: "<<ch;
}

*/

void fuzzyPSI(u64 keysize)
{
    //space to test out dependencies
    //fss_interval(1, 10);

    // Setup networking. Setting up channels for PSI, sender and recver according to OT not PSI 
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // key length = # of blocks needed to represent the FSS keys
    // determine this according to the FSS function we calling PSI on
    u64 fsskeySize =  keysize; 

    // The number of OTs - we need \kappa OTs
    int baseCount = 1; // for  BFSS - baseCount = hamming distance
    

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
        // delete below 
        for (int i = 0; i < recv_ciphertxt0.size(); i++)
            std::cout << "recv side cipher0 " << i << "   " << recv_ciphertxt0[i] << std::endl; 

        for (int i = 0; i < recv_ciphertxt1.size(); i++)
            std::cout << "recv side cipher1 " << i << "   " << recv_ciphertxt1[i] << std::endl; 

        for (int i = 0; i < fsskeysRecv.size(); i++)
            std::cout << "recv side fsskeys " << i << "   " << fsskeysRecv[i] << std::endl; 


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
    std::vector<std::array<osuCrypto::block, 2>> fsskeys(baseCount*fsskeySize);
    std::vector<block> keys_ciphertxt0, keys_ciphertxt1;
    s_prng.get(fsskeys.data(), fsskeys.size());

    AES aeskey0, aeskey1;
    block ciphertxt0, ciphertxt1;  
    //encrypt the fsskeys using baseSend(encryption keys)
    for (int i = 0; i < baseCount; i++){
        aeskey0.setKey(baseSend[i][0]);
        aeskey1.setKey(baseSend[i][1]);
        for (int j = 0; j < fsskeySize; j++){
            std::cout << "sender: fss key (i, j) " << fsskeys[(i*fsskeySize) + j][0] << " " << fsskeys[(i*fsskeySize) + j][1] << std::endl; 
            aeskey0.ecbEncBlock(fsskeys[(i*fsskeySize) + j][0], ciphertxt0);
            aeskey1.ecbEncBlock(fsskeys[(i*fsskeySize) + j][1], ciphertxt1);
            std::cout << "sender: fss cipher (i, j) " << ciphertxt0 << " " << ciphertxt1 << std::endl; 
            keys_ciphertxt0.push_back(ciphertxt0);
            keys_ciphertxt1.push_back(ciphertxt1);
            
        }
    }
    senderChl.asyncSend(keys_ciphertxt0);
    senderChl.asyncSend(keys_ciphertxt1);
    
    
    recverThread.join();


}


/*    std::cout << "testing chosen base OT " << std::endl; 
    for (u64 i = 0; i < (baseCount*fsskeySize); ++i)
        {
            
            if (neq(baseRecv[i], baseSend[i][choices[i]]))
            {
                std::cout << "failed " << i << " exp = m[" << int(choices[i]) << "], act = " << baseRecv[i] << " true = " << baseSend[i][0] << ", " << baseSend[i][1] << std::endl;
            }
            
        }

        //making fsskeys blocks of choices same value
        for (u64 i = 0; i < baseCount; i++){
            for (u64 j = 1; j < fsskeySize; j++)
                choices[(i*fsskeySize) + j] = choices[(i*fsskeySize) + 0];
        }
        std::cout << "I am in the receiver side sending my hash seed " << r_HashSeed << std::endl; 
        std::cout << "I am in the receiver side receiving my hash seed " << theirHashseed << std::endl; 
        std::cout << "receiver checking hashing " << ciphertxt << std::endl;



        senderChl.recv(recvd_outputs);
        for (int i = 0; i < recvd_outputs.size(); i++)
          std::cout << " hash values " << recvd_outputs[i] << std::endl; 

 for (int i = 0; i < baseCount; i++){
            for (int j = 0; j < fsskeySize; j++){
                hasher.ecbEncBlock(baseRecv[(i*fsskeySize) + j], ciphertxt);
                hashed_item.push_back(ciphertxt);
                //hasher.ecbEncBlocks(baseRecv.data(), baseRecv.size(), ciphertxts.data());
            }
            for (int k = 0; k < hashed_item.size(); k++) // equal to 
                hash = hash ^ hashed_item[k];
            hashed_inputs.push_back(hash);
            hashed_item.clear();
        }
        recverChl.asyncSend(hashed_inputs);      
        for (int i = 0; i < hashed_inputs.size(); i++)
          std::cout << " hash values " << hashed_inputs[i] << std::endl; 

        */