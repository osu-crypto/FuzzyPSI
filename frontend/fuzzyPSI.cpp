#include "fuzzyPSI.h"

using namespace osuCrypto;
//using namespace cryptoTools;

void fuzzyPSI(u64 keysize)
{
    // Setup networking. See cryptoTools\frontend_cryptoTools\Tutorials\Network.cpp
    IOService ios;
    Channel senderChl = Session(ios, "localhost:1212", SessionMode::Server).addChannel();
    Channel recverChl = Session(ios, "localhost:1212", SessionMode::Client).addChannel();

    // key length = # of blocks needed to represent the FSS keys
    u64 fsskeySize =  keysize; 

    // The number of OTs - we need \kappa OTs
    int baseCount = 3; // for  BFSS - baseCount = hamming distance 
    std::vector<osuCrypto::block> baseRecv(baseCount*fsskeySize);
    BitVector choices(baseCount*fsskeySize);

    // The code to be run by the OT receiver.
    auto recverThread = std::thread([&]() {
        PRNG r_prng(toBlock(14));
        MasnyRindal recver;

        // Choose which messages should be received.
        choices.randomize(r_prng);
        //making fsskeys blocks of choices same value
        for (u64 i = 0; i < baseCount; i++){
            for (u64 j = 1; j < fsskeySize; j++)
                choices[(i*fsskeySize) + j] = choices[(i*fsskeySize) + 0];
        }

        // Receive the messages
        recver.receiveChosen(choices, baseRecv, r_prng, recverChl);

        //AES.cpp has hash block, twoblocks, many blocks, woohoo 
        std::cout << "plain block " << baseRecv[0] << std::endl; 
        block ciphertxt;
        std::vector<block> ciphertxts;
        block rand = r_prng.get<block>();
        details::AESTypes::NI;
        AES hasher(rand);
        //hasher.AES(rand);
        //for (int i = 0; i < (baseCount*fsskeySize); i++){
          //  block cipertxt; 
        hasher.ecbEncBlock(baseRecv[0], ciphertxt);
        //hasher.ecbEncBlocks(baseRecv.data(), baseRecv.size(), ciphertxts.data());

       // }
                
        //for (int k = 0; k < (baseCount*fsskeySize); k++)
          std::cout << "plain text " << baseRecv[0] << "ciphertext " << ciphertxt << std::endl; 
        });

    PRNG s_prng(toBlock(12));
    MasnyRindal sender;

    // Choose which messages should be sent.
    std::vector<std::array<osuCrypto::block, 2>> baseSend(baseCount*fsskeySize);
    s_prng.get(baseSend.data(), baseSend.size());

    // Send the messages.
    sender.sendChosen(baseSend, s_prng, senderChl);
    recverThread.join();

    std::cout << "testing chosen base OT " << std::endl; 
    for (u64 i = 0; i < (baseCount*fsskeySize); ++i)
        {
            
            if (neq(baseRecv[i], baseSend[i][choices[i]]))
            {
                std::cout << "failed " << i << " exp = m[" << int(choices[i]) << "], act = " << baseRecv[i] << " true = " << baseSend[i][0] << ", " << baseSend[i][1] << std::endl;
            }
            
        }
}


/*void chosen_Base_OT(std::string ip, std::string tag, CLP&, BaseOT ot)
    {
        IOService ios;
        PRNG prng(sysRandomSeed());

        int totalOTS = 128; 

       // if (numThreads > 1)
        //    std::cout << "multi threading for the base OT example is not implemented.\n" << std::flush;

        Timer t;
        Timer::timeUnit s;
        if (role == Role::Receiver)
        {
            auto chl0 = Session(ios, ip, SessionMode::Server).addChannel();
            chl0.waitForConnection();
            BaseOT recv = ot;

            std::vector<block> msg(totalOTs);
            BitVector choice(totalOTs);
            choice.randomize(prng);


            s = t.setTimePoint("base OT start");

            recv.receive(choice, msg, prng, chl0);
        }
        else
        {

            auto chl1 = Session(ios, ip, SessionMode::Client).addChannel();
            chl1.waitForConnection();

            BaseOT send = ot;

            std::vector<std::array<block, 2>> msg(totalOTs);

            s = t.setTimePoint("base OT start");

            send.send(msg, prng, chl1);
        }

        auto e = t.setTimePoint("base OT end");
        auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(e - s).count();

        std::cout << tag << (role == Role::Receiver ? " (receiver)" : " (sender)")
            << " n=" << totalOTs << " " << milli << " ms" << std::endl;
    }


}
*/
