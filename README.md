

A fast and portable C++17 library for Oblivious Transfer extension (OTe). The 
primary design goal of this library to obtain *high performance* while being 
*easy to use*.  This library currently implements:
 
 
## Introduction
 
This library provides several different classes of OT protocols. First is the 
base OT protocol of [NP01, CO15, MR19, MRR21]. These protocol bootstraps all the other
OT extension protocols.  Within the OT extension protocols, we have 1-out-of-2,
1-out-of-N and K-out-of-N, both in the semi-honest and malicious settings. See The `frontend` or `libOTe_Tests` folder for examples.

All implementations are highly optimized using fast SSE instructions and vectorization
to obtain optimal performance both in the single and multi-threaded setting. See 
the **Performance** section for a comparison between protocols and to other libraries. 
 
Networking can be performed using both the sockets provided by the library and
external socket classes. See the [networking tutorial](https://github.com/ladnir/cryptoTools/blob/57220fc45252d089a7fd90816144e447a2ce02b8/frontend_cryptoTools/Tutorials/Network.cpp#L264) for an example.


## Build
 
The library is *cross platform* and has been tested on Windows, Mac and Linux. 
There is one mandatory dependency on [Boost 1.75](http://www.boost.org/) (networking),
and three **optional dependencies** on [libsodium](https://doc.libsodium.org/),
[Relic](https://github.com/relic-toolkit/relic), or
[SimplestOT](https://github.com/osu-crypto/libOTe/tree/master/SimplestOT) (Unix only)
for Base OTs.
CMake 3.15+ is required and the build script assumes python 3.
 
The library can be built as
```
git clone --recursive https://github.com/osu-crypto/libOTe.git
cd libOTe
python build.py --setup --boost --relic
python build.py -- -D ENABLE_RELIC=ON -D ENABLE_ALL_OT=ON
```
The main executable with examples is `frontend` and is located in the build directory, eg `out/build/linux/frontend/frontend.exe, out/build/x64-Release/frontend/Release/frontend.exe` depending on the OS. 

### Dependancies

Dependancies can be managed via the `build.py` script or or installed via an external tool. If an external tool is used install to system location or set  `-D CMAKE_PREFIX_PATH=path/to/install`.


## Install

libOTe can be installed and linked the same way as other cmake projects. By default the dependancies are not installed. To install all dependancies, run the following
```
python build.py --setup --boost --relic --sodium --install
```
You can also selectively install the dependancies. The install location can be specifying as `--install=path/to/installation`. Otherwise the system default is used.

The main library is similarly installed as
```
python build.py --install 
```

By default, sudo is not used. If installation requires sudo access, then add `--sudo` to the `build.py` script arguments. See `python build.py --help` for full details.



## Help
 
Contact Gayathri Garimella garimelg@oregonstate.edu for any assistance on building 
or running the library.

 
 ## License
 
This project has been placed in the public domain and/or MIT license. As such, you are unrestricted in how 
you use it, commercial or otherwise. However, no warranty of fitness is provided. If you 
found this project helpful, feel free to spread the word and cite us.

