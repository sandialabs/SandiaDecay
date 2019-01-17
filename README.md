# SandiaDecay
Simple, fast, and versatile C++ library to perform nuclear decay calculations, retrieve nuclide information, and compute decay products.

[ENSDF](https://www.nndc.bnl.gov/ensdf/) provides a wealth of raw nuclear information, the [Bateman equations](https://en.wikipedia.org/wiki/Bateman_equation) are well known, but tasks such as determining what gammas will be produced, or the activities of daughter products from a 20 year old sample of 226-Ra, can still be formidable.   **SandiaDecay** makes performing calculations like these easy.

An example of using this library to perform these calculations is:
```c++
#include "SandiaDecay.h"
...

//Initialize the database from the XML file that has all the nuclear info
SandiaDecay::SandiaDecayDataBase database( "sandia.decay.xml" );

//Get pointer to Ra226. "Ra226", "Ra-226", "226Ra", etc are all equiv.
const SandiaDecay::Nuclide * const ra226 = database.nuclide( "Ra226" );

//Use a SandiaDecay::NuclideMixture object to perform the decay
SandiaDecay::NuclideMixture mixture;

//Add 1 mCi of Ra226 to the mixture
mixture.addNuclideByActivity( ra226, 0.001*SandiaDecay::curie );

//Get activities of Ra226 and all of its descendants after 20 years of aging
vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( 20*SandiaDecay::year );
for( size_t i = 0; i < activities.size(); ++i )
  cout << activities[i].nuclide->symbol << " has activity "
       << activities[i].activity/SandiaDecay::curie << " curies\n";

//Get all of the x-rays and gammas the sample would produce
vector<SandiaDecay::EnergyRatePair> photons = mixture.photons( 20*SandiaDecay::year,
                                                      SandiaDecay::NuclideMixture::OrderByEnergy );
for( size_t i = 0; i < photons.size(); ++i )
  cout << photons[i].numPerSecond << ", "
       << photons[i].energy << " keV photons/second\n";

```

See the [sandia_decay_example.cpp](https://github.com/sandialabs/SandiaDecay/tree/master/examples/sandia_decay_example.cpp) file for more examples of things you can do with this library.

**SandiaDecay** was developed as part of [InterSpec](https://github.com/sandialabs/InterSpec), but has since proven useful in a number of other applications.

## Getting Started
If all you care about is performing decay calculations, grab a copy of [InterSpec](https://github.com/sandialabs/InterSpec/releases) and open up its *Nuclide Decay Info* tool from the *Tools* menu.  This will let you perform many types of decay calculations, visualize decay chains, and export activities or decay products to CSVs.

If you would like to do complex or specialized calculations, or you want to add nuclide decay capabilities to your application, you will need to compile and link to this library.

The only dependancies required to compile are [rapidxml](http://rapidxml.sourceforge.net/) (which is included in the repository), and a C++&ge;03 compiler.  There is a [CMake](https://cmake.org/) based build file included, but because there are no configurations and the library only consists of a single header and source file, it should be simple to integrate in with your program no matter the build system you use.

To include this library with a CMake based project, just add the following to your CmakeLists.txt file:
```CMake
add_subdirectory( path/to/SandiaDecay )
...
target_link_libraries( MyExe PRIVATE SandiaDecay )
```

Or to compile and run the example use code:
```bash
cd /path/to/SandiaDecay
mkdir build; cd build
cmake ..
make
./examples/sandia_decay_example ../sandia.decay.xml
```

Or if you're on a UNIX system you could compile the example-use code using a simple command like:
```bash
g++ -I. -I3rdparty SandiaDecay.cpp ./examples/sandia_decay_example.cpp -o example.exe
```

## Efficiency
On a 2017 macbook pro it takes about 100 ms to to initialize the a `SandiaDecayDataBase` object using
[sandia.decay.nocoinc.min.xml](https://github.com/sandialabs/SandiaDecay/blob/master/sandia.decay.nocoinc.min.xml), or about 225 ms using [sandia.decay.xml](https://github.com/sandialabs/SandiaDecay/tree/master/sandia.decay.xml), with the differences between the XML files being how human readable they are and the former not including gamma coincidence information (to save file size).
On average it takes about 200 &mu;s to calculate a nuclides decay equations, which only needs to be done once for every nuclide you are interested in, and then about 20 &mu;s to evaluate the decay equations and determine the decay product or activities for given time.
A `SandiaDecayDataBase` object allocates about 10 MB of memory during initialization.  After initialization all operations are `const` operations, meaning a single `SandiaDecayDataBase` can be used in multiple threads without issues; no global variables are used, and except during construction of a `SandiaDecayDataBase` object, no files are accessed.


## Source Data
Source data in [sandia.decay.xml](https://github.com/sandialabs/SandiaDecay/tree/master/sandia.decay.xml) is originally based on the [ENSDF](https://www.nndc.bnl.gov/ensdf/) datasheets, with further corrections manually applied to compensate for parsing issues, or discrepencies in observed data by users of this library.

Parsing the ENSDF datasheets is a challenge, so there are almost certainly mistakes or bugs in [sandia.decay.xml](https://github.com/sandialabs/SandiaDecay/tree/master/sandia.decay.xml), and so the values returned by this library should not be considered definitive.   If you do find any issues, please email wcjohns@sandia.gov or interspec@sandia.gov.

## Contact
Will Johnson (wcjohns@sandia.gov or interspec@sandia.gov).

## Credit
Code is written by Will Johnson (SNL), with [sandia.decay.xml](https://github.com/sandialabs/SandiaDecay/tree/master/sandia.decay.xml) 
prepared by John Mattingly (now NCSU), Steve Horne (SNL), Greg Thoreson (SNL), Will Johnson, and others.

## License
The source code for this project is licensed under the LGPL v2.1 License - see the [LICENSE.md](LICENSE.md) file for details

## Copyright
Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC (NTESS).
Under the terms of Contract DE-NA0003525 with NTESS, the U.S. Government retains certain rights in this software.

## Disclaimer
```
DISCLAIMER OF LIABILITY NOTICE:
The United States Government shall not be liable or responsible for any maintenance,
updating or for correction of any errors in the SOFTWARE or subsequent approved version
releases.


THE INTERSPEC (SOFTWARE) AND ANY OF ITS SUBSEQUENT VERSION
RELEASES, SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY WARRANTY OF
ANY KIND, EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT
NOT LIMITED TO, ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO
SPECIFICATIONS, ANY IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE, OR FREEDOM FROM INFRINGEMENT, ANY
WARRANTY THAT THE SOFTWARE WILL BE ERROR FREE, OR ANY WARRANTY
THAT THE DOCUMENTATION, IF PROVIDED, WILL CONFORM TO THE
SOFTWARE. IN NO EVENT SHALL THE UNITED STATES GOVERNMENT OR ITS
CONTRACTORS OR SUBCONTRACTORS BE LIABLE FOR ANY DAMAGES,
INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR
CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY
WAY CONNECTED WITH THE SOFTWARE OR ANY OTHER PROVIDED
DOCUMENTATION, WHETHER OR NOT BASED UPON WARRANTY, CONTRACT,
TORT, OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY
PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS
SUSTAINED FROM, OR AROSE OUT OF THE RESULT OF, OR USE OF, THE
SOFTWARE OR ANY PROVIDED DOCUMENTATION. THE UNITED STATES
GOVERNMENT DISCLAIMS ALL WARRANTIES AND LIABILITIES REGARDING
THIRD PARTY SOFTWARE, IF PRESENT IN THE SOFTWARE, AND DISTRIBUTES
IT "AS IS."
```

## Acknowledgement
This InterSpec Software was developed with funds from the Science and Technology Directorate of the Department of Homeland Security.

InterSpec SCR# 2173.1
