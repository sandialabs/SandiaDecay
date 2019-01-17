/* SandiaDecay: a library that provides nuclear decay info and calculations.
 
 Copyright 2018 National Technology & Engineering Solutions of Sandia, LLC
 (NTESS). Under the terms of Contract DE-NA0003525 with NTESS, the U.S.
 Government retains certain rights in this software.
 For questions contact William Johnson via email at wcjohns@sandia.gov, or
 alternative email of interspec@sandia.gov.
 
 This library is free software; you can redistribute it and/or
 modify it under the terms of the GNU Lesser General Public
 License as published by the Free Software Foundation; either
 version 2.1 of the License, or (at your option) any later version.
 
 This library is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 Lesser General Public License for more details.
 
 You should have received a copy of the GNU Lesser General Public
 License along with this library; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include <cmath>
#include <time.h>
#include <string>
#include <iostream>
#include <algorithm>

#include "SandiaDecay.h"

using namespace std;
using namespace SandiaDecay;

//decay_database_benchmark(): creates the database, and decays each nuclide
//  (at 5 seperate times), 10 times.  Prints out time to create the database,
//  and time to do the decay, for each trial, and then the best of the 10 at
//  the end.
//Marks from my '2.7 GHz Intel Core i7, 16 GB 1600 MHz DDR3, OS X 10.9.2'
//  computer, when the code is compiled as MinSizeRel as part of Interspec:
//--Code pre 20140426 (no multithread, doubles everywhere, independant
//  allocations of Nuclide and Elements objects):
//    Database takes up approx 12086238 bytes of memorry (11803 kb, 11.5263 MB)
//    Best trial: DB Open=79 ms, Decay=452 ms
//--Code after converting to mostly storing floats instead of doubles, and
//  short ints instead of ints (left half life as double).
//    Database takes up approx 9522918 bytes of memorry (9299.72 kb, 9.08176 MB)
//    Best trial: DB Open=79 ms, Decay=453 ms
//--Code after making it so Nuclides all stored in a single contiguos array
//  of memmory
//    Database takes up approx 9517966 bytes of memorry (9294.89 kb, 9.07704 MB)
//    Best trial: DB Open=0079 ms, Decay=0456 ms
//--After implemnting a niave mutlithreaded approach to initializing Nuclide,
//  Element, and Transitions, using 4 parrelel threads:
//    Database takes up approx 9517990 bytes of memorry (9294.91 kb, 9.07706 MB)
//    Best trial: DB Open=0085 ms, Decay=0459 ms
//--After going back to single threaded and removing some stringstream
//  operations:
//    Database takes up approx 9517990 bytes of memorry (9294.91 kb, 9.07706 MB)
//    Best trial: DB Open=0074 ms, Decay=0447 ms
void decay_database_benchmark( const char *nuclidedb );

int main( int argc, const char * argv[] )
{
  if( argc < 2 )
  {
    cerr << "Use: " << argv[0] << " <path to sandia.decay.xml>" << endl;
    return EXIT_FAILURE;
  }
  
  const char *dbfile = argv[1];
 
  decay_database_benchmark( dbfile );

  return 1;
}//int main( int argc, const char * argv[] )


bool has_decendant( const Nuclide *nuclide, string desc)
{
  const vector<const Nuclide *> kids = nuclide->descendants();
  for( size_t i = 0; i < kids.size(); ++i )
  {
    if( kids[i]->symbol == desc )
      return true;
  }
  return false;
};


void decay_database_benchmark( const char *nuclidedb )
{
  clock_t mindecay = std::numeric_limits<clock_t>::max();
  clock_t mindbopen = std::numeric_limits<clock_t>::max();
  
  for( int trial = 0; trial < 10; ++trial )
  {
    const clock_t dbopen_start = clock();
    
    SandiaDecayDataBase database( nuclidedb );
    
    const clock_t decay_start = clock();
    
    const vector<const Nuclide *> &nuclidesVec = database.nuclides();
    
    size_t num_nuclieds_decayed = 0;
    vector<const Nuclide *>::const_iterator iter;
    for( iter = nuclidesVec.begin(); iter != nuclidesVec.end(); ++iter )
    {
      const Nuclide *nuclide = *iter;
      const double halfLife = nuclide->halfLife;

      if( std::isinf(halfLife) )
        continue;
      
      ++num_nuclieds_decayed;
      NuclideMixture mixture;
      mixture.addNuclideByActivity( nuclide, 1.0*SandiaDecay::curie );
      for( size_t i = 0; i < 5; ++i )
      {
        const double time = i*7.0*halfLife/5.0;
        const vector<NuclideActivityPair> activities = mixture.activity( time );
        const vector<EnergyRatePair> gammas = mixture.gammas( time, NuclideMixture::OrderByAbundance, true );
        
        //feeble attempt to thowart compiler optimizations
        if( activities.size() > 9999999 || gammas.size() > 9999999 )
          cerr << "A dummy statment" << endl;  //wont ever printout
      }//for( test the decay for 7 halfLives )
    }//for( iterate over nuclides available, iter )
    
    const clock_t decay_done = clock();
    const clock_t opentime = decay_start - dbopen_start;
    const clock_t decaytime = decay_done - dbopen_start;
    mindecay = std::min( mindecay, decaytime );
    mindbopen = std::min( mindbopen, opentime );
    
    printf( "Trial %i: DB Open=%.4i ms, Decay=%.4i ms (%.3f us per nuclide deced, "
            "got gammas and activities for 5 different times for each decay)\n",
            trial, int(opentime*1000/CLOCKS_PER_SEC),
            int(decaytime*1000/CLOCKS_PER_SEC),
            double(decaytime*1000000.0/(CLOCKS_PER_SEC*num_nuclieds_decayed)));
  }//for( int trial = 0; trial < 10; ++trial )
  
  SandiaDecayDataBase database( nuclidedb );

  const size_t mem = database.memsize();
  std::cout << "Database takes up approx " << mem << " bytes of memorry ("
            << mem/1024.0 << " kb, " << mem/1048576.0 << " MB)\n" << std::endl;
  
  printf( "Best trial: DB Open=%.4i ms, Decay=%.4i ms\n",
          int(mindbopen*1000/CLOCKS_PER_SEC),
          int(mindecay*1000/CLOCKS_PER_SEC) );
}//decay_database_benchmark
