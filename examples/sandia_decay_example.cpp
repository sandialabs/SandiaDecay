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
#include <vector>
#include <string>
#include <iomanip>
#include <sstream>
#include <float.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <assert.h>
#include <algorithm>

#include "SandiaDecay.h"

using namespace std;

//To compile on UNIX like systems:
//  g++ sandia_decay_example.cpp SandiaDecay.cpp -o sandia_decay.exe
//
//Forward declarations
//

/** Print some basic infomration about a nuclide. */
void print_nuclide_info( const std::string &nuclideName );

/** Print the decay chain of the specified nuclide. */
void print_nuclide_decay_chain( const std::string &nuclideName );

/** Prints out the activity of the parent and decendant nuclides after the
   specified amount of time.
 */
void simple_decay_print_descendant_activities( const std::string &nuclideName,
                                               const double initialParentActivity,
                                               const double timeToAge );

/** Prints out the gamma and x-rays priduced the parent and its decendant
   nuclides after the specified amount of time.
 */
void simple_decay_print_photon_production( const std::string &nuclideName,
                                           const double initialParentActivity,
                                           const double timeToAge );

/** Prints the photons (gammas + xrays) produced by a nuclide with the specified
   age, and the specified current activity.
   This is useuful if you know the current parent nuclides activity and when it
   was produced, and would like to know what gammas+xrays if is producing.
 */
void print_photons_of_aged_sample( const std::string &nuclideName,
                                  const double currentParentActivity,
                                  const double parentAge );

/** Takes two nuclides that that have already aged, makes a mixture out of them,
   then prints out the photon yeild of the nuclides after aging for an
   additional amount of time.
 */
void example_nuclide_mixture_decay( const std::string &nuc1Name,
                                    const double nuc1_initial_age,
                                    const double nuc1_initial_activity,
                                    const std::string &nuc2Name,
                                    const double nuc2_initial_age,
                                    const double nuc2_initial_activity,
                                    const double timeToAgeMixture );

/** Prints the gammas produced in coincidence with eachother during decays.
 Only includes gammas of nuclide specified, and not any of its descendants.
 */
void print_coincident_gammas_of_nuclide( const std::string &nuclideName );

/** Prints out the nuclides that have a gamma in the specified energy range.
 Energies are in keV, and if allowAging is specified, then each nuclide will
 be aged by one half life before determining if it contains a gamma in range.
 */
void print_nuclides_with_gammas_in_energy_range( const double lower_energy,
                                                 const double upper_energy,
                                                 const bool allowAging );

/** Helper function for print_nuclides_with_gammas_in_energy_range(...) */
std::vector<const SandiaDecay::Nuclide *>
nuclidesWithGammaInRange( const float lowE, const float highE,
                         const std::vector<const SandiaDecay::Nuclide *> &candidates,
                         const bool allowaging );


std::string g_xml_filename = "sandia.decay.xml";


int main( int argc, const char * argv[] )
{

  if( argc > 1 )
    g_xml_filename = argv[1];

  cout << "Using '" << g_xml_filename << "' as input XML file." << endl;
  
  //Example of printing the activity of a nuclide and its descendants after some
  //  aging.
  string nuclide = "133Ba";
  double timeToAge = 100 * 24 * 3600 * SandiaDecay::second;  //100 days
  double initialParentActivity = 0.001 * SandiaDecay::curie; //1 mCi
  simple_decay_print_descendant_activities( nuclide, initialParentActivity, timeToAge );
  
  
  //Print the photon (xrays+gammas) energies and rates produced by a nuclide
  //  with the given initial activity, after it ages for the specified time.
  simple_decay_print_photon_production( nuclide, initialParentActivity, timeToAge );
  
  
  //Print the photon (xrays+gammas) energies and rates of a aged nuclide based
  //  on the parent isotopes current activity.
  print_photons_of_aged_sample( "238-U", 1.0E6*SandiaDecay::Bq, 20*SandiaDecay::year );
  
  
  //Example of aging a mixture of nuclides
  example_nuclide_mixture_decay( "60co", 1.0*SandiaDecay::year, 1.0E-6*SandiaDecay::Ci,
                                 "133 ba", 0.0*SandiaDecay::second, 1.0E-4*SandiaDecay::Ci,
                                 0.5*SandiaDecay::year );

  //Print the decay chain and branching ratios of a nuclide
  print_nuclide_decay_chain( "U-238" );

  
  //Print the coincident gammas of a nuclide
  print_coincident_gammas_of_nuclide( "60-Co" );
  
  
  //Print some basic information about a nuclide
  print_nuclide_info( "Co60" );
  
  
  //Find nuclides with gammas in a specified energy range.
  print_nuclides_with_gammas_in_energy_range( 1332.48, 1332.50, false );
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


void simple_decay_print_photon_production( const std::string &nuclideName,
                                              const double initialParentActivity,
                                              const double timeToAge )
{
  SandiaDecay::SandiaDecayDataBase database( g_xml_filename );
  
  const SandiaDecay::Nuclide * const nuclide = database.nuclide( nuclideName );
  
  if( !nuclide )
  {
    cerr << "'" << nuclideName << "' is not a valid nuclide - sorry\n";
    return;
  }
  
  assert( initialParentActivity > 0.0 );
  assert( timeToAge > 0.0 && timeToAge < 10*nuclide->halfLife );
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( nuclide, initialParentActivity );
  
  //Decide whether to include the 511 keV photons from positron annihilations
  //  included with the gammas.  If you want to get the the 511 keVs seperately
  //  call NuclideMixture::betaPlusses(timeToAge, ordering).
  const bool includeAnnihilations = true;
  
  //Decide how you want results ordered, either OrderByEnergy, or OrderByAbundance
  const SandiaDecay::NuclideMixture::HowToOrder ordering = SandiaDecay::NuclideMixture::OrderByEnergy;
  
  //Get the number of gammas you would expect to see, per second, at the desired age
  const vector<SandiaDecay::EnergyRatePair> gammas = mixture.gammas( timeToAge, ordering, includeAnnihilations );
  
  //Get the xrays produced by the decays that you would expect, per second, at the desired age
  const vector<SandiaDecay::EnergyRatePair> xrays = mixture.xrays( timeToAge, ordering );
  
  //If you just care about what your radiation detector might see, you could
  //  just call NuclideMixture::photons(timeToAge,ordering) to get all of the
  //  x-rays, gammas, and annihilation gammas mixed together
  //vector<SandiaDecay::EnergyRatePair> photons = mixture.photons( timeToAge, ordering );
  
  
  cout << "After aging " << timeToAge/SandiaDecay::year << " years, "
       << nuclide->symbol << " will produce:\n";
  
  for( size_t i = 0; i < xrays.size(); ++i )
  {
    const SandiaDecay::EnergyRatePair &energyRate = xrays[i];
    cout << "\txray " << energyRate.energy << " keV: " << energyRate.numPerSecond << "/second\n";
  }
  
  for( size_t i = 0; i < gammas.size(); ++i )
  {
    const SandiaDecay::EnergyRatePair &energyRate = gammas[i];
    cout << "\tgamma " << energyRate.energy << " keV: " << energyRate.numPerSecond << "/second\n";
  }
  
  cout << endl;
  
  //Note: do not delete or free `nuclide` as the object it points to is owned
  //      by SandiaDecayDataBase and will be cleaned up when `database`
  //      desctructs.
}//void simple_decay_print_photons( const std::string &nuclideName )


void print_photons_of_aged_sample( const std::string &nuclideName,
                                  const double currentParentActivity,
                                  const double parentAge )
{
  SandiaDecay::SandiaDecayDataBase database( g_xml_filename );
  
  const SandiaDecay::Nuclide * const nuclide = database.nuclide( nuclideName );
  
  if( !nuclide )
  {
    cerr << "'" << nuclideName << "' is not a valid nuclide - sorry\n";
    return;
  }
  
  assert( parentAge >= 0.0 );
  assert( currentParentActivity > 0.0 );
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addAgedNuclideByActivity( nuclide, currentParentActivity, parentAge );
  
  const double timeToAge = 0.0;
  const SandiaDecay::NuclideMixture::HowToOrder ordering = SandiaDecay::NuclideMixture::OrderByEnergy;
  
  const vector<SandiaDecay::EnergyRatePair> photons = mixture.photons( timeToAge, ordering );
  
  
  cout << currentParentActivity/SandiaDecay::Bq << " bq of " << nuclide->symbol
  << " that is " << parentAge/SandiaDecay::year << " years old will produce"
  << " the following photons:\n";
  
  for( size_t i = 0; i < photons.size(); ++i )
  {
    cout << "\t" << photons[i].energy << " keV: " << photons[i].numPerSecond << "/second\n";
    if( i > 35 )  //keep from printing out like a million energies.
    {
      cout << "... skipping " << (photons.size()-i) << " energies ..." << endl;
      break;
    }
  }
  cout << endl;
  
  //Note: do not delete or free `nuclide` as the object it points to is owned
  //      by SandiaDecayDataBase and will be cleaned up when `database`
  //      desctructs.
}//print_photons_of_aged_sample(...)



void simple_decay_print_descendant_activities( const std::string &nuclideName,
                                              const double initialParentActivity,
                                              const double timeToAge )
{
  SandiaDecay::SandiaDecayDataBase database( g_xml_filename );
  
  const SandiaDecay::Nuclide * const nuclide = database.nuclide( nuclideName );
  
  if( !nuclide )
  {
    cerr << "'" << nuclideName << "' is not a valid nuclide - sorry\n";
    return;
  }
  
  assert( initialParentActivity > 0.0 );
  assert( timeToAge > 0.0 && timeToAge < 10*nuclide->halfLife );
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addNuclideByActivity( nuclide, initialParentActivity );
  
  const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( timeToAge );
  
  cout << "After " << timeToAge/SandiaDecay::second << " seconds, "
  << nuclide->symbol << " and its descendants will have activity:\n";
  for( size_t i = 0; i < activities.size(); ++i )
  {
    const SandiaDecay::NuclideActivityPair &nucAct = activities[i];
    cout << "\t" << nucAct.nuclide->symbol << ": " << nucAct.activity/SandiaDecay::curie << " curie"
    << (nucAct.nuclide->isStable() ? " (stable)" : "") << endl;
  }
  
  cout << endl;
  
  //Note: do not delete or free `nuclide` as the object it points to is owned
  //      by SandiaDecayDataBase and will be cleaned up when `database`
  //      desctructs.
}//void simple_decay_print_photons( const std::string &nuclideName )


void example_nuclide_mixture_decay( const std::string &nuc1Name,
                                   const double nuc1_initial_age,
                                   const double nuc1_initial_activity,
                                   const std::string &nuc2Name,
                                   const double nuc2_initial_age,
                                   const double nuc2_initial_activity,
                                   const double timeToAgeMixture )
{
  assert( nuc1_initial_age >= 0.0 );
  assert( nuc2_initial_age >= 0.0 );
  assert( timeToAgeMixture >= 0.0 );
  
  SandiaDecay::SandiaDecayDataBase database( g_xml_filename );
  
  const SandiaDecay::Nuclide * const nuc1 = database.nuclide( nuc1Name );
  const SandiaDecay::Nuclide * const nuc2 = database.nuclide( nuc2Name );
  
  assert( nuc1 );
  assert( nuc2 );
  
  SandiaDecay::NuclideMixture mixture;
  mixture.addAgedNuclideByActivity( nuc1, nuc1_initial_activity, nuc1_initial_age );
  mixture.addAgedNuclideByActivity( nuc2, nuc2_initial_activity, nuc2_initial_age );
  
  cout << "After " << timeToAgeMixture/SandiaDecay::second << " seconds, "
  << nuc1->symbol << " (initially " << nuc1_initial_age << " seconds old with activity "
  << nuc1_initial_activity/SandiaDecay::Ci << " curie) and "
  << nuc2->symbol << " (initially " << nuc2_initial_age << " seconds old with activity "
  << nuc2_initial_activity/SandiaDecay::Ci << " curies)"
  << " and their descendants will produce the following photons (gammas+xrays):\n";
  
  //const vector<SandiaDecay::NuclideActivityPair> activities = mixture.activity( timeToAgeMixture );
  const SandiaDecay::NuclideMixture::HowToOrder ordering = SandiaDecay::NuclideMixture::OrderByEnergy;
  const vector<SandiaDecay::EnergyRatePair> photons = mixture.photons( timeToAgeMixture, ordering );
  
  for( size_t i = 0; i < photons.size(); ++i )
  {
    const SandiaDecay::EnergyRatePair &energyRate = photons[i];
    cout << "\t" << energyRate.energy << " kev, " << energyRate.numPerSecond << "/second\n";
  }
  cout << endl;
  
  //Note: do not delete or free `nuc1` or `nuc2` as the object it points to is
  //      owned by SandiaDecayDataBase and will be cleaned up when `database`
  //      desctructs.
}//void example_nuclide_mixture_decay(...)


void print_nuclide_info( const string &nuclideName )
{
  using namespace SandiaDecay;
  SandiaDecayDataBase database( g_xml_filename );

  const Nuclide * const nuclide = database.nuclide( nuclideName );

  if( !nuclide )
  {
    cerr << "Could not find isotope '" << nuclideName << "' in the database" << endl;
    return;
  }
  
  const vector<const Transition *> &decaysToChildren = nuclide->decaysToChildren;
  const vector<const Transition *> &decaysFromParents = nuclide->decaysFromParents;
  
  cout << nuclide->symbol
       << " Atomic Number " << nuclide->atomicNumber
       <<", Atomic Mass " << nuclide->massNumber
       << ", Isomer Number " << nuclide->isomerNumber << " "
       << nuclide->atomicMass << " AMU, HalfLife=" << nuclide->halfLife << " seconds";
  
  const size_t nParents = decaysFromParents.size();
  if( nParents )
    cout << "\n  Parent" << (nParents > 1?"s:":":") << (nParents?"":" none");
  
  for( size_t i = 0; i < nParents; ++i )
  {
    if( decaysFromParents[i]->parent )
      cout << (i?", ":"") << nuclide->decaysFromParents[i]->parent->symbol;
  }//for( loop over parents )
  
  const size_t nChilds = nuclide->decaysToChildren.size();
  
  for( size_t i = 0; i < nChilds; ++i )
  {
    cout << "\n    " << human_str_summary(*(decaysToChildren[i]));
  }

  //Instead of doing all the above you could have just called:
  //cout << SandiaDecay::human_str_summary(*nuclide) << endl;
  
  cout << endl << endl;
  
  //Note: do not delete or free `nuclide` or anything else returned from
  //      `database` as they are all owned by SandiaDecayDataBase and will be
  //      cleaned up when `database` desctructs.
}//int print_nuclide_info(int argc, const char * argv[])



void print_nuclide_decay_chain( const string &nuclideName )
{
  using namespace SandiaDecay;

  SandiaDecayDataBase database( g_xml_filename );

  const Nuclide * const nuclide = database.nuclide( nuclideName );

  if( !nuclide )
  {
    cerr << "Could not find isotope '" + nuclideName + "' in the database" << endl;
    return;
  }

  stringstream output; //use a stringstream so we dont mess up formatting in other functions
  output << setw(6) << nuclide->symbol << " decays through (% of time):" << endl;
  const vector<const Nuclide *> decendants = nuclide->descendants();

  for( size_t i = 0; i < decendants.size(); ++i )
  {
    const double br = nuclide->branchRatioToDecendant( decendants[i] );
    output << "    " << setw(6) << left << decendants[i]->symbol << " ("
         << fixed << setprecision(2) << setw(6) << right << 100.0*br
         << "%)" << endl;
  }//for( loop over decendants, i )

  cout << output.str() << endl;
  
  //Note: do not delete or free `nuclide` or anything else returned from
  //      `database` as they are all owned by SandiaDecayDataBase and will be
  //      cleaned up when `database` desctructs.
}//int print_nuclide_decay_chain( const char *nuclideName );


void print_coincident_gammas_of_nuclide( const std::string &nuclideName )
{
  SandiaDecay::SandiaDecayDataBase database( g_xml_filename );
  
  const SandiaDecay::Nuclide * const nuclide = database.nuclide( nuclideName );
  
  if( !nuclide )
    throw runtime_error( "Could not find isotope '" + nuclideName + "' in the database" );
  
  cout << nuclide->symbol << " has coincident gammas (does not include descendant nuclides):\n";
  
  for( size_t trans_index = 0; trans_index < nuclide->decaysToChildren.size(); ++trans_index )
  {
    const SandiaDecay::Transition * const trans = nuclide->decaysToChildren[trans_index];
    
    for( size_t prod_index = 0; prod_index < trans->products.size(); ++prod_index )
    {
      const SandiaDecay::RadParticle &particle = trans->products[prod_index];
      
      if( particle.coincidences.empty() )
        continue;
      
      cout << "\t" << particle.energy << " keV " << SandiaDecay::to_str(particle.type)
      << " (br=" << particle.intensity << ") coincident with:\n";
      for( size_t coinc_index = 0; coinc_index < particle.coincidences.size(); ++coinc_index )
      {
        const unsigned short int part_ind = particle.coincidences[coinc_index].first;
        const float fraction = particle.coincidences[coinc_index].second;
        const SandiaDecay::RadParticle &coinc_part = trans->products[part_ind];
        
        cout << "\t\t" << coinc_part.energy << " keV " << SandiaDecay::to_str(coinc_part.type)
             << " at " << fraction << endl;
      }//for( loop over coincidences )
    }//for( loop over transitions )
  }//for( loop over decay modes )
  
  cout << endl;
  
  //Note: do not delete or free `nuclide` or anything else returned from
  //      `database` as they are all owned by SandiaDecayDataBase and will be
  //      cleaned up when `database` desctructs.
}//void print_coincident_gammas_of_nuclide( const std::string &nuclideName );


void print_nuclides_with_gammas_in_energy_range( const double lower_energy,
                                                 const double upper_energy,
                                                 const bool allowaging )
{
  assert( upper_energy > lower_energy );
  
  using namespace SandiaDecay;
  SandiaDecayDataBase database( g_xml_filename );

  const vector<const Nuclide *> &all_nuclides = database.nuclides();

  vector<const Nuclide *> nuclides;
  nuclides = nuclidesWithGammaInRange( lower_energy, upper_energy,
                                       all_nuclides, allowaging );

  cout << "Nuclides with gammas in the range " << lower_energy/keV << " keV"
       << " to " << upper_energy/keV << " keV:\n\t";

  for( size_t i = 0; i < nuclides.size(); ++i )
    cout << (i ? ", " : "") << nuclides[i]->symbol;

  cout << endl << endl;
  
  //Note: do not delete or free the elements in `nuclides` as they are all owned
  //      by SandiaDecayDataBase and will be cleaned up when `database` desctructs.
}//int print_nuclides_with_gammas_in_energy_range()

std::vector<const SandiaDecay::Nuclide *>
nuclidesWithGammaInRange( const float lowE, const float highE,
                          const std::vector<const SandiaDecay::Nuclide *> &candidates,
                          const bool allowaging )
{
  using namespace SandiaDecay;
  
  std::vector<const Nuclide *> answer;
  
  if( allowaging )
  {
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      const Nuclide * const nuc = candidates[i];
      NuclideMixture mix;
      mix.addNuclideByActivity( nuc, 1.0E6 * becquerel );
      const vector<EnergyRatePair> gammas = mix.gammas( nuc->halfLife, NuclideMixture::OrderByAbundance, true );
      
      for( size_t i = 0; i < gammas.size(); ++i )
      {
        const EnergyRatePair &a = gammas[i];
        if( a.energy >= lowE && a.energy <= highE )
          if( find(answer.begin(), answer.end(), nuc) == answer.end() )
            answer.push_back( nuc );
      }//for( size_t i = 0; i < gammas.size(); ++i )
    }//for( size_t i = 0; i < candidates.size(); ++i )
  }else
  {
    for( size_t i = 0; i < candidates.size(); ++i )
    {
      const Nuclide * const nuc = candidates[i];
      const std::vector<const Transition *> &trans = nuc->decaysToChildren;
      for( size_t j = 0; j < trans.size(); ++j )
      {
        const Transition * const tran = trans[j];
        const std::vector<RadParticle> &particles = tran->products;
        for( size_t k = 0; k < particles.size(); ++k )
        {
          const RadParticle &particle = particles[k];
          if( particle.type != GammaParticle )
            continue;
          
          if( particle.energy >= lowE && particle.energy <= highE )
            if( find(answer.begin(),answer.end(),nuc) == answer.end() )
              answer.push_back( nuc );
        }//for( size_t k = 0; k < particles.size(); ++k )
      }//for( size_t j = 0; j < trans.size(); ++j )
    }//for( size_t i = 0; i < candidates.size(); ++i )
  }//if( allowaging ) / else
  
  return answer;
}//nuclidesWithGammaInRange(...)


