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

/** Print some basic information about a nuclide. */
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
   This is useful if you know the current parent nuclides activity and when it
   was produced, and would like to know what gammas+xrays if is producing.
 */
void print_photons_of_aged_sample( const std::string &nuclideName,
                                  const double currentParentActivity,
                                  const double parentAge );

/** Takes two nuclides that that have already aged, makes a mixture out of them,
   then prints out the photon yield of the nuclides after aging for an
   additional amount of time.
 */
void example_nuclide_mixture_decay( const std::string &nuc1Name,
                                    const double nuc1_initial_age,
                                    const double nuc1_initial_activity,
                                    const std::string &nuc2Name,
                                    const double nuc2_initial_age,
                                    const double nuc2_initial_activity,
                                    const double timeToAgeMixture );

/** Prints the gammas produced in coincidence with each other during decays.
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
                         const bool allow_aging );

/** Calculates the expected number of gammas during a measurement, taking into account
 decay of the nuclides during the measurement.
 */
void example_correct_for_decays_during_measurements();

/** A simple example of building up nuclides, for example, during
 neutron irradiation.
 
 Uses a simple niave integration to track nuclides created, as well
 as their decays and progeny.
 
 Using this technique is advantageous of niave analytical solution
 when some produced nuclides decay through long chains, possibly
 including through other produced nuclides, and you want to track
 everything.
 */
void example_buildup_calculation();

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
  
  //Example
  example_correct_for_decays_during_measurements();
  
  example_buildup_calculation();
  
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


void example_correct_for_decays_during_measurements()
{
  cout << "Will demonstrate correcting for a nuclides decay during a measurement." << endl;
  
  using namespace SandiaDecay;
  SandiaDecayDataBase database( g_xml_filename );

  {// begin In110 example
    const SandiaDecay::Nuclide * const nuclide = database.nuclide("In110"); //(hl=4.9h), over a 2.8 hour measurement,
    const double activity_at_meas_start = 1.0*SandiaDecay::curie;
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide, activity_at_meas_start );
    
    const double initial_age = 0.0;
    const double measurement_duration = 2.8*SandiaDecay::hour;
    const size_t num_timeslices = 500;
    const vector<SandiaDecay::EnergyCountPair> photons
      = mix.decayPhotonsInInterval( initial_age, measurement_duration,
                                       SandiaDecay::NuclideMixture::OrderByEnergy, num_timeslices );
    
    cout << "During a " << (measurement_duration/SandiaDecay::hour) << " hour measurement of "
    << nuclide->symbol << " with an activity at the start of measurement of "
    << activity_at_meas_start/SandiaDecay::curie << " ci, the number of photons emitted will be:"
    << endl;
    for( size_t i = 0; i < photons.size(); ++i )
    {
      cout << "\t" << setw(10) << photons[i].energy << " keV: "
      << photons[i].count << endl;
    }
  }// end In110 example
  
  
  {// begin Mn56 example
    const SandiaDecay::Nuclide * const nuclide = database.nuclide("Mn56"); //(hl=9283.8s)
    const double activity_at_meas_start = 1.0*SandiaDecay::curie;
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide, activity_at_meas_start );
    
    const double initial_age = 0.0;
    const double measurement_duration = 86423.86*SandiaDecay::second;
    const size_t num_timeslices = 500;
    const vector<SandiaDecay::EnergyCountPair> photons
      = mix.decayPhotonsInInterval( initial_age, measurement_duration,
                                       SandiaDecay::NuclideMixture::OrderByEnergy, num_timeslices );
  }// end Mn56 example
  
  {// begin U235 example
    const SandiaDecay::Nuclide * const nuclide = database.nuclide("U235");
    const double activity_at_meas_start = 1.0*SandiaDecay::curie;
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide, activity_at_meas_start );
    
    const double initial_age = 10.0*SandiaDecay::year;
    const double measurement_duration = 1.0*SandiaDecay::year;
    const size_t num_timeslices = 500;
    const vector<SandiaDecay::EnergyCountPair> photons
      = mix.decayPhotonsInInterval( initial_age, measurement_duration,
                                       SandiaDecay::NuclideMixture::OrderByEnergy, num_timeslices );
  }// end U235 example
  
  {// begin Tc99 example
    const SandiaDecay::Nuclide * const nuclide = database.nuclide("Tc99");
    const double activity_at_meas_start = 1.0*SandiaDecay::curie;
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide, activity_at_meas_start );
    
    const double initial_age = 3600*SandiaDecay::second;
    const double measurement_duration = 1.0*SandiaDecay::year;
    const size_t num_timeslices = 500;
    const vector<SandiaDecay::EnergyCountPair> photons
      = mix.decayPhotonsInInterval( initial_age, measurement_duration,
                                       SandiaDecay::NuclideMixture::OrderByEnergy, num_timeslices );
  }// end Tc99 example
  
  {// begin Tc99 example
    const SandiaDecay::Nuclide * const nuclide = database.nuclide("Tc99");
    const double activity_at_meas_start = 1.0*SandiaDecay::curie;
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide, activity_at_meas_start );
    
    const double initial_age = 3600*SandiaDecay::second;
    const double measurement_duration = 1.0*SandiaDecay::year;
    const size_t num_timeslices = 500;
    const vector<SandiaDecay::EnergyCountPair> photons
      = mix.decayPhotonsInInterval( initial_age, measurement_duration,
                                       SandiaDecay::NuclideMixture::OrderByEnergy, num_timeslices );
  }// end Tc99 example
  
  {// begin I125 example
    const SandiaDecay::Nuclide * const nuclide = database.nuclide("I125"); //t12=12.93d
    const double activity_at_meas_start = 1.0*SandiaDecay::curie;
    
    SandiaDecay::NuclideMixture mix;
    mix.addNuclideByActivity( nuclide, activity_at_meas_start );
    
    const double initial_age = 0;
    const double measurement_duration = 14*24*3600.0*SandiaDecay::year;
    const size_t num_timeslices = 500;
    const vector<SandiaDecay::EnergyCountPair> photons
      = mix.decayPhotonsInInterval( initial_age, measurement_duration,
                                       SandiaDecay::NuclideMixture::OrderByEnergy, num_timeslices );
  }// end I125 example
}//void example_correct_for_decays_during_measurements()


void example_buildup_calculation()
{
  cout << "Demonstrating calculating resulting nuclides from buildup during neutron activation.\n";
  using namespace SandiaDecay;
  SandiaDecayDataBase database( g_xml_filename );
  
  
  // Define how long the sample will be irradiated
  const double irradiation_time_seconds = 7*24*3600; //i.e., 1-week
  
  // Define how long after irradiation it will be before a measurement is started
  const double cool_off_time = 2*3600; //i.e., 2 hours
  
  // Define how many seconds each time-step should be.
  //  A smaller time step should be more accurate, but you should take into account the
  //  half-lives of the nuclides you care about are.  For example, below Cr51 (HL=27.7d)
  //  dictates the time delta we really care about.  Ni59 is too long (practically doesnt
  //  decay on our relevant time-scale), and Fe55 is too short (doesnt build-up on our
  //  timescale), but for Cr51 we will start seeing some small numerical accuracy issues
  //  if we make the time delta more than a few hours (but the computation is pretty fast,
  //  so we'll just use a minute).
  const double time_delta = 60;  //i.e., 60 seconds.
  
  
  // Define nuclides that will be built up, and their buildup rate, in atoms per second
  vector<pair<const SandiaDecay::Nuclide *,double> > nuclides_rates;
  
  // We'll only add a few nuclides being produced here, but you can add as many
  //  as you'de like.
  //  The below indicated our experiment is expected to generate 63.8 Ni59 atoms,
  //  per second, from neutron transmutation; if there are any decays to Ni59 from
  //  other activation products, we dont need to account for them here, they will
  //  taken care of accumulating during decay calculations.
  nuclides_rates.push_back( make_pair(database.nuclide("Ni59"), 63.8) );  // HL=7.6E4 years
  nuclides_rates.push_back( make_pair(database.nuclide("Cr51"), 29.1) );  // HL=27.7 days
  nuclides_rates.push_back( make_pair(database.nuclide("Fe55"), 24.0) );  // HL=0.15s
  
  // We will step through the build-up time frame, with what we want to know
  //  is how many, of each type of radioactive atom will be left at the end of
  //  buildup.  So we will define a `NuclideMixture` to calculate this for us
  //  at each time step.
  SandiaDecay::NuclideMixture ager;
  
  // Add our nuclides to it
  for( size_t i = 0; i < nuclides_rates.size(); ++i )
  {
    const SandiaDecay::Nuclide * const nuclide = nuclides_rates[i].first;
    const double rate = nuclides_rates[i].second;
    
    const double num_atoms = time_delta * rate;
    ager.addNuclideByAbundance( nuclide, num_atoms );
  }
  
  // Now we will crate a mixture that will represent, for its t=0, the end of
  //  build up, and after the looping over all the time steps, it will have
  //  all the information we want.
  SandiaDecay::NuclideMixture mixture;
  
  // This integration is very niave, and could be greatly improved.
  //  However, informal checks show for irradiation time of months, and common
  //  neutron activation products in metals, accuracy and numeric error didn't
  //  become notable issues (checked smaller time deltas, as well as using 128 bit,
  //  instead of 64 bit internal precisions in SandiaDecay, as well as exact
  //  expectations).
  for( double start_time = 0.0; start_time < irradiation_time_seconds; start_time += time_delta )
  {
    const double end_time = std::min( start_time + time_delta, irradiation_time_seconds );
    const double this_dt = end_time - start_time;
    const double mid_time = 0.5*(start_time + end_time);
    const double time_until_irad_end = irradiation_time_seconds - mid_time;
    
    // Get the number of atoms, for all activation products, and their progeny, we expect
    //  at the end of buildup time.
    const vector<SandiaDecay::NuclideNumAtomsPair> num_atoms = ager.numAtoms( time_until_irad_end );
    for( size_t index = 0; index < num_atoms.size(); ++index )
    {
      const SandiaDecay::NuclideNumAtomsPair &nuc_num = num_atoms[index];
      mixture.addNuclideByAbundance( nuc_num.nuclide, (this_dt / time_delta)*nuc_num.numAtoms );
    }
  }//for( loop over buildup time )
  
  
  const vector<NuclideActivityPair> irrad_end_activities = mixture.activity( 0.0 );
  cout << "At the end of irradiation, the activities are:\n";
  for( size_t index = 0; index < irrad_end_activities.size(); ++index )
    cout << "\t" << irrad_end_activities[index].nuclide->symbol
          << ": " << irrad_end_activities[index].activity << " bq" << endl;
  
  const vector<NuclideActivityPair> after_cool_off_activities = mixture.activity( cool_off_time );
  cout << "\n\nAfter cooling off for " << cool_off_time << " seconds the activities are:\n";
  for( size_t index = 0; index < after_cool_off_activities.size(); ++index )
    cout << "\t" << after_cool_off_activities[index].nuclide->symbol
    << ": " << after_cool_off_activities[index].activity << " bq" << endl;
  
  
  // We expect A = A_0 * (1 - exp(-lamda * t_activation), so lets check things, but
  //  please note that this is only a valid check if no other activation products
  //  decay through the activation nuclide of interest.
  cout << endl << endl;
  for( size_t index = 0; index < irrad_end_activities.size(); ++index )
  {
    const SandiaDecay::Nuclide * const output_nuc = irrad_end_activities[index].nuclide;
    const double out_act = irrad_end_activities[index].activity;
    for( size_t input_index = 0; input_index < nuclides_rates.size(); ++input_index )
    {
      const SandiaDecay::Nuclide * const input_nuc = nuclides_rates[input_index].first;
      const double input_rate = nuclides_rates[input_index].second;
      if( input_nuc != output_nuc )
        continue;
      
      const double lambda = input_nuc->decayConstant();
      const double expected_act = input_rate * (1.0 - exp( -lambda * irradiation_time_seconds) );
      cout << "For " << input_nuc->symbol << " analytically expected " << expected_act
            << " bq; our calculation is " << out_act << " bq" << endl;
    }//for( loop over input nuclides )
  }//for( loop over output nuclides )
  
}//void example_buildup_calculation()
