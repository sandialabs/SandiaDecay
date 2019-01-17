#include <cmath>
#include <vector>
#include <string>
#include <iomanip>
#include <float.h>
#include <cstring>
#include <sstream>
#include <cstdlib>
#include <iostream>
#include <assert.h>

#include "SandiaDecay.h"


//To compile on UNIX like systems:
//  g++ sandia_decay.cpp SandiaDecay.cpp -o sandia_decay.exe

/* Further sanity checks that could be done:
   -Set DecayCoefType to __float128 in SandiaDecay.cpp and check a bunch of
    decays to see how the answers differ from using double
   -Should probably manually add in some decays where sandia.decay.xml is
    unlikely to change, and "unit-test" those decays.
   -Add in numerical stability tests for BRs of decay products
   -Check that all float fields in the file have valid float values (right now
    SandiaDecay.cpp uses atof() which does not have a good error mechansim)
   -The coincidence range check should be improved
 */

using namespace std;

//Forward declarations
void check_sum_branching_ratios();
void sanity_check_all_decays();
void sanity_check_nuclides_and_transistions();
void check_bateman_denominator();

//global variables - whatever, I do what I want!
string g_xml_file;


int main( int argc, const char * argv[] )
{ 
  if( argc < 2 )
  {
    cerr << "Usage: " << argv[0] << " path/to/sandia.decay.xml" << endl;
    return EXIT_FAILURE;
  }
  
  g_xml_file = argv[1];
  
  try
  {
    sanity_check_nuclides_and_transistions();
    sanity_check_all_decays();
    check_bateman_denominator();
    check_sum_branching_ratios();
  }catch( std::exception &e )
  {
    cerr << "Caught error doign tests: " << e.what() << endl;
    return EXIT_FAILURE;
  }//try . catch
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )


namespace
{
  //Makeup for not having c++11
  template <typename Pod>
  std::string to_str( Pod val )
  {
    std::stringstream strm;
    strm << val;
    return strm.str();
  }
  
  
  /** Checks for breakdown of Bateman equations when two nuclides in a decay
    chain have really close halflives.
   
   Note that we should probably also check TimeEvolutionTerm::eval()
   */
  void check_bateman_stability( const SandiaDecay::Nuclide * const chain_parent,
                         std::vector<const SandiaDecay::Transition *> decay_path )
  {
    const SandiaDecay::Transition *trans = decay_path.back();
    const SandiaDecay::Nuclide *child = trans->child;
    
    const size_t row = decay_path.size();
    for( size_t col = 0; col < row; ++col )
    {
      //I *think* the only quantities that get subtracted, were we might get
      //  into real trouble, are these next two.  Everything else is division,
      //  multiplication, or sums where loosing some terms due to presicion
      //  shouldnt matter much... I'm sure I'm missing something...
      const double lambda_i = decay_path.at(row-1)->child->decayConstant();
      const double lambda_j = decay_path.at(col)->parent->decayConstant();
      
      //Using the criteria that the difference between the lambdas has to be 1000
      //  times larger than the epsilon (the smallest change to the value that can
      //  be made and not have the result equal to the initial value) of doubles
      //  for the values of lamda.  The 1.0E06 is totally arbitrary - and this
      //  check is probably totally niavly or stupidly done - but its something!
      //This essentually checks for nuclides in a decay chain with the same HL
      const double lamda_diff = fabs(lambda_i - lambda_j);
      if( ((lambda_i + (lamda_diff/1.0E06)) == lambda_i)
         || ((lambda_j + (lamda_diff/1.0E06)) == lambda_j) )
      {
        stringstream msg;
        msg << setprecision(17) << "Got lambda_i=" << lambda_i
            << ", lambda_j=" << lambda_j
            << " is unstable for lamda_diff=" << lamda_diff
            << " for lamda_i->child=" << decay_path.at(row-1)->child->symbol
            << " and lambda_j->parent=" << decay_path.at(col)->parent->symbol
            << " (chain parent=" << chain_parent->symbol << ")";
        throw runtime_error( msg.str() );
      }
    }//for( loop over 'col' of matrix A )
    
    vector<const SandiaDecay::Transition *> decays;
    for( size_t t = 0; t < child->decaysToChildren.size(); ++t )
      if( child->decaysToChildren[t]->child )
        decays.push_back( child->decaysToChildren[t] );
    
    if( decays.empty() )
      return;
    
    const size_t ndecays = decays.size();
    for( size_t t = 0; t < ndecays; ++t )
    {
      vector<const SandiaDecay::Transition *> new_path = decay_path;
      new_path.push_back( decays[t] );
      check_bateman_stability( chain_parent, new_path );
    }//for( loop over decaysToChildren )
  }//void check_bateman_stability(...)
}


void check_bateman_denominator()
{
  //Checks that no decendant nuclides have nearly same half lives so the bateman
  //  equations wont break down.
  using namespace SandiaDecay;
  SandiaDecayDataBase database( g_xml_file );
  
  const vector<const Nuclide *> &nucs = database.nuclides();
  for( size_t i = 0; i < nucs.size(); ++i )
  {
    const Nuclide *const parent = nucs[i];
    
    for( size_t decayNum = 0; decayNum < parent->decaysToChildren.size(); ++decayNum )
    {
      const Transition *trans = parent->decaysToChildren[decayNum];
      if( trans->child )
        check_bateman_stability( parent, vector<const Transition *>(1, trans) );
    }//for( loop over transition, decayNum )
  }//for( size_t i = 0; i < nucs.size(); ++i )
}//void check_bateman_denominator()


void sanity_check_nuclides_and_transistions()
{
  using namespace SandiaDecay;
  SandiaDecayDataBase database( g_xml_file );
  
  const vector<const Nuclide *> &nucs = database.nuclides();
  for( size_t i = 0; i < nucs.size(); ++i )
  {
    const Nuclide *const nuc = nucs[i];
    if( !nuc )
      throw runtime_error( "Found nullptr nuclide - shouldnt happen!" );
    
    if( nuc->symbol.empty() )
      throw runtime_error( "Found nuclide with no symbol" );
    
    if( std::isinf(nuc->halfLife) && nuc->decaysToChildren.size() )
      throw runtime_error( "Nuclide " + nuc->symbol + " has infinite halflife (i.e., stable) but has "
                           + to_str(nuc->decaysToChildren.size()) + " children." );
    
    if( nuc->isStable() )
      continue;
    
    if( nuc->halfLife <= 0.0 )
      throw runtime_error( "Halflife for " + nuc->symbol + " is " + std::to_string(nuc->halfLife) );
    
    if( nuc->isomerNumber < 0 || nuc->isomerNumber > 5 )
      throw runtime_error( "Isomer number for " + nuc->symbol + " is invalid." );
    
    const Element * const el = database.element(nuc->atomicNumber);
    if( !el )
      throw runtime_error( "Could not get element (atomic number=" + to_str(nuc->atomicNumber) + ") for nuclide " + nuc->symbol );
    
    string elstr = nuc->symbol;
    const size_t numpos = elstr.find_first_of( "0123456789" );
    elstr = elstr.substr(0,numpos);
    
    const Element * const el1 = database.element(elstr);
    if( nuc->symbol!="n1" && el1 != el )  //"n1" gets a pass...
      throw runtime_error( "Getting element for " + nuc->symbol + " via label and atomic number did not get same element." );
    
    
    const std::vector<const Transition *> &transitions = nuc->decaysToChildren;
    for( size_t j = 0; j < transitions.size(); ++j )
    {
      const Transition * const trans = transitions[j];
      if( !trans )
        throw runtime_error( "Found invalid transition for nuclide " + nuc->symbol );
      
      if( trans->parent != nuc )
        throw runtime_error( "Found transistion for " + nuc->symbol + " whose parent wasnt as expected" );
      
      const string childstr = trans->child ? trans->child->symbol : string("SF");
      
      if( trans->branchRatio < 0.0 || trans->branchRatio > 1.0 )
        throw runtime_error( "Transition for " + nuc->symbol + " to " + childstr + " was " + to_str(trans->branchRatio) );
      
      if( trans->child == trans->parent )
        throw runtime_error( "Transistion for " + nuc->symbol + " has a decay to itself!" );
      
      if( trans->mode == SandiaDecay::SpontaneousFissionDecay )
      {
        if( trans->child )
          throw runtime_error( "Spontaneous Fission for " + nuc->symbol + " has a child specified" );
      }else if( !trans->child )
      {
        throw runtime_error( "Null child nuclide for " + nuc->symbol );
      }
      
      switch( trans->mode )
      {
        case SandiaDecay::AlphaDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+2) )
            throw runtime_error( "Alpha decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+4) )
            throw runtime_error( "Alpha decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaDecay:
          if( (nuc->atomicNumber+1) != trans->child->atomicNumber )
            throw runtime_error( "Beta decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != trans->child->massNumber )
            throw runtime_error( "Beta decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaPlusDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+1) )
            throw runtime_error( "Beta+ decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != trans->child->massNumber )
            throw runtime_error( "Beta+ decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::ProtonDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+1) )
            throw runtime_error( "Proton decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+1) )
            throw runtime_error( "Proton decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::IsometricTransitionDecay:
          if( nuc->atomicNumber != trans->child->atomicNumber )
            throw runtime_error( "Isometric decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != trans->child->massNumber )
            throw runtime_error( "Isometric decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaAndNeutronDecay:
          if( (nuc->atomicNumber+1) != trans->child->atomicNumber )
            throw runtime_error( "Beta & neutron decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+1) )
            throw runtime_error( "Beta & neutron decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaAndTwoNeutronDecay:
          if( (nuc->atomicNumber+1) != trans->child->atomicNumber )
            throw runtime_error( "Beta & two neutron decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+2) )
            throw runtime_error( "Beta & two neutron decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::ElectronCaptureDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+1) )
            throw runtime_error( "Electron capture decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != trans->child->massNumber )
            throw runtime_error( "Electron capture decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::ElectronCaptureAndProtonDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+2) )
            throw runtime_error( "Electron capture and proton decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+1) )
            throw runtime_error( "Electron capture and proton capture decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::ElectronCaptureAndAlphaDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+3) )
            throw runtime_error( "e- capture and alpha decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+4) )
            throw runtime_error( "e- capture and alpha decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::ElectronCaptureAndTwoProtonDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+3) )
            throw runtime_error( "e- capture and double proton decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+2) )
            throw runtime_error( "e- capture and double proton decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaAndAlphaDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+1) )
            throw runtime_error( "Beta and alpha decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+4) )
            throw runtime_error( "Beta and alpha decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaPlusAndProtonDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+2) )
            throw runtime_error( "Beta+ and proton decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+1) )
            throw runtime_error( "Beta+ and proton decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaPlusAndTwoProtonDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+3) )
            throw runtime_error( "Beta+ and 2 proton decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+2) )
            throw runtime_error( "Beta+ and 2 proton decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaPlusAndThreeProtonDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+4) )
            throw runtime_error( "Beta+ and 3 proton decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+3) )
            throw runtime_error( "Beta+ and 3 proton decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::BetaPlusAndAlphaDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+3) )
            throw runtime_error( "Beta+ and alpha decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+4) )
            throw runtime_error( "Beta+ and alpha decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::DoubleBetaDecay:
          if( (nuc->atomicNumber+2) != trans->child->atomicNumber )
            throw runtime_error( "Double Beta decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != trans->child->massNumber )
            throw runtime_error( "Double Beta decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::DoubleElectronCaptureDecay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+2) )
            throw runtime_error( "Double electron capture decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != trans->child->massNumber )
            throw runtime_error( "Double electron capture decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::Carbon14Decay:
          if( nuc->atomicNumber != (trans->child->atomicNumber+6) )
            throw runtime_error( "C14 decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+14) )
            throw runtime_error( "C14 decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::SpontaneousFissionDecay:
          if( trans->child )
            throw runtime_error( "Spontaneous Fission for " + nuc->symbol + " has a child specified" );
          break;
          
        case SandiaDecay::ClusterDecay:
          //doesnt appear to be used
          break;
          
        case SandiaDecay::DoubleProton:
          if( nuc->atomicNumber != (trans->child->atomicNumber+2) )
            throw runtime_error( "Double proton decay of " + nuc->symbol + " violates atomic number" );
          if( nuc->massNumber != (trans->child->massNumber+2) )
            throw runtime_error( "Double proton decay of " + nuc->symbol + " violates mass number" );
          break;
          
        case SandiaDecay::UndefinedDecay:
          throw runtime_error( "Undefined decay mode for " + nuc->symbol );
          
        default:
          throw runtime_error( "Unrecognized decay mode" );
      }//switch( trans->mode )
      
      
      //Check to make sure coincidences are reasonable (this isnt really a
      //  check on the XML file, but more the logic to associate coincidences)
      for( size_t k = 0; k < trans->products.size(); ++k )
      {
        const RadParticle &part = trans->products[k];
        for( size_t l = 0; l < part.coincidences.size(); ++l )
        {
          if( part.coincidences[l].first >= trans->products.size() )
            throw runtime_error( "A coincidence index for " + nuc->symbol + " is too large" );
          if( part.coincidences[l].second < 0.0f /* || part.concidences[l].second > 1.0f */ )
          {
            throw runtime_error( "A coincidence intensity for " + nuc->symbol
                                 + " is out of range (" + to_str(part.coincidences[l].second) + ")");
          }
        }
      }//for( size_t k = 0; k < trans->products.size(); ++k )
    }//for( size_t j = 0; j < transitions.size(); ++j )
  }//for( size_t i = 0; i < nucs.size(); ++i )
}//void sanity_check_nuclides_and_transistions()


void sanity_check_all_decays()
{
  //This function just makes sure that the decays for all nuclides can be
  //  evaluated; does not check to makesure answers are necassarily correct.
  try
  {
    using namespace SandiaDecay;
    SandiaDecayDataBase database( g_xml_file );
    
    const vector<const Nuclide *> &nuclidesVec = database.nuclides();
    
    size_t num_nuclieds_decayed = 0;
    vector<const Nuclide *>::const_iterator iter;
    for( iter = nuclidesVec.begin(); iter != nuclidesVec.end(); ++iter )
    {
      const Nuclide *nuclide = *iter;
      const double halfLife = nuclide->halfLife;
      
      if( nuclide->isStable() )
        continue;
      
      ++num_nuclieds_decayed;
      NuclideMixture mixture;
      mixture.addNuclideByActivity( nuclide, 1.0*SandiaDecay::curie );
      for( size_t i = 0; i < 20; ++i )
      {
        const double time = i*7.0*halfLife/20.0;
        
        const vector<NuclideActivityPair> activities = mixture.activity( time );
        const vector<EnergyRatePair> gammas = mixture.gammas( time, SandiaDecay::NuclideMixture::OrderByEnergy, true );
      }//for( test the decay for 7 halfLives )
    }//for( iterate over nuclides available, iter )
    
    cout << "All " << num_nuclieds_decayed << " radioactive nuclides have been "
    << " decayed and evaluated at 20 different times, for both daughter "
    << " products, and gammas produced" << endl;
  }catch( std::exception &e )
  {
    throw std::runtime_error( "Problem teting decays: " + string(e.what()) );
  }//try / catch
}//int sanity_check_all_decays();





void check_sum_branching_ratios()
{
  //This is a scratchwork function used to look for issues in sandia.decay.xml
  //  or the SandiaDecay code.
  using namespace SandiaDecay;
  
  SandiaDecayDataBase database( g_xml_file );
  
  const vector<const Nuclide *> &nuclides = database.nuclides();
  
  for( size_t i = 0; i < nuclides.size(); ++i )
  {
    const Nuclide * const nuclide = nuclides[i];
    
    if( nuclide->isStable() )
      continue;
    
    double totalBr = 0.0;
    for( size_t i = 0; i < nuclide->decaysToChildren.size(); ++i )
      totalBr += nuclide->decaysToChildren[i]->branchRatio;
    
    if( (fabs(1.0-totalBr)>0.001) && (totalBr!=0.0) /* && nuclide->halfLife>300.0*/ )
      cout << nuclide->symbol << " (hl=" << nuclide->halfLife << "s) has a BR=" << totalBr << endl;;
  }//for( iterate over nuclides available, iter )
}//int check_sum_branching_ratios();



