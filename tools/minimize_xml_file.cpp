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
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <iostream>
#include <iterator>
#include <algorithm>

#include "SandiaDecay.h"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"


using namespace std;

/** Program to remove the edit history of sandia.decay.xml.
 
 Some items TODO:
 - Floating point values could be limited to reasonable accuracies. E.g., only coincidence fractions
   are currently limited to reasonable accuracy.
 - Decay particles with zero intensity could be removed.
 */

//forward declarations
void remove_id_attrib( rapidxml::xml_node<char> *node );
void remove_cascade_detail_attribs( rapidxml::xml_node<char> *node );
void check_revision_attrib( rapidxml::xml_attribute<char> *att );
void shrink_xml_elements( rapidxml::xml_node<char> *doc_node, bool cascade_detail_present );
void test_db_still_same( const string orig_xml_filename,
                         const string final_xml_filename,
                         const bool coincidences_removed,
                         const int coinc_nsig_fig );
template <typename pod>
std::string pod_to_str( const pod &val );
string print_compact(const double val, const size_t precision);
bool make_float_val_compact(rapidxml::xml_attribute<char>* attrib, const size_t precision);


//Usuage:
//XML_FOREACH_DAUGHTER( child_node_variable, parent_node, "ChildElementName" ){
//  assert( child_node_variable->name() == string("ChildElementName") );
// }
#define XML_FOREACH_DAUGHTER( nodename, parentnode, daughternamestr ) \
for( rapidxml::xml_node<char> *nodename = parentnode->first_node(daughternamestr); \
  nodename; \
  nodename = nodename->next_sibling(daughternamestr) )



int main( int argc, char **argv )
{
  vector<string> args;
  
  /* Removing ID attributes shrinks file from 54MB to 32MB.
     Removing concidences ad ID attribs shrinks to 9.6MB
     Further shrinking tag names gets you down to 6MB
   */
  bool do_remove_id_attrib = false; //setting to true reduces output size from 54MB to 32MB
  bool do_remove_coincidence_gamma = false; //setting to true reduces file size from 32MB to 9.6M
  bool do_shrink_coincidence = false; //reduces ID attribute value to minimum size needed, removes unecassary info from XML.
  bool do_shrink_tag_names = false; //Shrinks XML element and attribute names to one or two character names, ex "nuclide"->"n", "intensity"->"i", etc.  Also removes trailing zeroes from float values.
  const int num_coinc_shrunk_sig_figs = 4;
  bool do_limit_coinc_precision = false; //Shrings coincidence intensity value to have `num_coinc_shrunk_sig_figs` values

  for( int i = 1; i < argc; ++i )
  {
    const string arg = argv[i];
    if( arg == "--remove-id-attrib" )
      do_remove_id_attrib = true;
    else if( arg == "--remove-coinc" )
      do_remove_coincidence_gamma = true;
    else if( arg == "--shrink-coinc" )
      do_shrink_coincidence = true;
    else if (arg == "--limit-coinc-precision")
      do_limit_coinc_precision = true;
    else if( arg == "--shrink-tag-names" )
      do_shrink_tag_names = true;
    else
      args.push_back( arg );
  }//for( loop over the
  
  
  if( args.size() != 2 )
  {
    fprintf( stderr, "Usage: %s <optional flags> <input sandia.decay.xml> <output file>\n"
            "\tWhere optional flags are '--remove-id-attrib', '--remove-coinc', '--shrink-coinc', '--limit-coinc-precision' and/or '--shrink-tag-names'.\n"
            "\t e.x.: ./%s --remove-id-attrib input.sandia.decay.xml output.sandia.decay.xml\n",
            argv[0], argv[0] );
    return EXIT_FAILURE;
  }//if( argc != 3 )
  
  if( do_remove_id_attrib && !do_remove_coincidence_gamma )
  {
    fprintf( stderr, "You cannot remove ID attribute without removing conicidence partivles.\n" );
    return EXIT_FAILURE;
  }
  
  try
  {
    rapidxml::file<char> xmlfile( args[0].c_str() );
    
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_full>( xmlfile.data() );
    rapidxml::xml_node<char> *doc_node = doc.first_node();
  
    if( doc_node->first_attribute("version") )
      doc_node = doc_node->next_sibling();
    
    int ntrans = 0, nodes_removed = 0, num_missing_uuid = 0;
    
    vector<rapidxml::xml_node<char> *> nuc_to_delete;
    for( rapidxml::xml_node<char> *nuc = doc_node->first_node( "nuclide" );
        nuc;
        nuc = nuc->next_sibling("nuclide") )
    {
      rapidxml::xml_attribute<char> *revision = nuc->first_attribute( "revision" );
      check_revision_attrib( revision );
      
      if( revision && rapidxml::internal::compare(revision->value(), revision->value_size(), "deleted", 7, false ) )
      {
        nuc_to_delete.push_back( nuc );
        continue;
      }
      
      if( revision )
        nuc->remove_attribute( revision );
      
      vector<rapidxml::xml_node<char> *> to_delete;
      for( rapidxml::xml_node<char> *child = nuc->first_node(); child; child = child->next_sibling() )
      {
        to_delete.push_back( child );
      }
      for( size_t i = 0; i < to_delete.size(); ++i )
        nuc->remove_node( to_delete[i] );
    }//for( loop over nuclides )
    
    for( size_t i = 0; i < nuc_to_delete.size(); ++i )
      doc_node->remove_node( nuc_to_delete[i] );
    
    for( rapidxml::xml_node<char> *xmltrans = doc_node->first_node( "transition" );
        xmltrans;
        xmltrans = xmltrans->next_sibling("transition") )
    {
      rapidxml::xml_attribute<char> *child_label_attrib = xmltrans->first_attribute("child");
      rapidxml::xml_attribute<char> *parent_label_attrib = xmltrans->first_attribute("parent");
      
      const string child_nuc = child_label_attrib ? child_label_attrib->value() : "";
      const string parent_nuc = parent_label_attrib ? parent_label_attrib->value() : "";
      
      vector<bool> gamma_uuid_used;
      vector<string> gamma_uuid;
      
      // This is a bit of a micro-optimization, but we will set the id element value equal to
      //  the child decay particle index, and then in SandiaDecay::Transition::set we will
      //  optimistically check if the id value is the correct index, and if so, avoid looping
      //  over all decay particles to match id values.
      //  Saves 20 ms initializing DB for a minimized XML with this (micro-)optimization.
      //  (this is about 20% total time when using sandia.decay.min.xml)
      size_t child_index = 0;
      vector<size_t> uuid_particle_index;
      
      vector<rapidxml::xml_node<char> *> gamma_uuid_nodes;
      for( rapidxml::xml_node<char> *child = xmltrans->first_node();
          child;
          child = child->next_sibling() )
      {
        if( rapidxml::internal::compare( child->name(), child->name_size(),
                                         "levelFeed", 9, true ) )
          continue; // metadata, not a particle: never shifts coincidence indices
        if( rapidxml::internal::compare( child->name(), child->name_size(), "gamma", 5, true ) )
        {
          rapidxml::xml_attribute<char> *attrib = child->first_attribute("id");
          if( attrib && attrib->value_size() )
          {
            gamma_uuid.push_back( attrib->value() );
            gamma_uuid_used.push_back( false );
            gamma_uuid_nodes.push_back( child );
            uuid_particle_index.push_back( child_index );
          }
        }//if( <gamma> node )
        
        child_index += 1;
      }//for( loop over xmltrans's children )
      
      
      vector<rapidxml::xml_node<char> *> to_delete;
      for( rapidxml::xml_node<char> *child = xmltrans->first_node();
          child;
          child = child->next_sibling() )
      {
        if( rapidxml::internal::compare( child->name(), child->name_size(),
                                         "levelFeed", 9, true ) )
        {
          // Direct alpha/beta-minus source physics, not coincidence-only data.
          // It must survive --remove-coinc.
          continue;
        }
        if( do_remove_id_attrib )
          remove_id_attrib( child );
        
        rapidxml::xml_attribute<char> *id_sttriv = child->first_attribute( "id" );
        rapidxml::xml_attribute<char> *revision = child->first_attribute( "revision" );
        
        check_revision_attrib( revision );
        
        const string child_uuid = id_sttriv ? id_sttriv->value() : "";
        
        if( rapidxml::internal::compare(child->name(), child->name_size(), "insertion", 9, false )
           || rapidxml::internal::compare(child->name(), child->name_size(), "edit", 4, false ) )
        {
          to_delete.push_back( child );
        }else if( revision && rapidxml::internal::compare(revision->value(), revision->value_size(), "deleted", 7, false ) )
        {
          to_delete.push_back( child );
        }else if( child->type() == rapidxml::node_comment )
        {
          to_delete.push_back( child );
        }else
        {
          if( rapidxml::internal::compare(child->name(),child->name_size(),"beta",4,false)
             || rapidxml::internal::compare(child->name(),child->name_size(),"positron",8,false)
             || rapidxml::internal::compare(child->name(),child->name_size(),"electronCapture",15,false)
             || rapidxml::internal::compare(child->name(),child->name_size(),"alpha",5,false) )
          {
            rapidxml::xml_attribute<char> *attrib = child->first_attribute("id");
            if( attrib )
              child->remove_attribute( attrib );
          }

          // With coincidences removed there is no cascade summing, so the
          // per-gamma / per-electronCapture cascade-detail attributes are unused;
          // strip them so the nocoinc file stays lean.
          if( do_remove_coincidence_gamma )
            remove_cascade_detail_attribs( child );

          vector<rapidxml::xml_node<char> *> gchild_to_delete;
          for( rapidxml::xml_node<char> *gchild = child->first_node();
              gchild;
              gchild = gchild->next_sibling() )
          {
            if( do_remove_id_attrib )
              remove_id_attrib( gchild );
            
            if( rapidxml::internal::compare(gchild->name(), gchild->name_size(), "edit", 4, false )
               || rapidxml::internal::compare(gchild->name(), gchild->name_size(), "insertion", 9, false )
               || rapidxml::internal::compare(gchild->name(), gchild->name_size(), "deletion", 8, false )
               || (do_remove_coincidence_gamma && rapidxml::internal::compare(gchild->name(),gchild->name_size(),"coincidentGamma",15,false) )
               || gchild->type() == rapidxml::node_comment )
            {
              gchild_to_delete.push_back( gchild );
            }else if( gchild->name_size() < 1 )
            {
              printf( "Found un-named node\n" );
            }else if( rapidxml::internal::compare(gchild->name(), gchild->name_size(), "coincidentgamma", 15, false ) )
            {
              if( do_shrink_coincidence )
              {
                //gchild->name("coinc");
                rapidxml::xml_attribute<char> *attrib = gchild->first_attribute("energy");
                if( attrib )
                  gchild->remove_attribute( attrib );
                
                attrib = gchild->first_attribute("id");
                if( attrib )
                {
                  //child_uuid
                  const string gchild_uuid = attrib->value();
                  assert( gchild_uuid != child_uuid );
                  vector<string>::const_iterator pos = std::find( gamma_uuid.begin(), gamma_uuid.end(), gchild_uuid );
                  //assert( pos != gamma_uuid.end() );
                  if( pos == gamma_uuid.end() )
                  {
                    num_missing_uuid += 1;
                    gchild->remove_attribute( attrib );
                    
                    fprintf( stderr, "Warning: Could not find UUID='%s' for %s->%s - removing coincidence.\n",
                            gchild_uuid.c_str(), parent_nuc.c_str(), child_nuc.c_str() );
                  }else
                  {
                    const int64_t index = pos - gamma_uuid.begin();
                    assert( index < uuid_particle_index.size() );
                    const size_t child_index = uuid_particle_index[index];
                    gamma_uuid_used[index] = true;
                    const string newuuid = pod_to_str(child_index);
                    const char *newuuidstr = doc.allocate_string( newuuid.c_str() );
                    attrib->value( newuuidstr );
                  }
                }
              }//if( do_shrink_coincidence )

              if (do_limit_coinc_precision)
              {
                rapidxml::xml_attribute<char>* attrib = gchild->first_attribute("intensity");
                if (attrib)
                {
                  if( !make_float_val_compact(attrib, num_coinc_shrunk_sig_figs) )
                    printf("Error shrinking coincidence precision.\n");
                }
                else
                {
                  printf("Coincidence node with no intensity value.\n");
                }
              }//if (do_limit_coinc_precision)

            }else //if( !rapidxml::internal::compare(gchild->name(), gchild->name_size(), "coincidentgamma", 15, false ) )
            {
              throw runtime_error( "Unknown gchild node: '"
                                   + std::string(gchild->name(), gchild->name() + gchild->name_size()) + "'" );
            }
          }
          for( size_t i = 0; i < gchild_to_delete.size(); ++i )
            child->remove_node( gchild_to_delete[i] );
        
          if( revision )
            child->remove_attribute( revision );
        }
      }//for( loop over children )
    
      rapidxml::xml_attribute<char> *attrib = xmltrans->first_attribute("remarks");
      if( attrib )
      {
        string val = attrib->value();
        if( (val.find("warning:") == string::npos) && (val.find("audit:") == string::npos) )
          xmltrans->remove_attribute( attrib );
      }
      
      
      if( do_shrink_coincidence )
      {
        for( size_t i = 0; i < gamma_uuid.size(); ++i )
        {
          rapidxml::xml_attribute<char> *attrib = gamma_uuid_nodes[i]->first_attribute("id");
          
          if( !gamma_uuid_used[i] )
          {
            if( attrib )
              gamma_uuid_nodes[i]->remove_attribute(attrib);
          }else
          {
            assert( attrib );
            const size_t child_index = uuid_particle_index[i];
            const char *newuuidstr = doc.allocate_string( pod_to_str(child_index).c_str() );
            attrib->value( newuuidstr );
          }
        }
      }//if( do_shrink_coincidence )
      
      
      nodes_removed += to_delete.size();
      for( size_t i = 0; i < to_delete.size(); ++i )
        xmltrans->remove_node( to_delete[i] );
      
      ++ntrans;
    }//loop over transitions.
    
  
    for( rapidxml::xml_node<char> *node = doc_node->first_node(); node;  )
    {
      rapidxml::xml_node<char> *next_node = node->next_sibling();
      if( node->type() == rapidxml::node_comment )
      {
        doc_node->remove_node( node );
        ++nodes_removed;
      }
      
      node = next_node;
    }
    
    if( do_shrink_tag_names )
      shrink_xml_elements( doc_node, !do_remove_coincidence_gamma );
      
    printf( "There were %i transitions; removed %i xml elements, and %i missing coincidences.\n",
           ntrans, nodes_removed, num_missing_uuid );
    
    {
      ifstream test( args[1].c_str() );
      if( test.is_open() )
        throw runtime_error( "File '" + args[1] + " already exists - not overwriting!" );
    }
    
    {
      ofstream outputfile( args[1].c_str(), ios::out | ios::binary );
      if( !outputfile )
        throw runtime_error( "Could not open output file: " + string(args[1]) );
    
      rapidxml::print<char>( outputfile, doc, 0 );
    }
    
    printf( "Saved output file '%s'\n", args[1].c_str() );
    
    const int coinc_nsig_fig = do_limit_coinc_precision ? num_coinc_shrunk_sig_figs : 0;
    
    test_db_still_same( args[0], args[1], do_remove_coincidence_gamma, coinc_nsig_fig );
  }catch( std::exception &e )
  {
    fprintf( stderr, "Error: %s\n", e.what() );
    // Validation and output errors are delivery failures. Falling through to
    // EXIT_SUCCESS lets automation accept a bad generated database.
    return EXIT_FAILURE;
  }//try / catch
  
  
  
  return EXIT_SUCCESS;
}//int main( int argc, char **argv )



void test_db_still_same( const string orig_xml_filename, const string final_xml_filename,
                        const bool coincidences_removed,
                        const int coinc_nsig_fig )
{
  using namespace SandiaDecay;
  
  SandiaDecayDataBase olddb( orig_xml_filename );
  SandiaDecayDataBase newdb( final_xml_filename );
  
  const vector<const Nuclide *> &oldnucs = olddb.nuclides();
  const vector<const Nuclide *> &newnucs = newdb.nuclides();
  
  if( oldnucs.size() != newnucs.size() )
    throw runtime_error( "Number of nuclides changed" );
  
  for( size_t i = 0; i < oldnucs.size(); ++i )
  {
    const Nuclide * const oldnuc = oldnucs[i];
    const Nuclide * const newnuc = newnucs[i];
    
    if( oldnuc->symbol != newnuc->symbol )
      throw runtime_error( "old and new symbols dont match" );
    
    if( oldnuc->atomicNumber != newnuc->atomicNumber )
      throw runtime_error( "old and new  dont match" );
    
    if( oldnuc->massNumber != newnuc->massNumber )
      throw runtime_error( "old and new  dont match" );
    
    if( oldnuc->isomerNumber != newnuc->isomerNumber )
      throw runtime_error( "old and new  dont match" );
    
    if( oldnuc->atomicMass != newnuc->atomicMass )
      throw runtime_error( "old and new  dont match" );
    
    if( oldnuc->halfLife != newnuc->halfLife )
      throw runtime_error( "old and new  dont match" );
    
    if( oldnuc->decaysFromParents.size() != newnuc->decaysFromParents.size() )
      throw runtime_error( "old and new number of transitions from parents dont match for " + oldnuc->symbol
                          + "(" + pod_to_str(oldnuc->decaysFromParents.size()) + " vs "
                           + pod_to_str(newnuc->decaysFromParents.size()) + ")" );
    
    if( oldnuc->decaysToChildren.size() != newnuc->decaysToChildren.size() )
      throw runtime_error( "old and new number of transitions to children dont match for " + oldnuc->symbol
                          + "(" + pod_to_str(oldnuc->decaysFromParents.size()) + " vs "
                          + pod_to_str(newnuc->decaysFromParents.size()) + ")" );
  
    std::vector<const Transition *> oldtrans = oldnuc->decaysFromParents;
    oldtrans.insert( oldtrans.end(), oldnuc->decaysToChildren.begin(), oldnuc->decaysToChildren.end() );
    
    std::vector<const Transition *> newtrans = newnuc->decaysFromParents;
    newtrans.insert( newtrans.end(), newnuc->decaysToChildren.begin(), newnuc->decaysToChildren.end() );
    
    assert( oldtrans.size() == newtrans.size() );
    
    for( size_t trans_index = 0; trans_index < oldtrans.size(); ++trans_index )
    {
      const Transition * const oldtran = oldtrans[trans_index];
      const Transition * const newtran = newtrans[trans_index];
      
      if( (!oldtran->parent) != (!newtran->parent) )
        throw runtime_error( "Availability of parent on old/new transition not equal" );
      
      if( (!oldtran->child) != (!newtran->child) )
        throw runtime_error( "Availability of parent on old/new transition not equal" );
      
      if( oldtran->parent && (oldtran->parent->symbol != newtran->parent->symbol) )
        throw runtime_error( "Transition parent not equal" );
      
      if( oldtran->child && (oldtran->child->symbol != newtran->child->symbol) )
        throw runtime_error( "Transition child not equal" );
      
      if( oldtran->mode != newtran->mode )
        throw runtime_error( "Transition decay mode not equal" );
      
      if( oldtran->branchRatio != newtran->branchRatio )
        throw runtime_error( "Transition branchRatio not equal" );

      if( oldtran->level_feeds.size() != newtran->level_feeds.size() )
        throw runtime_error( "Transition level-feed count changed" );
      for( size_t feed_index = 0; feed_index < oldtran->level_feeds.size(); ++feed_index )
      {
        const LevelFeed &oldfeed = oldtran->level_feeds[feed_index];
        const LevelFeed &newfeed = newtran->level_feeds[feed_index];
        const double pscale = std::max( std::fabs(static_cast<double>(oldfeed.probability)),
                                        std::fabs(static_cast<double>(newfeed.probability)) );
        const double qscale = std::max( std::fabs(static_cast<double>(oldfeed.q_keV)),
                                        std::fabs(static_cast<double>(newfeed.q_keV)) );
        if( oldfeed.level != newfeed.level ||
            oldfeed.forbiddenness != newfeed.forbiddenness ||
            (oldfeed.probability != newfeed.probability &&
             (pscale == 0.0 ||
              std::fabs(static_cast<double>(oldfeed.probability) -
                        static_cast<double>(newfeed.probability)) > 1.1e-6 * pscale)) ||
            (oldfeed.q_keV != newfeed.q_keV &&
             (qscale == 0.0 ||
              std::fabs(static_cast<double>(oldfeed.q_keV) -
                        static_cast<double>(newfeed.q_keV)) > 1.1e-6 * qscale)) )
          throw runtime_error( "Transition level-feed data changed" );
      }
      
      if( oldtran->products.size() != newtran->products.size() )
        throw runtime_error( "Transition number of products not the same size" );
      
      for( size_t prod_index = 0; prod_index < oldtran->products.size(); ++prod_index )
      {
        const RadParticle &oldpart = oldtran->products[prod_index];
        const RadParticle &newpart = newtran->products[prod_index];
        
        if( oldpart.type != newpart.type )
          throw runtime_error( "Particle type does not match" );
        
        if( oldpart.energy != newpart.energy )
          throw runtime_error( "Particle energy does not match" );
        
        if( oldpart.intensity != newpart.intensity )
          throw runtime_error( "Particle intensity does not match" );
        
        if( oldpart.hindrance != newpart.hindrance )
          throw runtime_error( "Particle hindrance does not match" );
        
        if( oldpart.logFT != newpart.logFT )
          throw runtime_error( "Particle logFT does not match" );
        
        if( oldpart.forbiddenness != newpart.forbiddenness )
          throw runtime_error( "Particle forbiddenness does not match" );

        // Verified pure-E0 provenance controls whether a photon physically
        // exists; it must survive every minimized form, including --remove-coinc.
        if( oldpart.e0_verified != newpart.e0_verified )
          throw runtime_error( "Particle verified-E0 provenance does not match" );

        // ICC and EC coefficients are compacted to nominally four significant
        // digits below. The legacy float-based compactor can truncate one
        // additional unit near a rounding edge, so validate their round trip
        // with the same bounded tolerance in both minimized output modes.
        const auto compact_equal = []( const float a, const float b ) {
          if( a == b ) return true;
          const double scale = std::max( std::fabs(static_cast<double>(a)),
                                         std::fabs(static_cast<double>(b)) );
          return scale > 0.0 &&
                 std::fabs(static_cast<double>(a) - static_cast<double>(b))
                     <= 1.1e-3 * scale;
        };

        // Cascade-detail fields are part of the parsed database contract too.
        // In a coincidence-bearing minimized file they must survive long->short
        // renaming (notably positron/electronCapture lvl->lv). In a --remove-coinc
        // file its topology must be absent/defaulted, while physical gamma ICC
        // and EC shell fractions must remain available to fallback consumers.
        if( !coincidences_removed )
        {
          if( oldpart.from_level != newpart.from_level
             || oldpart.to_level != newpart.to_level
             || oldpart.fed_level != newpart.fed_level
             || oldpart.multipole != newpart.multipole
             || oldpart.mixing_ratio != newpart.mixing_ratio
             || !compact_equal(oldpart.icc_total, newpart.icc_total)
             || !compact_equal(oldpart.icc_K, newpart.icc_K)
             || !compact_equal(oldpart.icc_L1, newpart.icc_L1)
             || !compact_equal(oldpart.icc_L2, newpart.icc_L2)
             || !compact_equal(oldpart.icc_L3, newpart.icc_L3)
             || !compact_equal(oldpart.ec_K, newpart.ec_K)
             || !compact_equal(oldpart.ec_L, newpart.ec_L)
             || !compact_equal(oldpart.ec_M, newpart.ec_M)
             || !compact_equal(oldpart.ec_N, newpart.ec_N) )
            throw runtime_error( "Particle cascade detail does not match" );
        }
        else if( newpart.from_level != -1 || newpart.to_level != -1
                 || newpart.fed_level != -1 || newpart.mixing_ratio != 0.0f )
        {
          throw runtime_error( "Cascade detail survived --remove-coinc" );
        }

        if( coincidences_removed )
        {
          if( oldpart.type == GammaParticle )
          {
            if( oldpart.multipole != newpart.multipole
                || !compact_equal(oldpart.icc_total, newpart.icc_total)
                || !compact_equal(oldpart.icc_K, newpart.icc_K)
                || !compact_equal(oldpart.icc_L1, newpart.icc_L1)
                || !compact_equal(oldpart.icc_L2, newpart.icc_L2)
                || !compact_equal(oldpart.icc_L3, newpart.icc_L3) )
            {
              std::ostringstream msg;
              msg << "Gamma physical-emission metadata changed in --remove-coinc for "
                  << (oldtran->parent ? oldtran->parent->symbol : string("?"))
                  << "->" << (oldtran->child ? oldtran->child->symbol : string("?"))
                  << " at " << oldpart.energy << " keV: mult "
                  << oldpart.multipole << "->" << newpart.multipole
                  << ", ICC " << oldpart.icc_total << "->" << newpart.icc_total
                  << ", K " << oldpart.icc_K << "->" << newpart.icc_K
                  << ", L1 " << oldpart.icc_L1 << "->" << newpart.icc_L1
                  << ", L2 " << oldpart.icc_L2 << "->" << newpart.icc_L2
                  << ", L3 " << oldpart.icc_L3 << "->" << newpart.icc_L3;
              throw runtime_error( msg.str() );
            }
          }
          else if( oldpart.type == CaptureElectronParticle )
          {
            if( !compact_equal(oldpart.ec_K, newpart.ec_K)
                || !compact_equal(oldpart.ec_L, newpart.ec_L)
                || !compact_equal(oldpart.ec_M, newpart.ec_M)
                || !compact_equal(oldpart.ec_N, newpart.ec_N) )
              throw runtime_error(
                  "Electron-capture shell fractions changed in --remove-coinc" );
          }
          else if( newpart.multipole != 0 || newpart.icc_total != 0.0f
                   || newpart.icc_K != 0.0f || newpart.icc_L1 != 0.0f
                   || newpart.icc_L2 != 0.0f || newpart.icc_L3 != 0.0f
                   || newpart.ec_K != 0.0f || newpart.ec_L != 0.0f
                   || newpart.ec_M != 0.0f || newpart.ec_N != 0.0f )
          {
            throw runtime_error(
                "Gamma physical-emission metadata leaked onto another particle" );
          }
        }
        
        if( !coincidences_removed )
        {
          const auto &oldcoinc = oldpart.coincidences;
          const auto &newcoinc = newpart.coincidences;
          
          if( oldcoinc.size() != newcoinc.size() )
          {
            const string parent = oldtran->parent ? oldtran->parent->symbol : string();
            const string child = oldtran->child ? oldtran->child->symbol : string();
            const string nold = to_string(oldcoinc.size());
            const string nnew = to_string(newcoinc.size());
            throw runtime_error( "Number of coincidences does not match for "
                                + parent + " -> " + child + ": from " + nold + " to " + nnew
                                + " for gamma " + std::to_string(oldpart.energy) + " keV" );
          }//if( oldcoinc.size() != newcoinc.size() )
          
          //num_coinc_shrunk_sig_figs
          
          for( size_t coinc_index = 0; coinc_index < oldcoinc.size(); ++coinc_index )
          {
            if( oldcoinc[coinc_index].first != newcoinc[coinc_index].first )
              throw runtime_error( "Coincidence index does not match" );
            
            if( coinc_nsig_fig > 0 )
            {
              const double diff = fabs( oldcoinc[coinc_index].second - newcoinc[coinc_index].second);
              
              if( diff > (oldcoinc[coinc_index].second * std::pow(10.0,-1.0*(coinc_nsig_fig-1))) )
                throw runtime_error( "Coincidence rate does not match (old="
                                    + std::to_string(oldcoinc[coinc_index].second)
                                    + ", new=" + std::to_string(newcoinc[coinc_index].second)
                                    + ", diff=" + std::to_string(diff)
                                    + ")" );
            }else
            {
              //We didnt limit precision
              if( oldcoinc[coinc_index].second != newcoinc[coinc_index].second )
                throw runtime_error( "Coincidence rate does not match" );
            }
          }//for( loop over coincidences )
        }//if( coincidenses shouldnt hav echanged )
      }//for( loop over decay particles )
    }//for( loop over transisitons )
  }//for( loop over nuclides )
  
  
  const vector<const Element *> &oldels = olddb.elements();
  const vector<const Element *> &newels = newdb.elements();
  if( oldels.size() != newels.size() )
    throw runtime_error( "Number of elements does not match" );
  
  for( size_t el_index = 0; el_index < oldels.size(); ++el_index )
  {
    const Element *oldel = oldels[el_index];
    const Element *newel = newels[el_index];
    
    if( oldel->symbol != newel->symbol )
      throw runtime_error( "Element symbols do not match" );
    
    if( oldel->name != newel->name )
      throw runtime_error( "Element name do not match" );
    
    if( oldel->atomicNumber != newel->atomicNumber )
      throw runtime_error( "Element atomicNumber do not match" );
    
    if( oldel->isotopes.size() != newel->isotopes.size() )
      throw runtime_error( "Element isotopes size is different" );
    
    if( oldel->xrays.size() != newel->xrays.size() )
      throw runtime_error( "Element xrays size is different" );
    
    for( size_t iso_index = 0; iso_index < oldel->isotopes.size(); ++iso_index )
    {
      const NuclideAbundancePair &oldiso = oldel->isotopes[iso_index];
      const NuclideAbundancePair &newiso = newel->isotopes[iso_index];
      
      if( oldiso.nuclide->symbol != newiso.nuclide->symbol )
        throw runtime_error( "Element isos nuclide do not match" );
      
      if( oldiso.abundance != newiso.abundance )
        throw runtime_error( "Element isos abundance do not match" );
    }//for( loop over isotopes )
    
    for( size_t x_index = 0; x_index < oldel->xrays.size(); ++x_index )
    {
      const EnergyIntensityPair &oldxray = oldel->xrays[x_index];
      const EnergyIntensityPair &newxray = newel->xrays[x_index];
      
      if( oldxray.energy != newxray.energy )
        throw runtime_error( "Element xray energy do not match" );
      
      if( oldxray.intensity != newxray.intensity )
        throw runtime_error( "Element intensity energy do not match" );
    }//for( loop over xrays )
  }//for( loop over elements )
  
  printf( "Old and new XML agree over all validated SandiaDecayDataBase fields\n" );
}//void test_db_still_same()


template <typename pod>
std::string pod_to_str( const pod &val ){
  std::stringstream srm;
  srm << val;
  return srm.str();
}


void rename_attrib( rapidxml::xml_node<char> *node, const char *oldname, const char *newname )
{
  rapidxml::xml_attribute<char> *att = node->first_attribute(oldname);
  if( att )
    att->name( newname );
}

void remove_cascade_detail_attribs( rapidxml::xml_node<char> *node )
{
  // Coincidence-only enrichment (level topology and multipole mixing ratio) is
  // dead weight when coincidences are removed, so strip it to keep the nocoinc
  // file lean. Gamma ICC, multipole, e0Verified, and electron-capture shell
  // fractions are deliberately NOT in this list: together they determine
  // whether a photon physically exists and the physical IC/EC vacancy yield,
  // even for a consumer using the fallback rather than a level graph. Handles
  // long names (this runs before shrink_xml_elements).
  static const char * const kAttribs[] = {
    "fl", "tl", "delta",                                      // <gamma>
    "lvl"                                              // <electronCapture>/<positron>
  };
  for( const char *name : kAttribs )
  {
    rapidxml::xml_attribute<char> *att = node->first_attribute( name );
    if( att )
      node->remove_attribute( att );
  }
}

void remove_trailing_zeroes( rapidxml::xml_attribute<char> *n )
{
  //"1.0000E-03" -> "1E-03"
  //"2.7500E-04" -> "2.75E-04"
  //"1.1000" -> "1.1"
  //"1.0000" -> "1"
  //"1.00E+00" -> "1"
  if( !n )
    return;
  
  const string origval = n->value();
  string smallval = origval;
  
  const size_t ppos = smallval.find( "." );
  const size_t epos = smallval.find_first_of( "Ee" );
  
  if( ppos == string::npos && epos == string::npos )
    return;
  
  string epart;
  if( epos != string::npos )
  {
    assert( epos != 0 );
    assert( (ppos == string::npos) || (ppos < epos) );
    epart = smallval.substr(epos);
    smallval = smallval.substr(0,epos);
    if( epart == "E+00" || epart == "e+00" || epart == "E00" )
      epart = "";
  }//if( epos != string::npos )

  if( ppos != string::npos )
  {
    while( smallval[smallval.size()-1] == '0' )
    {
      smallval = smallval.substr(0,smallval.size()-1);
      assert( smallval.size() );
    }
    
    if( smallval[smallval.size()-1] == '.' )
      smallval = smallval.substr(0,smallval.size()-1);
    if( smallval.empty() ) //original number was like "0.0000"
      smallval = "0";
  }//if( ppos != string::npos )
  
  smallval += epart;
  
  const double orig = strtod( origval.c_str(), NULL );
  const double modded = strtod( smallval.c_str(), NULL );
  if( orig != modded )
  {
    fprintf( stderr, "Removing trailing zeros of '%s' let to %f->%f\n", n->value(), orig, modded );
  }
  assert( orig == modded );
  
  n->value( n->document()->allocate_string(smallval.c_str()) );
}//void remove_trailing_zeroes(...)


void shrink_xml_elements( rapidxml::xml_node<char> *doc_node, bool cascade_detail_present )
{
  XML_FOREACH_DAUGHTER( nuc, doc_node, "nuclide" )
  {
    nuc->name( "n" );
    
    remove_trailing_zeroes( nuc->first_attribute("atomicMass") );
    remove_trailing_zeroes( nuc->first_attribute("halfLife") );
    
    rename_attrib( nuc, "atomicMass", "am" );
    rename_attrib( nuc, "atomicNumber", "an" );
    rename_attrib( nuc, "halfLife", "hl" );
    rename_attrib( nuc, "isomerNumber", "iso" );
    rename_attrib( nuc, "massNumber", "mn" );
    rename_attrib( nuc, "symbol", "s" );
  }//foreach( nuclide )
  
  
  XML_FOREACH_DAUGHTER( trans, doc_node, "transition" )
  {
    trans->name( "t" );
    
    remove_trailing_zeroes( trans->first_attribute("branchRatio") );
    
    rename_attrib( trans, "branchRatio", "br" );
    rename_attrib( trans, "child", "c" );
    rename_attrib( trans, "mode", "m" );
    rename_attrib( trans, "parent", "p" );

    XML_FOREACH_DAUGHTER( feed, trans, "levelFeed" )
    {
      feed->name( "lf" );
      make_float_val_compact( feed->first_attribute("probability"), 8 );
      make_float_val_compact( feed->first_attribute("q"), 8 );
      rename_attrib( feed, "lvl", "l" );
      rename_attrib( feed, "probability", "p" );
      rename_attrib( feed, "forbiddenness", "f" );
    }
    
    XML_FOREACH_DAUGHTER( child, trans, "beta" )
    {
      child->name( "b" );
      
      remove_trailing_zeroes( child->first_attribute("energy") );
      remove_trailing_zeroes( child->first_attribute("intensity") );
      remove_trailing_zeroes( child->first_attribute("logFT") );
      remove_trailing_zeroes( child->first_attribute("forbiddenness") );
      
      rename_attrib( child, "energy", "e" );
      rename_attrib( child, "intensity", "i" );
      rename_attrib( child, "logFT", "lft" );
      rename_attrib( child, "forbiddenness", "f" );
    }
    
    XML_FOREACH_DAUGHTER( child, trans, "positron" )
    {
      child->name( "p" );
      
      remove_trailing_zeroes( child->first_attribute("energy") );
      remove_trailing_zeroes( child->first_attribute("intensity") );
      remove_trailing_zeroes( child->first_attribute("logFT") );
      remove_trailing_zeroes( child->first_attribute("forbiddenness") );
      
      rename_attrib( child, "energy", "e" );
      rename_attrib( child, "intensity", "i" );
      rename_attrib( child, "logFT", "lft" );
      rename_attrib( child, "forbiddenness", "f" );
      rename_attrib( child, "lvl", "lv" );  // fed level index (integer)
    }
    
    XML_FOREACH_DAUGHTER( child, trans, "electronCapture" )
    {
      child->name( "ec" );

      remove_trailing_zeroes( child->first_attribute("energy") );
      remove_trailing_zeroes( child->first_attribute("intensity") );
      remove_trailing_zeroes( child->first_attribute("logFT") );
      remove_trailing_zeroes( child->first_attribute("forbiddenness") );

      rename_attrib( child, "energy", "e" );
      rename_attrib( child, "intensity", "i" );
      rename_attrib( child, "logFT", "lft" );
      rename_attrib( child, "forbiddenness", "f" );

      // Cascade-detail: electron-capture shell fractions + fed level index.
      make_float_val_compact( child->first_attribute("ecK"), 4 );
      make_float_val_compact( child->first_attribute("ecL"), 4 );
      make_float_val_compact( child->first_attribute("ecM"), 4 );
      make_float_val_compact( child->first_attribute("ecN"), 4 );
      rename_attrib( child, "ecK", "ek" );
      rename_attrib( child, "ecL", "el" );
      rename_attrib( child, "ecM", "em" );
      rename_attrib( child, "ecN", "en" );
      rename_attrib( child, "lvl", "lv" );  // fed level index (integer)
    }
    
    XML_FOREACH_DAUGHTER( child, trans, "alpha" )
    {
      child->name( "a" );
      
      remove_trailing_zeroes( child->first_attribute("energy") );
      remove_trailing_zeroes( child->first_attribute("intensity") );
      remove_trailing_zeroes( child->first_attribute("hindrance") );
      remove_trailing_zeroes( child->first_attribute("forbiddenness") );
      
      rename_attrib( child, "energy", "e" );
      rename_attrib( child, "intensity", "i" );
      rename_attrib( child, "hindrance", "h" );
      rename_attrib( child, "forbiddenness", "f" );
    }
    
    XML_FOREACH_DAUGHTER( child, trans, "xray" )
    {
      child->name( "x" );
      
      remove_trailing_zeroes( child->first_attribute("energy") );
      remove_trailing_zeroes( child->first_attribute("intensity") );
      
      rename_attrib( child, "energy", "e" );
      rename_attrib( child, "intensity", "i" );
      rename_attrib( child, "assignment", "l" );
      
      rename_attrib( child, "l1intensity", "l1i" );
      rename_attrib( child, "l2intensity", "l2i" );
      rename_attrib( child, "l3intensity", "l3i" );
      
    }
    
    XML_FOREACH_DAUGHTER( gamma, trans, "gamma" )
    {
      gamma->name( "g" );
      
      remove_trailing_zeroes( gamma->first_attribute("energy") );
      remove_trailing_zeroes( gamma->first_attribute("intensity") );
      
      rename_attrib( gamma, "energy", "e" );
      rename_attrib( gamma, "intensity", "i" );

      // Cascade-detail: internal-conversion coefficients. 4 significant figures
      // is far below the physical uncertainty on these. Level indices fl/tl are
      // small integers -> left as-is.
      make_float_val_compact( gamma->first_attribute("iccTotal"), 4 );
      make_float_val_compact( gamma->first_attribute("iccK"), 4 );
      make_float_val_compact( gamma->first_attribute("iccL1"), 4 );
      make_float_val_compact( gamma->first_attribute("iccL2"), 4 );
      make_float_val_compact( gamma->first_attribute("iccL3"), 4 );
      rename_attrib( gamma, "iccTotal", "ic" );
      rename_attrib( gamma, "iccK", "ik" );
      rename_attrib( gamma, "iccL1", "il1" );
      rename_attrib( gamma, "iccL2", "il2" );
      rename_attrib( gamma, "iccL3", "il3" );
      rename_attrib( gamma, "e0Verified", "e0" );

      XML_FOREACH_DAUGHTER( coinc, gamma, "coincidentGamma" )
      {
        coinc->name("c");

        remove_trailing_zeroes( coinc->first_attribute("energy") );
        remove_trailing_zeroes( coinc->first_attribute("intensity") );

        rename_attrib( coinc, "energy", "e" );
        rename_attrib( coinc, "intensity", "i" );

        // Cascade-detail: gamma-gamma angular-correlation coefficients (names a2,
        // a4 are already minimal). 4 sig figs is ample.
        make_float_val_compact( coinc->first_attribute("a2"), 4 );
        make_float_val_compact( coinc->first_attribute("a4"), 4 );
      }
    }
  }//foreach( transition )
  
  
  XML_FOREACH_DAUGHTER( el, doc_node, "element" )
  {
    el->name( "el" );
    
    rename_attrib( el, "atomicNumber", "an" );
    rename_attrib( el, "name", "n" );
    rename_attrib( el, "symbol", "s" );
    
    XML_FOREACH_DAUGHTER( x, el, "xray" )
    {
      x->name( "x" );
      remove_trailing_zeroes( x->first_attribute("energy") );
      remove_trailing_zeroes( x->first_attribute("intensity") );
    
      rename_attrib( x, "relintensity", "i" );
      rename_attrib( x, "energy", "e" );
    }
    
    XML_FOREACH_DAUGHTER( x, el, "isotope" )
    {
      x->name( "iso" );
      rename_attrib( x, "abundance", "a" );
      rename_attrib( x, "symbol", "s" );
    }
  }//foreach( element )
  
  // Header comment: name mappings + a plain-language field glossary so a reader
  // who is not a nuclear-data expert can interpret the minimized file. (XML
  // comments may not contain a double hyphen, so only single hyphens are used.)
  std::string desc =
  "\nXML element/attribute names have been renamed to minimize file size, with the following mappings:\n"
  "\t'transition'      to 't'\n"
  "\t'branchRatio'     to 'br'\n"
  "\t'child'           to 'c'\n"
  "\t'mode'            to 'm'\n"
  "\t'parent'          to 'p'\n"
  "\t'productIntensityBasis' kept as-is (reviewed product-intensity provenance)\n"
  "\t'beta'            to 'b'\n"
  "\t'energy'          to 'e'\n"
  "\t'intensity'       to 'i'\n"
  "\t'logFT'           to 'lft'\n"
  "\t'forbiddenness'   to 'f'\n"
  "\t'electronCapture' to 'ec'\n"
  "\t'alpha'           to 'a'\n"
  "\t'gamma'           to 'g'\n"
  "\t'coincidentGamma' to 'c'\n"
  "\t'element'         to 'el'\n"
  "\t'atomicNumber'    to 'an'\n"
  "\t'name'            to 'n'\n"
  "\t'symbol'          to 's'\n"
  "\t'xray'            to 'x'\n"
  "\t'relintensity'    to 'i'\n"
  "\t'abundance'       to 'a'\n";

  if( cascade_detail_present )
  {
    desc +=
    "  Cascade-detail attributes (present only in this coincidence-bearing file).\n"
    "  ICC, topology and angular data are from GEANT4 11.4.0\n"
    "  PhotonEvaporation/RadioactiveDecay; verified-E0 provenance is separately\n"
    "  curated from the primary NNDC ENSDF evaluation:\n"
    "\t'iccTotal'        to 'ic'\n"
    "\t'iccK'            to 'ik'\n"
    "\t'iccL1'           to 'il1'\n"
    "\t'iccL2'           to 'il2'\n"
    "\t'iccL3'           to 'il3'\n"
    "\t'ecK'             to 'ek'\n"
    "\t'ecL'             to 'el'\n"
    "\t'lvl'             to 'lv' (on electronCapture and positron)\n"
    "\t'levelFeed'       to 'lf'; lvl/probability/forbiddenness to l/p/f\n"
    "\t'e0Verified'      to 'e0' (reviewed pure-E0 provenance on gamma)\n"
    "\t('fl','tl','a2','a4' are kept as-is)\n";
  }
  else
  {
    // Compact physical gamma-emission metadata remains in the
    // coincidence-free file even though topology/correlation detail is gone.
    desc +=
    "  Physical gamma-emission attributes retained in this file:\n"
    "\t'iccTotal'        to 'ic'\n"
    "\t'iccK'            to 'ik'\n"
    "\t'iccL1'           to 'il1'\n"
    "\t'iccL2'           to 'il2'\n"
    "\t'iccL3'           to 'il3'\n"
    "\t'e0Verified'      to 'e0' (reviewed pure-E0 provenance on gamma)\n"
    "\t'mult'            kept as-is (advisory GEANT4 multipolarity)\n";
  }

  desc +=
  "\nField meanings (for readers not familiar with nuclear decay data):\n"
  "  e   = energy in keV.\n"
  "  i   = intensity conditional on this decay transition; multiply by br for\n"
  "        probability per parent decay (<x> is relative to the strongest line).\n"
  "  hl  = half-life in seconds (\"INF\" = stable); iso = isomer number\n"
  "        (0 = ground state, 1 = first metastable isomer).\n"
  "  an  = atomic number Z; mn = mass number A.\n"
  "  br  = branch ratio: fraction of parent decays taking this transition;\n"
  "  productIntensityBasis=transitionConditional means reviewed EC/beta-plus/x-ray i values\n"
  "        are conditional on this transition and are multiplied by br once.\n"
  "  m   = decay mode (b-, b+, ec, alpha, it).\n"
  "  lft = logFT: log10 of the beta/EC comparative half-life (ft); a smaller\n"
  "        value means a more allowed (faster) branch.\n"
  "  f   = forbiddenness class (\"\" allowed, \"1U\" first forbidden unique).\n";

  if( cascade_detail_present )
  {
    desc +=
    "  Cascade detail (supports true-coincidence summing and post-decay x-rays):\n"
    "  ic  = total internal-conversion coefficient (IC electrons)/(gammas);\n"
    "        ik = K-shell part; il1,il2,il3 = L1/L2/L3 subshell parts (on <g>).\n"
    "  fl,tl = 'from level' / 'to level' indices in the GEANT4 progeny level\n"
    "        scheme that a gamma connects (the cascade topology, on <g>).\n"
    "  e0  = 1 only when a reviewed primary evaluation verifies a pure E0; it\n"
    "        is independent of the advisory GEANT4 multipolarity code.\n"
    "  ek,el = electron-capture K/L shell fractions; lv = progeny level fed by\n"
    "        an EC or beta-plus branch (on <ec> or <p>).\n"
    "  <lf l=... p=... q=... f=...> = complete alpha/beta-minus daughter-level\n"
    "        population law, conditional on the transition; entries include\n"
    "        ground l=0. q is the RDM final Q (beta endpoint; not alpha KE), and\n"
    "        optional f is the beta forbiddenness/shape class.\n"
    "  a2,a4 = gamma-gamma angular-correlation coefficients, on <c>, in\n"
    "        W(theta) = 1 + a2*P2(cos theta) + a4*P4(cos theta) (0 = isotropic).\n"
    "The full human-readable sandia.decay.xml has the long-form glossary and the\n"
    "per-transition \"cascade-match:\" audit comments (recording where the GEANT4\n"
    "cascade match was incomplete); those audit comments are stripped here.\n";
  }
  else
  {
    desc +=
    "  ic  = total internal-conversion coefficient (IC electrons)/(gammas);\n"
    "        ik = K-shell part; il1,il2,il3 = L1/L2/L3 subshell parts (on <g>).\n"
    "  mult = advisory GEANT4 multipolarity; it is not verified E0 provenance.\n"
    "  e0  = 1 only when a reviewed primary evaluation verifies a pure E0; it\n"
    "        is independent of the advisory GEANT4 multipolarity code.\n"
    "Gamma ICC, mult and e0 are retained because they control physical emission.\n"
    "Direct alpha/beta-minus level-feed laws are retained for the same reason.\n"
    "EC fed-level data, gamma topology/angular detail, and per-transition\n"
    "cascade-match audit comments are omitted from this coincidence-free file;\n"
    "EC shell fractions are retained because they control physical x-ray yield.\n"
    "those topology/audit details are available in the full sandia.decay.xml.\n";
  }

  const char *doc_desc = doc_node->document()->allocate_string(desc.c_str());
  rapidxml::xml_node<char> *desc_node = doc_node->document()->allocate_node(rapidxml::node_comment, "", doc_desc );
  doc_node->insert_node( doc_node->first_node(), desc_node );
}//void shrink_xml_elements( )


void remove_id_attrib( rapidxml::xml_node<char> *node )
{
  rapidxml::xml_attribute<char> *attrib = node->first_attribute("id");
  if( attrib )
    node->remove_attribute( attrib );
}


void check_revision_attrib( rapidxml::xml_attribute<char> *att )
{
  if( !att )
    return;
  const string val = att->value();
  if( val != "deleted" && val != "edited" && val != "inserted" && val != "edit" )
    throw runtime_error( "Invalid revision attribute value: '" + val + "'" );
}


string round_at_pos( string input, size_t pos )
{
  //pos indicates the first digit we do not want
  if (pos >= input.size())
    return input;
  
  assert((input[pos] >= '0') && (input[pos] <= '9'));
  
  for (size_t i = pos + 1; i < input.size(); ++i)
    input[i] = '0';
  
  // We'll attempt to do some rounding - man, is this horrible
  bool no_round = (input[pos] < '5');
  if (input[pos] == '5')
  {
    // round to nearest and ties to even
    char prevchar = '\0';
    if ((pos > 0) && (input[pos - 1] != '.'))
      prevchar = input[pos - 1];
    else if (pos > 1)
      prevchar = input[pos - 2];
    
    const int digval = prevchar - '0';
    no_round = !(digval % 2);
  }//if( input[pos] == '5' )
  
  
  if (no_round)
  {
    input = input.substr(0, pos);
  }
  else
  {
    input = input.substr(0, pos);
    
    const bool is_negative = (input[0] == '-');
    if (is_negative)
      input = input.substr(1);
    
    for (size_t index = input.size() - 1; index > 0; --index)
    {
      if (input[index] != '.')
      {
        assert(input[index] >= '0' && input[index] <= '9');
        
        if (input[index] != '9')
        {
          input[index] += 1;
          break;
        }
        
        assert(input[index] == '9');
        
        input[index] = '0';
      }//if( input[index] != '.' )
      
      if (index == 1)
      {
        if (input[0] == '.')
        {
          input = "1" + input;
        }
        else if (input[0] != '9')
        {
          input[0] += 1;
        }
        else
        {
          input[0] = '0';
          input = "1" + input;
          break;
        }
      }
    }//for( size_t index = decimal.size() - 1; index > 0; --index )
    
    if (is_negative)
      input = "-" + input;
  }//if( we need to round our last digit up )
  
  return input;
}//string round_at_pos( string input, size_t pos )


// A function to convert floats to their shortest text representation with at least the
//  specified precision (e.g., number of sig figs),
// Definitely not a very efficient function, or even optimal, but better than nothing.
string printCompact_trial(double val, size_t precision, const bool use_scientific)
{
  assert(precision > 0);
  if (!precision)
    precision = 1;

  //Note that using %.{precision}G will cause snprintf to use the
  //  IEEE 754 rounding rule of "round to nearest and ties to even"...

  char buffer[256] = { '\0' };
  if (use_scientific)
  {
    snprintf(buffer, sizeof(buffer), ("%." + to_string(precision) + "E").c_str(), val);
  }
  else
  {
    // We add 1 to precision to avoid double-rounding that leads to errors (I think 1 is
    //  probably enough, but we're not exactly being efficient anyway)
    //  Ex. 1.2345 to precision=3, if "%.3f" then will be "1.235", which we would then round to
    //      "1.24", which is not correct
    size_t ndec = precision + 2;
    if ((fabs(val) < 1.0) && (fabs(val) > std::numeric_limits<float>::min()))
      ndec = 2 + precision + std::ceil(fabs(std::log10(fabs(val))));

    snprintf(buffer, sizeof(buffer), ("%." + to_string(ndec) + "f").c_str(), val);
  }

  string power;
  string decimal = buffer;

  const string::size_type e_pos = decimal.find('E');
  if (e_pos != string::npos)
  {
    power = decimal.substr(e_pos);
    decimal = decimal.substr(0, e_pos);
  }//if( scientific notation )

  if (decimal.find('.') != string::npos)
  {
    while (decimal.size() && (decimal.back() == '0'))
      decimal.resize(decimal.size() - 1);

    if (decimal.size() && (decimal.back() == '.'))
      decimal.resize(decimal.size() - 1);
  }//if( decimal.find('.') != string::npos )



  bool hit_dec = false, hit_nonzero = false;
  size_t num_digit = 0;
  for (size_t i = 0; i < decimal.size(); ++i)
  {
    if (decimal[i] == '-')
      continue;

    if (decimal[i] == '.')
    {
      hit_dec = true;

      if ((num_digit >= precision) && hit_nonzero && ((i + 1) < decimal.size()))
      {
        decimal = round_at_pos(decimal, i + 1);
        break;
      }

      continue;
    }//if( decimal[i] == '.' )

    assert((decimal[i] >= '0') && (decimal[i] <= '9'));

    if (decimal[i] != '0')
      hit_nonzero = true;

    if (!hit_nonzero)
      continue;

    num_digit += 1;

    if (num_digit >= precision)
    {
      // If were at the last character, theres nothing to do
      if ((i + 1) == decimal.size())
        break;

      // If we're already past the decimal, get rid of everything past current character
      if (hit_dec)
      {
        decimal = round_at_pos(decimal, i + 1);
        break;
      }//if( hit_dec )
    }//if( we already have enough digits )
  }//for( loop over digits before "E" )



  // Remove trailing zeros after the decimal
  const bool has_per = (decimal.find('.') != string::npos);
  if (has_per)
  {
    while (decimal.size() && (decimal.back() == '0'))
      decimal.resize(decimal.size() - 1);

    // Dont leave a dangling decimal
    if (decimal.size() && decimal.back() == '.')
      decimal.resize(decimal.size() - 1);
  }//if( has_per )


  if (!power.empty())
  {
    power = power.substr(1);
    int intpow = atoi(power.c_str());

    if ((decimal.find('.') == string::npos) && (decimal.size() >= 2))
    {
      // decimal could now be "10", so we could increment the power, and remove a zero
      size_t ndecdigit = decimal.size();
      if (ndecdigit && (decimal[0] == '-'))
        ndecdigit -= 1;

      if ((ndecdigit > 0) && (decimal.back() == '0'))
      {
        intpow += 1;
        decimal.resize(decimal.size() - 1);
      }
    }//if( decimal is maybe reducable )

    power = "E" + std::to_string(intpow);
  }//if( !power.empty() )

  return decimal + power;
};//printCompact


string print_compact(const double val, const size_t precision)
{
  // From PhysicalUnits.cpp in InterSpec
  const string as_science = printCompact_trial(val, precision, true);
  const string as_fixed = printCompact_trial(val, precision, false);

  return ((as_fixed.size() <= as_science.size()) ? as_fixed : as_science);
}//string printCompact( double val, const size_t precision )


bool make_float_val_compact(rapidxml::xml_attribute<char>* attrib, const size_t precision)
{
  if (!attrib)
    return false;

  const string orig_val = attrib->value();
  if (orig_val.empty())
    return false;

  if (orig_val.size() <= precision)
    return true;

  double val;
  const int nread = sscanf(orig_val.c_str(), "%lf", &val);
  if (nread != 1)
    return false;

  if ((std::isnan)(val) || (std::isinf)(val))
    return false;

  string compacted = print_compact(val, precision);
  if (compacted.size() >= orig_val.size())
    return true;

  assert(attrib->document());
  const char* newval = attrib->document()->allocate_string(compacted.c_str());
  attrib->value(newval);
  
  return true;
}//make_float_val_compact(...)
