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

#include "SandiaDecay.h"

#include "rapidxml/rapidxml.hpp"
#include "rapidxml/rapidxml_utils.hpp"
#include "rapidxml/rapidxml_print.hpp"


using namespace std;


//Usuage:
//XML_FOREACH_DAUGHTER( child_node_variable, parent_node, "ChildElementName" ){
//  assert( child_node_variable->name() == string("ChildElementName") );
// }
#define XML_FOREACH_DAUGHTER( nodename, parentnode, daughternamestr ) \
for( rapidxml::xml_node<char> *nodename = parentnode->first_node(daughternamestr); \
  nodename; \
  nodename = nodename->next_sibling(daughternamestr) )

/** This code is just boiler-plate code to make mass-edits to sandia.decay.xml.
 Currently divides the gamma intensities of Th227 by 1.3565891473 (which our previous values were consitently off relative to JENDL-5)
 e.g., see https://www-nds.iaea.org/exfor/servlet/E4sShowDecayData?db=e4&op=showDecayData&req=3313&Sect=8937092,9150214,14660924&plus=3
 (the alphas were not off)
 */

int main( int argc, char **argv )
{
  if( argc != 3 )
  {
    cerr << "Usage: " << argv[0] << " <input xml> <output xml>" << endl;
    return 1;
  }
  
  try
  {
    rapidxml::file<char> xmlfile( argv[1] );
    
    rapidxml::xml_document<char> doc;
    doc.parse<rapidxml::parse_full>( xmlfile.data() );
    rapidxml::xml_node<char> *doc_node = doc.first_node();
  
    if( doc_node->first_attribute("version") )
      doc_node = doc_node->next_sibling();
    
    int ntrans = 0, nodes_removed = 0;
    
    vector<rapidxml::xml_node<char> *> nuc_to_delete;
    for( rapidxml::xml_node<char> *nuc = doc_node->first_node( "nuclide" );
        nuc;
        nuc = nuc->next_sibling("nuclide") )
    {
      // ...
    }//for( loop over nuclides )
    
    
    for( rapidxml::xml_node<char> *xmltrans = doc_node->first_node( "transition" );
        xmltrans;
        xmltrans = xmltrans->next_sibling("transition") )
    {
      rapidxml::xml_attribute<char> *child_label_attrib = xmltrans->first_attribute("child");
      rapidxml::xml_attribute<char> *parent_label_attrib = xmltrans->first_attribute("parent");
      
      const string child_nuc = child_label_attrib ? child_label_attrib->value() : "";
      const string parent_nuc = parent_label_attrib ? parent_label_attrib->value() : "";
      
      if( (parent_nuc != "Th227") || (child_nuc != "Ra223") )
        continue;
      

      for( rapidxml::xml_node<char> *child = xmltrans->first_node("gamma");
          child;
          child = child->next_sibling("gamma") )
      {
        rapidxml::xml_attribute<char> *attrib = child->first_attribute("intensity");
        const string valstr = attrib->value();
        
        double val;
        const int nread = sscanf(valstr.c_str(), "%lf", &val);
        
        if (nread != 1)
          throw runtime_error( "Error reading intensity value '" + valstr + "'" );
        
        const double fixed_val = val / 1.3565891473;
        char buffer[64] = { '\0' };
        snprintf( buffer, sizeof(buffer), "%.6G", fixed_val );
        
        cout << "Setting " << val << " --> " << fixed_val << endl;
        const char *newvalstr = doc.allocate_string( buffer );
        attrib->value( newvalstr );
        
        rapidxml::xml_node<char> *edit_node = doc.allocate_node( rapidxml::node_element, "edit" );
        rapidxml::xml_attribute<char> *edit_attrib = doc.allocate_attribute("dateTime", "2023-02-10T23:52:00Z");
        edit_node->append_attribute( edit_attrib );
        
        edit_attrib = doc.allocate_attribute("newValue", newvalstr);
        edit_node->append_attribute( edit_attrib );
        
        const char *oldvalstr = doc.allocate_string( valstr.c_str() );
        edit_attrib = doc.allocate_attribute("oldValue", oldvalstr);
        edit_node->append_attribute( edit_attrib );
        
        edit_attrib = doc.allocate_attribute("property", "intensity");
        edit_node->append_attribute( edit_attrib );
        
        edit_attrib = doc.allocate_attribute("remarks", "Fix incorrect Th232 scale factor");
        edit_node->append_attribute( edit_attrib );
        
        edit_attrib = doc.allocate_attribute("userUri", "mailto:wcjohns@sandia.gov");
        edit_node->append_attribute( edit_attrib );
        
        child->append_node( edit_node );
      }//for( loop over <gamma> elements
      
  
      for( rapidxml::xml_node<char> *child = xmltrans->first_node();
          child;
          child = child->next_sibling() )
      {
        //...
      }//for( loop over children )
    
      rapidxml::xml_attribute<char> *attrib = xmltrans->first_attribute("remarks");
      if( attrib )
      {
        string val = attrib->value();
        if( (val.find("warning:") == string::npos) && (val.find("audit:") == string::npos) )
          xmltrans->remove_attribute( attrib );
      }
      
      ++ntrans;
    }//loop over transitions.
    
  
    
    {
      ifstream test( argv[2] );
      if( test.is_open() )
        throw runtime_error( "File '" + string(argv[2]) + " already exists - not overwriting!" );
    }
    
    {
      ofstream outputfile( argv[2], ios::out | ios::binary );
      if( !outputfile )
        throw runtime_error( "Could not open output file: " + string(argv[2]) );
    
      rapidxml::print<char>( outputfile, doc, 0 );
    }
    
    printf( "Saved output file '%s'\n", argv[2] );
    
    
    
  }catch( std::exception &e )
  {
    fprintf( stderr, "Error: %s\n", e.what() );
  }//try / catch
  
  
  
  return EXIT_SUCCESS;
}//int main( int argc, char **argv )
