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
#include <cctype>
#include <vector>
#include <memory>
#include <string>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <float.h>
#include <cstring>
#include <cstdlib>
#include <iostream>
#include <assert.h>

#include "SandiaDecay.h"

using namespace std;

//To compile on UNIX like systems:
//  g++ batch_decay_example.cpp SandiaDecay.cpp -o batch_decay_example.exe
//
// This program decays all nuclides either given on the command line, or input 
//  in a CSV/TSV file for a given time period, and outputs their activities
//  and optionally their proginies activities, either to stdout or a file, 
//  either for the time specified, or for the specified number of time steps


// poor-persons string trim
void trim( std::string &s )
{
  while( s.size() && std::isspace(s[0]) )
    s = s.substr(1);
  while( s.size() && std::isspace(s[s.size()-1]) )
    s = s.substr(0,s.size()-1);
}//trim(...)


void remove_chars( string &str, const string &chars )
{
  for( size_t i = 0; i < chars.size(); ++i )
  {
    const char c = chars[i];
    string::size_type pos = str.find(c);
    while( pos != string::npos )
      str.erase( str.begin() + pos );
  }
}//remove_chars(...)


bool starts_with( const std::string &line, const char *label )
{
  const size_t len1 = line.size();
  const size_t len2 = strlen(label);
  
  if( len1 < len2 )
    return false;
  
  for( size_t i = 0; i < len2; ++i )
    if( line[i] != label[i] )
      return false;
  
  return true;
}//starts_with(...)


void get_arg_value( const int argc, const char * argv[], const int index,
                    string &action, string &value )
{
  assert( index < argc );
  string arg = argv[index];
  if( arg.length() > 2 && arg[0]=='-' && arg[1]=='-' )
    arg = arg.substr(2);
  
  const string::size_type pos = arg.find( "=" );
  if( pos == string::npos )
  {
    action = arg;
    value = "";
  }else
  {
    action = arg.substr(0,pos);
    value = arg.substr(pos+1);
  }
}//get_arg_value(...)


string print_usage( const string exename )
{
  stringstream msg;
  msg << "Example batch decay calculator for SandiaDecay.\n"
         "\n"
         "Code licensed under the LGPL-2.1 license\n"
         "http://github.com/sandialabs/SandiaDecay\n"
         "\n"
         "Gives activities for nuclides, and optionally their progeny for a specified amount\n"
         " of decay time, and optionally for a number of steps between zero and then.\n"
         "\n"
         "Available command line arguments:" << endl;
  msg << "\tinput: The name of the CSV giving nuclide (first column), and initial activity (second column)\n"
       "\t\tInput may be either comma or tab delimited."<< endl;
  msg << "\toutput: The name of the CSV file to put output into.  If not specified, will output to stdout" << endl;
  msg << "\tdecay-data: The nuclear decay XML file to use.  Defaults to \"sandia.decay.xml\" in either\n"
       "\t\tthe current working directory (tried first), or the executables directory." << endl;
  msg << "\ttime: The time to decay nuclides too.  If no units are given, defaults to seconds, but will\n"
       "\t\talso accept units s (seconds), m (minutes), h (hours), d (days), y (years)." << endl;
  msg << "\tsteps: The number of steps between zero and time to print info out for.  Defaults to only\n"
       "\t\tprinting out the final activities (steps=1)." << endl;
  msg << "\tmix-input: if this argument is specified (no value can be given), then all input nuclides will\n"
       "\t\tbe mixed together, and co-decayed, and activities will be given for all children of the mixture.\n" << endl;
  msg << "\tshow-progeny: if this argument is specified (no value can be given), then activities of all\n"
       "\t\tprogeny will be given for each input nuclide." << endl;
  msg << "\tnuclide: Can be used to supply parent nuclides in addition to, or instead of CSV file.\n"
         "\t\tIf you follow the nuclide with a comma and then number, the number gives activity, or\n"
         "\t\telse defaults to 1" << endl;
  msg << "\tno-header: do not print the column header as first line of output\n";
  msg << "\thelp: show this message" << endl;
  
  msg << "\n"
         "\n"
         "Example uses:\n"
         "\t" << exename << " --input=input_nuclide.csv --output=output.csv --time=13\n"
         "\t\t(Produces output.csv with one line for each input line in input_nuclide.csv\n"
         "\t\t with the first column the input nuclide, and the second column the nuclides\n"
         "\t\t activity after 13 seconds)\n"
         "\n"
         "\t" << exename << " --input=input_nuclide.csv --output=\"output file.csv\" --time=13 --steps=20\n"
         "\t\t(similar to previous example, except will have 21 columns, the first giving the\n"
         "\t\t nuclide, and the following columns 20 giving the activity at 20 time points up\n"
         "\t\t through 14 seconds)\n"
         "\n"
         "\t" << exename << " --input=input_nuclide.csv --time=13y --show-progeny\n"
         "\t\t(For every input nuclide, will give the activity of input nuclide, and also for\n"
         "\t\t each descendant nuclide it decays to, after 13 years.  Output goes to stdout.)\n"
         "\n"
         "\t" << exename << " --input=input_nuclide.csv --time=13y --mix-input\n"
         "\t\t(Sums all input nuclides together into a mixture and outputs to stdout the\n"
         "\t\t activity of each input nuclide and their progeny after 13 years)\n"
         "\n"
         "\t" << exename << " --nuclide=I133,2.6 --nuclide=I131,5.9 --time=5d --mix-input\n"
         "\t\t(Creates a mixture of I-133 with activity 2.6 and I-131 with activity 5.9, and\n"
         "\t\t outputs their, and their progeny nuclides activities after 5 days)\n"
  << endl;
  
  return msg.str();
}//void print_usage()


struct InputLine
{
  string nuc_str;
  const SandiaDecay::Nuclide *nuclide;
  double start_activity;
};//struct InputLine


int main( int argc, const char * argv[] )
{
  try
  {
    // Very niavely get the executable name and path to make printing out a little nicer, and also
    //  to help search for sandia.decay.xml
    string exe_name = argv[0];
    string exe_path;
    const string::size_type sep_pos = exe_name.find_last_of( "/\\" );
    if( sep_pos != string::npos )
    {
      exe_path = exe_name.substr(0,sep_pos);
      exe_name = exe_name.substr(sep_pos+1);
    }
    
    if( argc < 2 )
    {
      cout << print_usage(exe_name) << endl;
      return EXIT_SUCCESS;
    }
    
    // Options that can be set from the command line
    size_t numhistories = 1;
    bool mix_input = false;
    bool show_children = false;
    bool print_header = true;
    string inputcsv, outputcsv;
    double decay_time = -1.0;
    
    
    // The input nuclides
    vector<InputLine> inputs;
    
    SandiaDecay::SandiaDecayDataBase database;
    SandiaDecay::NuclideMixture sum_mix;
    
    // First, we'll search for sandia.decay.xml, just in case user is setting it.
    string sandia_decay_xml = "sandia.decay.xml";
    for( int i = 1; i < argc; ++i )
    {
      string arg, value;
      get_arg_value( argc, argv, i, arg, value );
      if( arg == "decay-data" )
        sandia_decay_xml = value;
    }//for( int i = 1; i < argc; ++i )
    
    {//begin test if we can open sandia_decay_xml
      ifstream test( sandia_decay_xml.c_str() );
      if( !test.is_open() )
      {
#ifdef _WIN32
        const string testpath = exe_path + "\\" + sandia_decay_xml;
#else
        const string testpath = exe_path + "/" + sandia_decay_xml;
#endif
        cout << "Will test path " << testpath << endl;
        ifstream test2( testpath.c_str() );
        if( !test2.is_open() )
          throw runtime_error( "Could not open input nuclear data file (e.g. sandia.decay.xml)" );
        
        sandia_decay_xml = testpath;
      }//if( !test.is_open() )
    }//end test if we can open sandia_decay_xml
    
    // We have to initialize the database before parsing the rest of the command line arguments.
    database.initialize( sandia_decay_xml );
    
    // Parse command line arguments in a very niave way
    for( int i = 1; i < argc; ++i )
    {
      string arg, value;
      get_arg_value( argc, argv, i, arg, value );
      
      if( arg == "input" )
      {
        inputcsv = value;
        remove_chars( inputcsv, "\"'" );
      }else if( arg == "output" )
      {
        outputcsv = value;
        remove_chars( outputcsv, "\"'" );
      }else if( arg == "time" )
      {
        if( value.empty() )
          throw runtime_error( "No value for time was provided" );
        
        double time_mult = 1.0;
        const char last_char = value[value.size()-1];
        switch( last_char )
        {
          case 's': time_mult = 1.0; break;
          case 'm': time_mult = 60.0; break;
          case 'h': time_mult = 60.0*60.0; break;
          case 'd': time_mult = 24.0*60.0*60.0; break;
          case 'y': time_mult = 365.2425*24.0*60*60.0; break;
          case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
          case '8': case '9':
            time_mult = 1.0;
            break;
            
          default:
            throw runtime_error( "Invalid time given; only units allowed are s, m, h, d, y" );
        }//switch( last_char )
        
        if( value.size() >=2 )
        {
          switch( value[value.size()-2] )
          {
            case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7':
            case '8': case '9':
              break;
            default:
              throw runtime_error( "Invalid time units given; only units allowed are s, m, h, d, y" );
          }
        }//if( value.size() >=2 )
        
        
        if( !std::isdigit( value[value.size()-1] ) )
          value = value.substr(0, value.size()-1);
        
        if( !(stringstream(value) >> decay_time) )
          throw runtime_error( "Could not interpret time number" );
        
        decay_time *= time_mult;
      }else if( arg == "steps" )
      {
        if( !(stringstream(value) >> numhistories) )
          throw runtime_error( "Unable to convert argument for 'steps' to an integer." );
        if( numhistories < 1 )
          throw runtime_error( "You must specify at least one step." );
      }else if( arg == "mix-input" )
      {
        mix_input = true;
      }else if( arg == "show-progeny" )
      {
        show_children = true;
      }else if( arg == "nuclide" || arg == "isotope" || arg == "iso" )
      {
        InputLine line;
        line.start_activity = 1.0;
        
        const string::size_type pos = value.find( "," );
        if( pos == string::npos )
        {
          line.nuc_str = value;
        }else
        {
          line.nuc_str = value.substr(0,pos);
          if( !(stringstream(value.substr(pos+1)) >> line.start_activity) )
            throw runtime_error( "Invalid activity for nuclide '" + value + "'" );
        }
        
        line.nuclide = database.nuclide(line.nuc_str);
        if( !line.nuclide )
          throw runtime_error( "Nuclide '" + line.nuc_str + "' is not a valid nuclide" );
        
        inputs.push_back(line);
      }else if( arg == "help" )
      {
        cout << print_usage(exe_name) << endl;
        return EXIT_SUCCESS;
      }else if( arg == "no-header" )
      {
        print_header = false;
      }else if( arg == "decay-data" )
      {
        // already dealt with
      }else
      {
        throw runtime_error( "Unrecognized command line argument: '" + string(argv[i]) + "'" );
      }
    }//for( int i = 1; i < argc; ++i )
    
    if( inputs.empty() && inputcsv.empty() )
      throw runtime_error( "No input nuclides or input CSV file specified" );
    
    if( mix_input )
      show_children = true;
    
    if( decay_time < 0 )
      throw runtime_error( "No decay time was given." );
    
    std::auto_ptr<ifstream> input_file;
    std::auto_ptr<ofstream> output_file;
    if( !inputcsv.empty() )
    {
      input_file.reset( new ifstream( inputcsv.c_str() ) );
      if( !input_file->is_open() )
        throw runtime_error( "Could not open input file '" + inputcsv + "'" );
    }
    
    if( !outputcsv.empty() )
    {
      output_file.reset( new ofstream( outputcsv.c_str() ) );
      if( !output_file->is_open() )
        throw runtime_error( "Could not open final activity output file '" + outputcsv + "'" );
    }
  
    assert( numhistories > 0 );
    
    if( input_file.get() )
    {
      string line;
      while( getline(*input_file, line))
      {
        trim( line );
        
        if( line.size() < 2 ) //Allow (nearly) empty lines
          continue;
        
        if( line[0] == '#' )  //Allow lines to be commented with the hash character
          continue;
        
        const string::size_type delim_pos = line.find_first_of( ",\t" );
        if( delim_pos == string::npos )
          throw runtime_error( "Could not find a comma or tab for line '" + line + "'" );
        
        //Make sure only two columns
        const string::size_type second_delim = line.find_first_of( ",\t", delim_pos + 1 );
        if( second_delim != string::npos )
          throw runtime_error( "Line '" + line + "' had more than two columns" );
        
        string nuc_str = line.substr(0,delim_pos);
        string act_str = line.substr(delim_pos+1);
        
        trim( nuc_str );
        trim( act_str );
        
        InputLine input;
        input.nuc_str = nuc_str;
        input.nuclide = database.nuclide( nuc_str );
        //if( !input.nuclide )
        //  throw runtime_error( "Invalid nuclide specified: '" + nuc_str + "'" );
        
        if( !(stringstream( act_str ) >> input.start_activity) )
          throw runtime_error( "Unable to convert '" + act_str + "' to a number for activity" );
        
        if( input.start_activity < 0.0 )
          throw runtime_error( "Activity less than zero: '" + act_str + "'" );
        
        inputs.push_back( input );
      }//while( get lines )
    }//if( input_file.get() )
    
    ostream &output = outputcsv.empty() ? cout : *output_file;
    
    if( print_header )
    {
      if( mix_input || !show_children )
      {
        output << "Nuclide";
        for( size_t step = 0; step < numhistories; ++step )
        {
          const double t = ((numhistories == 1) ? decay_time : (step * decay_time / (numhistories - 1)));
          output << ",Act " << t << " seconds";
        }
        output << endl;
      }else
      {
        assert( show_children );
        output << "ParentNuclide,ProgenyNuclide";
        for( size_t step = 0; step < numhistories; ++step )
        {
          const double t = ((numhistories == 1) ? decay_time : (step * decay_time / (numhistories - 1)));
          output << ",Act " << t << " seconds";
        }
        output << endl;
      }//if( mix input or show children ) / else
    }//if( print_header )
    
    for( size_t i = 0; i < inputs.size(); ++i )
    {
      const InputLine &input = inputs[i];
      
      if( !input.nuclide )
      {
        cerr << "Warning: " << input.nuc_str << " does not appear to be a valid nuclide" << endl;
        if( !mix_input )
          output << input.nuc_str << ",InvalidNuclide" << endl;
        continue;
      }
      
      if( input.start_activity == 0.0 )
      {
        cerr << "Warning: zero initial activity for " << input.nuc_str << endl;
        if( !mix_input )
          output << input.nuc_str << "," << input.start_activity << endl;
        continue;
      }
      
      if( input.nuclide->isStable() )
      {
        if( !mix_input )
          output << input.nuc_str << ",Stable" << endl;
        continue;
      }
      
      if( mix_input )
      {
        sum_mix.addNuclideByActivity( input.nuclide, input.start_activity*SandiaDecay::curie );
        continue;
      }//if( mix_input )
      
      
      const double fake_activity = 0.01*SandiaDecay::curie;
      const double scale = input.start_activity / fake_activity;
      
      SandiaDecay::NuclideMixture mix;
      mix.addNuclideByActivity( input.nuclide, fake_activity );
      
      if( show_children )
      {
        output << input.nuc_str;
        const int nanswer = mix.numSolutionNuclides();
        for( int nuc_index = nanswer-1; nuc_index >= 0; --nuc_index )
        {
          const SandiaDecay::Nuclide *nuc = mix.solutionNuclide( nuc_index );
          
          output << "," << nuc->symbol;
          for( size_t step = 0; step < numhistories; ++step )
          {
            const double t = ((numhistories == 1) ? decay_time : (step * decay_time / (numhistories - 1)));
            const double fake_decayed_act = mix.activity(t, nuc_index);
            output << "," << fake_decayed_act*scale;
          }
          output << endl;
        }//for( int nuc_index = 0; nuc_index < nanswer; ++nuc_index )
      }else
      {
        output << input.nuc_str;
        for( size_t step = 0; step < numhistories; ++step )
        {
          //We could actually do the simple thing of
          //const double expCoeff = input.nuclide->decayConstant();
          //const double decayact = input.start_activity * exp( -decay_time * expCoeff );
          //cout << input.nuclide->symbol << endl;
          //output << input.nuc_str << "," << decayact << endl;
          
          const double t = ((numhistories == 1) ? decay_time : (step * decay_time / (numhistories - 1)));
          const double fake_decayed_act = mix.activity( t*SandiaDecay::second, input.nuclide );
          output << "," << fake_decayed_act*scale;
        }//for( size_t step = 0; step < numhistories; ++step )
        
        output << endl;
      }//if( show_children ) / else
      
    }//for( size_t i = 0; i < inputs.size(); ++i )
    
    
    if( mix_input )
    {
      const int nanswer = sum_mix.numSolutionNuclides();
      for( int nuc_index = nanswer-1; nuc_index >= 0; --nuc_index )
      {
        const SandiaDecay::Nuclide *nuc = sum_mix.solutionNuclide( nuc_index );
        
        output << nuc->symbol;
        for( size_t step = 0; step < numhistories; ++step )
        {
          const double t = ((numhistories == 1) ? decay_time : (step * decay_time / (numhistories - 1)));
          const double act = sum_mix.activity(t, nuc_index);
          output << "," << act/SandiaDecay::curie;
        }
        output << endl;
      }//for( loop over nuclides to list )
    }//if( numhistories > 0 )
  }catch( std::exception &e )
  {
    cerr << "Error: " << e.what() << endl;
    
    return EXIT_FAILURE;
  }//try / catch
  
  
  return EXIT_SUCCESS;
}//int main( int argc, const char * argv[] )

