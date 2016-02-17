//
// file smg.cc
// Dave Cosgrove
// AstraZeneca
// 11th September 2006
//
// smg generates 2D pphore fingerprints as a batch process, with code ripped
// out of spiv.

#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>

#include "FileExceptions.H"
#include "SMARTSExceptions.H"
#include "PharmPoint.H"
#include "SpivMolecule.H"
#include "spiv_nogr_bits.H"

using namespace boost;
using namespace std;

typedef enum { SMG_UNDEFINED , SMG_SITES , SMG_PAIRS ,
	       SMG_TRIPLETS } SMG_OUTPUT_TYPE;
typedef enum { SMG_BITSTRINGS , SMG_LABELS } SMG_OUTPUT_FORMAT;

namespace DACLIB {
  void read_smarts_file( const string &smarts_file ,
			 vector<pair<string,string> > &input_smarts ,
			 vector<pair<string,string> > &smarts_sub_defn );
  void apply_daylight_aromatic_model( OEMolBase &mol );
}

// in spiv_nogr_bits.cc
void build_short_feature_names( const map<string,int> &uniq_names ,
				int min_occur , char feat_label ,
				vector<pair<string,string> > &short_names );
string hash_feature_name( const string &fn );

extern string BUILD_TIME; // in build_time.cc

// ***************************************************************************
void print_usage( ostream &os ) {

  os << "smg -mo[lecule_file] <string>" << endl
     << "    -sm[arts_file] <string>" << endl
     << "    -po[ints_file] <string>" << endl
     << "    -ou[tput_file] <string>" << endl
     << "    [-si[tes]" << endl
     << "    [-pa[irs]" << endl
     << "    [-t[riplets]" << endl
     << "    [-a[wk] <int>]" << endl
     << "    [-or[c] <int>]"
     << "    [-mi[n_dist] <int>]"
     << "    [-ma[x_dist] <int>]" << endl
     << "    [-b[itstrings]]" << endl
     << "    [-l[abels]]" << endl;

}

// ***************************************************************************
void parse_args( int argc , char **argv , string &mol_filename ,
		 string &smarts_filename , string &points_filename ,
		 string &output_filename , SMG_OUTPUT_TYPE &output_type ,
		 SMG_OUTPUT_FORMAT &output_format ,
		 int &min_occur , int &min_dist , int &max_dist ) {

  if( 1 == argc ) {
    print_usage( cout );
    exit( 0 );
  }
  output_type = SMG_UNDEFINED;
  output_format = SMG_BITSTRINGS;
  min_occur = -1;
  min_dist = 0;
  max_dist = 100;

  for( int i = 1 ; i < argc ; ++i ) {
    if( !strncmp( argv[i] , "-molecule_file" , 3 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-molecule_file requires a second argument.";
	exit( 1 );
      }
      mol_filename = argv[i];
    } else if( !strncmp( argv[i] , "-points_file" , 3 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-points_file requires a second argument.";
	exit( 1 );
      }
      points_filename = argv[i];
    } else if( !strncmp( argv[i] , "-smarts_file" , 3 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-smarts_file requires a second argument.";
	exit( 1 );
      }
      smarts_filename = argv[i];
    } else if( !strncmp( argv[i] , "-output_file" , 3 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-output_file requires a second argument.";
	exit( 1 );
      }
      output_filename = argv[i];
    } else if( !strncmp( argv[i] , "-orc" , 3 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-orc requires a second argument.";
	exit( 1 );
      }
      try {
	min_occur = lexical_cast<int>( argv[i] );
      } catch( bad_lexical_cast &e ) {
	cerr << "-orc requires an integer argument." << endl;
	exit( 1 );
      }
    } else if( !strncmp( argv[i] , "-awk" , 2 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-awk requires a second argument.";
	exit( 1 );
      }
      try {
	min_occur = lexical_cast<int>( argv[i] );
      } catch( bad_lexical_cast &e ) {
	cerr << "-awk requires an integer argument." << endl;
	exit( 1 );
      }
    } else if( !strncmp( argv[i] , "-min_dist" , 3 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-min_dist requires a second argument.";
	exit( 1 );
      }
      try {
	min_dist = lexical_cast<int>( argv[i] );
      } catch( bad_lexical_cast &e ) {
	cerr << "-min_dist requires an integer argument." << endl;
	exit( 1 );
      }
    } else if( !strncmp( argv[i] , "-max_dist" , 3 ) ) {
      ++i;
      if( i == argc ) {
	cerr << "-max_dist requires a second argument.";
	exit( 1 );
      }
      try {
	max_dist = lexical_cast<int>( argv[i] );
      } catch( bad_lexical_cast &e ) {
	cerr << "-max_dist requires an integer argument." << endl;
	exit( 1 );
      }
    } else if( !strncmp( argv[i] , "-help" , 2 ) ) {
      print_usage( cout );
      exit( 0 );
    } else if( !strncmp( argv[i] , "-sites" , 3 ) ) {
      output_type = SMG_SITES;
    } else if( !strncmp( argv[i] , "-pairs" , 3 ) ) {
      output_type = SMG_PAIRS;
    } else if( !strncmp( argv[i] , "-triplets" , 2 ) ) {
      output_type = SMG_TRIPLETS;
    } else if( !strncmp( argv[i] , "-bitstrings" , 2 ) ) {
      output_format = SMG_BITSTRINGS;
    } else if( !strncmp( argv[i] , "-labels" , 2 ) ) {
      output_format = SMG_LABELS;
    }
  }

  if( mol_filename.empty() ) {
    cerr << "No molecule file specfied." << endl;
    print_usage( cerr );
    exit( 1 );
  }
  if( smarts_filename.empty() ) {
    cerr << "No SMARTS file specfied." << endl;
    print_usage( cerr );
    exit( 1 );
  }
  if( points_filename.empty() ) {
    cerr << "No points file specfied." << endl;
    print_usage( cerr );
    exit( 1 );
  }
  if( output_filename.empty() ) {
    cerr << "No output file specfied." << endl;
    print_usage( cerr );
    exit( 1 );
  }
  if( SMG_UNDEFINED == output_type ) {
    cerr << "No output format specified (sites, pairs or triplets)." << endl;
    print_usage( cerr );
    exit( 1 );
  }

}

// ***************************************************************************
void read_molecules( const string &mol_filename ,
		     vector<SpivMolecule *> &all_mols ) {

  oemolistream ims( mol_filename.c_str() );
  if( !ims ) {
    throw string( "File " + mol_filename + " could not be read." );
  }
  OEMol oemol;
  while( ims >> oemol ) {
    OEAssignAromaticFlags( oemol , OEAroModelDaylight );
    all_mols.push_back( new SpivMolecule( oemol ) );
    oemol.Clear();
    if( !(all_mols.size() % 1000 ) ) {
      cout << "Read " << all_mols.size() << " molecules." << endl;
    }
  }

}

// ***************************************************************************
void write_feature_bits( char feat_label , const string &output_filename ,
			 const vector<vector<string> > &feature_names ,
			 int min_occur , SMG_OUTPUT_FORMAT output_format ,
			 map<string,int> &uniq_names ) {

  // build a list of unique names for the features found, with a count of each
  build_unique_names( feature_names , uniq_names );

  // build the short names from the unique names for the feature types, using
  // SuperFastHash, warning of any collisions. Only do it for this names that
  // match the min_occur criterion, as they'll be the only ones we're using.
  // The returned vector has the short name first in each pair, the long name
  // second.
  vector<pair<string,string> > short_names;
  build_short_feature_names( uniq_names , min_occur , feat_label , short_names );

  // write a file that decodes the bit headings, so as not to have the headings
  // unfeasibly long. Do this by generating hash codes. This may not give
  // unique names, so need to warn of collisions.
  string decode_filename = output_filename + ".name_decode";
  write_name_decode_file( decode_filename , feat_label , short_names );

#ifdef NOTYET  
  for( s = uniq_names.begin() , ss = uniq_names.end() ; s != ss ; ++s )
    cout << s->first << " : " << s->second << endl;
#endif

  if( output_format == SMG_BITSTRINGS ) {
    write_bits_file( output_filename , feat_label , short_names , feature_names );
  } else if( output_format == SMG_LABELS ) {
    write_labels_file( output_filename , feat_label , short_names ,
		       feature_names );
  }

}
      
// ***************************************************************************
void extract_feature_names( SpivMolecule &mol , SMG_OUTPUT_TYPE output_type ,
			    int min_dist , int max_dist ,
			    vector<string> &feat_names ) {

  feat_names.push_back( mol.GetTitle() );
  if( SMG_SITES == output_type ) {
    feat_names.insert( feat_names.end() ,
		       mol.pphore_site_labels().begin() ,
		       mol.pphore_site_labels().end() );
  } else if( SMG_PAIRS == output_type ) {
    const vector<SPIV_PAIR> &pairs = mol.pphore_pairs();
    for( int j = 0 , js = pairs.size() ; j < js ; ++j ) {
      if( pairs[j].dist_ >= min_dist && pairs[j].dist_ <= max_dist ) {
	feat_names.push_back( pairs[j].label_ );
      }
    }
  } else if( SMG_TRIPLETS == output_type ) {
    const vector<SPIV_TRIPLET> &trips = mol.pphore_triplets();
    for( int j = 0 , js = trips.size() ; j < js ; ++j ) {
      if( trips[j].dists_[0] >= min_dist && trips[j].dists_[0] <= max_dist &&
	  trips[j].dists_[1] >= min_dist && trips[j].dists_[1] <= max_dist &&
	  trips[j].dists_[2] >= min_dist && trips[j].dists_[2] <= max_dist ) {
	feat_names.push_back( trips[j].label_ );
      }
    }
  }

  // first entry is the molecule name, which needs to stay put
  sort( feat_names.begin() + 1 , feat_names.end() );
  feat_names.erase( unique( feat_names.begin() + 1 , feat_names.end() ) ,
		    feat_names.end() );

}

// ***************************************************************************
void write_output( SMG_OUTPUT_FORMAT output_format ,
		   SMG_OUTPUT_TYPE output_type ,
		   const string &output_filename ,
		   const vector<vector<string> > &feature_names ,
		   int min_occur ,
		   map<string,int> &unique_names ) {

  if( SMG_SITES == output_type ) {
    write_feature_bits( 'S' , output_filename , feature_names , min_occur ,
			output_format , unique_names );
    
  } else if( SMG_PAIRS == output_type ) {
    write_feature_bits( 'P' , output_filename , feature_names , min_occur ,
			output_format , unique_names );
  } else if( SMG_TRIPLETS == output_type ) {
    write_feature_bits( 'T' , output_filename , feature_names , min_occur ,
			output_format , unique_names );
  }

}

// ***************************************************************************
// final check of unique_names for the defitive collision check. We're not
// interested in the final answer, just what appears along the way.
void check_hash_collisions( const map<string,int> &unique_names ) {

  vector<pair<string,string> > short_names;
  cout << "Final check of collisions in all bit label names." << endl;
  build_short_feature_names( unique_names , 0 , 'X' , short_names );

}

// ***************************************************************************
int main( int argc , char **argv ) {

  string mol_filename , smarts_filename , points_filename;
  string output_filename;
  SMG_OUTPUT_TYPE output_type;
  SMG_OUTPUT_FORMAT output_format;

  int    min_occur; /* set by -awk or -orc, minimum number of instances of feature
		       for it to be written */
  int    min_dist , max_dist; /* min and max bond distances for output */

  cerr << "smg : "
       << BUILD_TIME << " using OEToolits version "
       << OEChem::OEChemGetRelease() << "." << endl;

  parse_args( argc , argv , mol_filename , smarts_filename , points_filename ,
	      output_filename , output_type , output_format ,
	      min_occur , min_dist , max_dist );

  vector<pair<string,string> > input_smarts , smarts_sub_defn , exp_smarts;

  try {
    DACLIB::read_smarts_file( smarts_filename , input_smarts , smarts_sub_defn );
  } catch( DACLIB::FileReadOpenError &e ) {
    cout << e.what() << endl;
    cerr << e.what() << endl;
  } catch( DACLIB::SMARTSSubDefnError &e ) {
    cout << e.what() << endl;
    cerr << e.what() << endl;
  } catch( DACLIB::SMARTSFileError &e ) {
    cout << e.what() << endl;
    cerr << e.what() << endl;
  }

  expand_smarts_defs( input_smarts , smarts_sub_defn , exp_smarts );

  PharmPoint pharm_points;
  try {
    pharm_points.read_points_file( points_filename );
  } catch( DACLIB::FileReadOpenError &e ) {
    cout << e.what() << endl;
    cerr << e.what() << endl;
    exit( 1 );
  }

  map<string,OESubSearch *> subs;
  build_oesubsearches( pharm_points , exp_smarts , subs );

  vector<vector<string> > feature_names;

  oemolistream ims( mol_filename.c_str() );
  if( !ims ) {
    throw string( "File " + mol_filename + " could not be read." );
  }

  OEMol oemol;
  int mol_count = 0;
  int file_num = 0;
  map<string,int> unique_names; // count of all long bit labels found.
  while( ims >> oemol ) {
    DACLIB::apply_daylight_aromatic_model( oemol );
    boost::scoped_ptr<SpivMolecule> spiv_mol( new SpivMolecule( oemol ) );
    try {
      spiv_mol->make_pphore_sites( pharm_points , subs );
    } catch( string msg ) {
      cout << msg << endl;
      exit( 1 );
    }
    vector<string> feat_names;
    if( SMG_PAIRS == output_type ) {
      spiv_mol->make_pphore_pairs();
    } else if( SMG_TRIPLETS == output_type ) {
      spiv_mol->make_pphore_triplets();      
    }
    extract_feature_names( *spiv_mol , output_type , min_dist , max_dist ,
			   feat_names );
    feature_names.push_back( feat_names );
    ++mol_count;
    if( ( ( mol_count < 5000 && !( mol_count % 100 ) ) ||
	  ( mol_count < 50000 && !( mol_count % 1000 ) ) ||
	  ( mol_count > 50000 && !( mol_count % 10000 ) ) ) )
      cerr << "Processed " << mol_count << " molecules." << endl;
    // if doing labels output, dump the results out every 200000 molecules.
    // Can't do same for bitstrings, so it will probably run out of memory at
    // some point for large data sets.
    if( !( mol_count % 200000 ) && SMG_LABELS == output_format ) {
      string tmp_file_name = output_filename + string( "." ) +
	boost::lexical_cast<string>( file_num );
      write_output( output_format , output_type , tmp_file_name ,
		    feature_names , min_occur , unique_names );
      feature_names.clear();
      ++file_num;
    }
  }

  if( SMG_BITSTRINGS == output_format || !file_num ) {
    write_output( output_format , output_type , output_filename ,
		  feature_names , min_occur , unique_names );
  } else if( SMG_LABELS == output_format && file_num ) {
    // finish off last ones
    string tmp_file_name = output_filename + string( "." ) +
      boost::lexical_cast<string>( file_num );
    cout << "Writing final part of output to " << tmp_file_name << endl;
    write_output( output_format , output_type , tmp_file_name ,
		  feature_names , min_occur , unique_names );
    // final check of collisions for unique names. If the file is written
    // out in bits, the interim reports may not be complete.
    check_hash_collisions( unique_names );
  }

}
