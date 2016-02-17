//
// file spiv_nogr_bits.cc
// David Cosgrove
// AstraZeneca
// 11th September 2006
//
// This is a collection of functions used by spiv and smg, ripped out of the
// original Spiv.cc

#include <string>
#include <vector>

#include <boost/bind.hpp>

#include <oechem.h>

#include "crash.H"
#include "stddefs.H"
#include "PharmPoint.H"

using namespace std;
using namespace OEChem;

// on a bigendian machine, spells Dave.
// used in MurmurHash2
static const int MAGIC_INT = 0x65766144;

namespace DACLIB {
  unsigned long SuperFastHash( const char *str , int l );
}
unsigned int MurmurHash2 ( const void * key, int len, unsigned int seed );

// ****************************************************************************
// do the expansion of any vector bindings
void expand_smarts_defs( const vector<pair<string,string> > &input_smarts ,
			 vector<pair<string,string> > &smarts_sub_defn ,
			 vector<pair<string,string> > &exp_smarts ) {

  string exp_smts;
  for( int i = 0 , is = input_smarts.size() ; i < is ; ++i ) {
    exp_smts = input_smarts[i].second;
    if( string::npos != exp_smts.find( "$" ) &&
	!OESmartsLexReplace( exp_smts , smarts_sub_defn ) ) {
      string msg = "Can't expand vector bindings in SMARTS string "
	+ input_smarts[i].second + " label " + input_smarts[i].first;
      exp_smarts.clear();
      throw( msg );
    }
    exp_smarts.push_back( make_pair( input_smarts[i].first , exp_smts ) );
  }

}

// ****************************************************************************
// build a list of unique names for the features found, with a count of each
void build_unique_names( const vector<vector<string> > &feature_names ,
			 map<string,int> &uniq_names ) {

  vector<vector<string> >::const_iterator r , rs;
  vector<string>::const_iterator q , qs;
  map<string,int>::iterator p;
  vector<string> mol_uniq_names;
  vector<string>::iterator t , ts;

  // build the list of all sites that appear in any molecule, counting how many
  // each occurs, just counting 1 occurrence per mol.
  for( r = feature_names.begin() , rs = feature_names.end() ; r != rs ; ++r ) {
    q = r->begin() + 1; // first name in vector is molecule name

    mol_uniq_names = vector<string>( r->begin() + 1 , r->end() );
    sort( mol_uniq_names.begin() , mol_uniq_names.end() );
    mol_uniq_names.erase( unique( mol_uniq_names.begin() , mol_uniq_names.end() ) ,
			  mol_uniq_names.end() );
    for( t = mol_uniq_names.begin() , ts = mol_uniq_names.end() ; t != ts ; ++t ) {
      p = uniq_names.find( *t );      
      if( p == uniq_names.end() )
	uniq_names.insert( make_pair( *t , 1 ) );
      else
	p->second++;
    }
  }

}

// ****************************************************************************
// convert a long feature name into a number, by hashing
string hash_feature_name( const string &fn ) {

#ifdef NOTYET
  unsigned int key = DACLIB::SuperFastHash( fn.c_str() , fn.length() );
#endif
  unsigned int key = MurmurHash2( fn.c_str() , fn.length() , MAGIC_INT );
  ostringstream oss;
  oss << hex << key;

  //  cout << oss.str() << " : " << key << " : " << fn << endl;

  return oss.str();

}

// ****************************************************************************
// build the short names from the unique names for the feature types, using
// SuperFastHash, warning of any collisions. Only do it for this names that
// match the min_occur criterion, as they'll be the only ones we're using.
// The returned vector has the short name first in each pair, the long name
// second.
void build_short_feature_names( const map<string,int> &uniq_names ,
				int min_occur , char feat_label ,
				vector<pair<string,string> > &short_names ) {

  short_names.clear();
  map<string,int>::const_iterator s , ss;
  for( s = uniq_names.begin() , ss = uniq_names.end() ; s != ss ; ++s ) {
    if( s->second >= min_occur ) {
      short_names.push_back( make_pair( hash_feature_name( s->first ) ,
					s->first ) );
    }
  }

  // now check for any collisions. Do it on a copy, we don't want to change the
  // order of the output vector
  vector<pair<string,string> > tmp_short_names( short_names );
  sort( tmp_short_names.begin() , tmp_short_names.end() ,
	boost::bind( less<string>() ,
		     boost::bind( &pair<string,string>::first , _1 ) ,
		     boost::bind( &pair<string,string>::first , _2 ) ) );

  vector<pair<string,string> >::iterator where = tmp_short_names.begin();
  vector<pair<string,string> >::iterator stop = tmp_short_names.end();

#ifdef NOTYET
  vector<pair<string,string> >::iterator r , rs;
  for( r = tmp_short_names.begin() , rs = tmp_short_names.end() ; r != rs ; ++r )
    cout << r->first << " : " << r->second << endl;
#endif

  while( ( where = adjacent_find( where , stop ,
				  boost::bind( equal_to<string>() ,
					       boost::bind( &pair<string,string>::first , _1 ) ,
					       boost::bind( &pair<string,string>::first , _2 ) ) ) ) != stop ) {
    cout << "AWOOGA Collision : " << where->second << " and " << (where-1)->second
	 << " have the same short name " << feat_label << where->first << endl;
    ++where;
  }

}

// ****************************************************************************
// write a file that decodes the bit headings, so as not to have the headings
// unfeasibly long
void write_name_decode_file( const string &decode_filename , char feat_label ,
			     const vector<pair<string,string> > &short_names ) {

  ofstream ofs1( decode_filename.c_str() );
  vector<pair<string,string> >::const_iterator s , ss;
  for( s = short_names.begin() , ss = short_names.end() ; s != ss ; ++s ) {
    ofs1 << feat_label << s->first << " " << s->second << endl;
  }

}

// ****************************************************************************
void write_bits_file( const string &output_filename , char feat_label ,
		      const vector<pair<string,string> > &short_names ,
		      const vector<vector<string> > &feature_names ) {

  // write the bit headings
  ofstream ofs2( output_filename.c_str() );
  ofs2 << "Molecule";
  vector<pair<string,string> >::const_iterator s , ss;
  for( s = short_names.begin() , ss = short_names.end() ; s != ss ; ++s ) {
    ofs2 << " " << feat_label << s->first;
  }
  ofs2 << endl;

  // and the bits
  vector<vector<string> >::const_iterator r , rs;
  vector<string>::const_iterator q , qs;
  for( r = feature_names.begin() , rs = feature_names.end(); r != rs ; ++r ) {
    ofs2 << r->front();
    vector<string> sorted_feature_names( r->begin() + 1 , r->end() );
    sort( sorted_feature_names.begin() , sorted_feature_names.end() );
    for( s = short_names.begin() , ss = short_names.end() ; s != ss ; ++s ) {
      if( std::binary_search( sorted_feature_names.begin() ,
			      sorted_feature_names.end() , s->second ) ) {
	ofs2 << " 1";
      } else {
	ofs2 << " 0";
      }
    }
    ofs2 << endl;
  }

}

// ****************************************************************************
void write_labels_file( const string &output_filename , char feat_label ,
			const vector<pair<string,string> > &short_names ,
			const vector<vector<string> > &feature_names ) {

  ofstream ofs2( output_filename.c_str() );

  vector<vector<string> >::const_iterator r , rs;
  vector<string>::const_iterator q , qs;
  vector<pair<string,string> >::const_iterator s , ss;
  for( r = feature_names.begin() , rs = feature_names.end(); r != rs ; ++r ) {
    cout << "Molecule name : " << r->front() << endl;
    ofs2 << r->front();
    vector<string> sorted_feature_names( r->begin() + 1 , r->end() );
    sort( sorted_feature_names.begin() , sorted_feature_names.end() );
    sorted_feature_names.erase( unique( sorted_feature_names.begin() ,
					sorted_feature_names.end() ) ,
				sorted_feature_names.end() );
    typedef vector<pair<string,string> >::const_iterator SNIter;
    typedef pair<SNIter,SNIter> SNIterPair;
    for( q = sorted_feature_names.begin() , qs = sorted_feature_names.end() ;
	 q != qs ; ++q ) {
      pair<string,string> search_pair = make_pair( "GARBAGE" , *q );
      SNIterPair snp = equal_range( short_names.begin() , short_names.end() ,
				    search_pair ,
				    boost::bind( less<string>() ,
						 boost::bind( &pair<string,string>::second , _1 ) ,
						 boost::bind( &pair<string,string>::second , _2 ) ) );
      if( snp.first == snp.second ) {
	cerr << "A major bad karma event - " << *q
	     << " didn't have a corresponding short name. Early bath indicated."
	     << endl;
	cout << "A major bad karma event - " << *q
	     << " didn't have a corresponding short name. Early bath indicated."
	     << endl;
	DACLIB::crash();
      }
      // probably if I had more faith in my understanding of the STL, I wouldn't
      // bother with this - it seems to come out first every time.
      if( snp.first->second == *q ) {
	ofs2 << " " << feat_label << snp.first->first;
      } else {
	ofs2 << " " << feat_label << snp.second->first;
      }
    }
    ofs2 << endl;
  }

}

// ****************************************************************************
void build_hash_codes( const vector<vector<string> > &all_names ,
		       char hash_prepend ,
		       vector<pair<string,string> > &names_to_hash ,
		       vector<pair<string,string> > &hash_to_names ) {

  set<string> uniq_names;
  vector<vector<string> >::const_iterator p , ps;
  vector<string>::const_iterator q , qs;
  for( p = all_names.begin() , ps = all_names.end() ; p != ps ; ++p ) {
    for( q = p->begin() , qs = p->end() ; q != qs ; ++q ) {
      uniq_names.insert( *q );
    }
  }
 
  set<string>::iterator r , rs;
  string shash_prepend;
  shash_prepend += hash_prepend;
  for( r = uniq_names.begin() , rs = uniq_names.end() ; r != rs ; ++r ) {
    string hash_code = shash_prepend + hash_feature_name( *r );
    names_to_hash.push_back( make_pair( *r , hash_code ) );
    hash_to_names.push_back( make_pair( hash_code , *r ) );
  }

  // names_to_hash should be in the right order, hash_to_names needs sorting
  std::sort( hash_to_names.begin() , hash_to_names.end() ,
	     boost::bind( less<string>() ,
			  boost::bind( &pair<string,string>::first , _1 ) ,
			  boost::bind( &pair<string,string>::first , _2 ) ) );

}

// ****************************************************************************
string name_from_pair_lookup( const string &search_name ,
			      const vector<pair<string,string> > &name_pairs ) {

  pair<string,string> search_pair = make_pair( search_name , "GARBAGE" );
  vector<pair<string,string> >::const_iterator lb =
    lower_bound( name_pairs.begin() , name_pairs.end() , search_pair ,
		 boost::bind( less<string>() ,
			      boost::bind( &pair<string,string>::first , _1 ) ,
			      boost::bind( &pair<string,string>::first , _2 ) ) );
  if( lb != name_pairs.end() && lb->first == search_name ) {
    return lb->second;
  } else {
    return "";
  }
			      
}

