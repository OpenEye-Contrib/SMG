//
// file spiv_nogr_bits.H
// David Cosgrove
// AstraZeneca
// 11th September 2006
//
// Declarations of functions in spiv_nogr_bits.cc

#include <string>
#include <vector>

#include <oechem.h>

class PharmPoint;

// do the expansion of any vector bindings
void expand_smarts_defs( const std::vector<std::pair<std::string,std::string> > &input_smarts ,
			 std::vector<std::pair<std::string,std::string> > &smarts_sub_defn ,
			 std::vector<std::pair<std::string,std::string> > &exp_smarts );

// make a map of OESubSearch objects corresponding to the input SMARTS,
// keyed on the SMARTS name.  The name is the first string in the pair,
// the SMARTS definition (with all vector bindings expanded) in the second.
void build_oesubsearches( PharmPoint &pharm_points ,
			  const std::vector<std::pair<std::string,std::string> > &smarts_defs ,
			  std::map<std::string,OEChem::OESubSearch *> &subs );

// build a list of unique names for the features found, with a count of each
void build_unique_names( const std::vector<std::vector<std::string> > &feature_names ,
			 std::map<std::string,int> &uniq_names );

// convert a long feature name into a number, by hashing
std::string hash_feature_name( const std::string &fn );

// build the short names from the unique names for the feature types, using
// SuperFastHash, warning of any collisions. Only do it for this names that
// match the min_occur criterion, as they'll be the only ones we're using.
// The returned vector has the short name first in each pair, the long name
// second.
void build_short_feature_names( const std::map<std::string,int> &uniq_names ,
				int min_occur , char feat_label ,
				std::vector<std::pair<std::string,std::string> > &short_names );

// write a file that decodes the bit headings, so as not to have the headings
// unfeasibly long
void write_name_decode_file( const std::string &decode_filename , char feat_label ,
			     const std::vector<std::pair<std::string,std::string> > &short_names );

void write_bits_file( const std::string &output_filename , char feat_label ,
		      const std::vector<std::pair<std::string,std::string> > &short_names ,
		      const std::vector<std::vector<std::string> > &feature_names );
void write_labels_file( const std::string &output_filename , char feat_label ,
			const std::vector<std::pair<std::string,std::string> > &short_names ,
			const std::vector<std::vector<std::string> > &feature_names );

// create lookup tables for hash codes for the feature names
void build_hash_codes( const std::vector<std::vector<std::string> > &all_names ,
		       char hash_prepend ,
		       std::vector<std::pair<std::string,std::string> > &names_to_hash ,
		       std::vector<std::pair<std::string,std::string> > &hash_to_names );

// find the search name in the pair, which is assumed to be sorted in order
// of something that search_name is in the set of. For looking up feature names
// given the hash code or vice versa.
std::string name_from_pair_lookup( const std::string &search_name ,
				   const std::vector<std::pair<std::string,std::string> > &name_pairs );
