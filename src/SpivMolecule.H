//
// file SpivMolecule.H
// David Cosgrove
// AstraZeneca
// 23rd August 2006
//
// This is the interface for the class SpivMolecule, which derives from an
// OEMol and holds the info for the pharmacophore sites, pairs and triplets

#ifndef DAC_SPIV_MOLECULE__
#define DAC_SPIV_MOLECULE__

#include <functional>
#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <oechem.h>

using namespace std;
using namespace OEChem;

// **********************************************************************

typedef struct {
  int site1_ , site2_;
  string site_label1_ , site_label2_;
  int dist_;
  string label_;
} SPIV_PAIR;

typedef struct {
  int sites_[3];
  string site_labels_[3];
  int dists_[3];
  string label_;
} SPIV_TRIPLET;

class SpivPairIsLess : public binary_function<SPIV_PAIR , SPIV_PAIR , bool> {
public :
  result_type operator()( first_argument_type a ,
			  second_argument_type b ) const {
    if( a.site_label1_ == b.site_label1_ ) {
      if( a.site_label2_ == b.site_label2_ )
	return a.dist_ < b.dist_;
      else
	return a.site_label2_ < b.site_label2_;
    } else
      return a.site_label1_ < b.site_label1_;
  }

};

class SpivTripletIsLess : public binary_function<SPIV_TRIPLET , SPIV_TRIPLET , bool> {
public :
  result_type operator()( first_argument_type a ,
			  second_argument_type b ) const {
    pair<string *,string *> res1 =
      mismatch( a.site_labels_ , a.site_labels_ + 3 , b.site_labels_ );
    if( res1.first == a.site_labels_ + 3 ) {
      pair<int *,int *> res2 = mismatch( a.dists_ , a.dists_ + 3 , b.dists_ );
      if( res2.first == a.dists_ + 3 )
	return false; // they're equal
      else
	return *(res2.first) < *(res2.second);
    } else
      return *(res1.first) < *(res1.second);
  }
};

class SpivPairIsSame : public binary_function<SPIV_PAIR , SPIV_PAIR , bool> {
public :
  result_type operator()( first_argument_type a ,
			  second_argument_type b ) const {
    return( a.site1_ == b.site1_ && a.site2_ == b.site2_ );
  }

};

class SpivTripletIsSame : public binary_function<SPIV_TRIPLET , SPIV_TRIPLET , bool> {
public :
  result_type operator()( first_argument_type a ,
			  second_argument_type b ) const {
    return equal( a.sites_ , a.sites_ + 3 , b.sites_ );
  }
};

// function to decide if a is longer than b, used when sorting the edges of
// a triplet.

class SpivPairIsLonger : public binary_function<SPIV_PAIR , SPIV_PAIR , bool> {
public :
  result_type operator()( first_argument_type a ,
			  second_argument_type b ) const {
    if( a.dist_== b.dist_ ) {
      if( a.site_label1_ == b.site_label1_ )
	return a.site_label2_ > b.site_label2_;
      else
	return a.site_label1_ > b.site_label1_;
    } else
      return a.dist_ > b.dist_;
  }

};

bool spiv_triplet_matches_criteria( const SPIV_TRIPLET &spiv_triplet ,
				    string *site_labels , int *min_dists ,
				    int *max_dists );

class PharmPoint;

class SpivMolecule : public OEMol {

public :

  SpivMolecule( OEMolBase &mol );
  virtual ~SpivMolecule();

  // make the 2D pharmacophore sites.
  void make_pphore_sites( PharmPoint &pharm_points ,
			  map<string,OESubSearch *> &oe_subs );
  void report_pphore_sites( ostream &os );
  void make_pphore_pairs();
  void make_pphore_triplets();

  const vector<vector<unsigned int> > &pphore_site_atoms() {
    return pphore_site_atoms_;
  }
  const vector<string> &pphore_site_labels() {
    return pphore_site_labels_;
  }
  const vector<SPIV_PAIR> &pphore_pairs() {
    return pphore_pairs_;
  }
  const vector<SPIV_TRIPLET> &pphore_triplets() {
    return pphore_triplets_;
  }

  // get the atoms that define the named feature. Empty vectors will be returned
  // if not relevant, e.g. if it's a Pairs feature, atoms3 will be empty. 
  void get_feature_atoms( const string &feature_type ,
			  const string &feature_name ,
			  vector<unsigned int> &atoms1 ,
			  vector<unsigned int> &atoms2 ,
			  vector<unsigned int> &atoms3 ) const;

  void make_atom_atom_dists_matrix();
  // find the shortest distance between an atom in site1 and an atom in
  // site2.
  int shortest_site_site_dist( int site1 , int site2 );

protected :

  vector<vector<unsigned int> > pphore_site_atoms_;
  vector<string>                pphore_site_labels_;
  vector<SPIV_PAIR>             pphore_pairs_;
  vector<SPIV_TRIPLET>          pphore_triplets_;

  int **atom_atom_dists_; // the distances between all atoms, shortest bond paths

  SPIV_PAIR make_spiv_pair( int site1 , int site2 );

  void get_sites_atoms( const string &feature_name ,
			vector<unsigned int> &atoms1 ) const;
  void get_pairs_atoms( const string &feature_name ,
			vector<unsigned int> &atoms1 ,
			vector<unsigned int> &atoms2 ) const;
  void get_triplets_atoms( const string &feature_name ,
			   vector<unsigned int> &atoms1 ,
			   vector<unsigned int> &atoms2 ,
			   vector<unsigned int> &atoms3 ) const;

};

#endif
