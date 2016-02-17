//
// file SpivMolecule.cc
// David Cosgrove
// AstraZeneca
// 23rd August 2006
//
// Implementation of SpivMolecule

#include <algorithm>

#include "stddefs.H"
#include "PharmPoint.H"
#include "SpivMolecule.H"

using namespace OESystem;
using namespace OEPlatform;

// ***********************************************************************
SpivMolecule::SpivMolecule( OEMolBase &mol ) : OEMol( mol ) {

  atom_atom_dists_ = 0;

}

// ***********************************************************************
SpivMolecule::~SpivMolecule() {

  DACLIB::destroy_square_matrix( atom_atom_dists_ );

}

// ***********************************************************************
// make the 2D pharmacophore sites.
void SpivMolecule::make_pphore_sites( PharmPoint &pharm_points ,
				      map<string,OESubSearch *> &oe_subs ) {

  pphore_site_atoms_.clear();
  pphore_site_labels_.clear();

  map<string,vector<string> > &points_defs = pharm_points.points_defs();
  map<string,vector<string> >::iterator p , ps;
  vector<string>::iterator q;
  map<string,OESubSearch *>::iterator r;

  OESubSearch *subs;
  OEIter<OEMatchBase> match;
  vector<unsigned int> next_ats;
  for( p = points_defs.begin() , ps= points_defs.end() ; p != ps ; ++p ) {
    //    cout << "Point type " << p->first << endl;
    if( p->second.empty() )
      continue; // point defined by key word (e.g. ITMOC, ITMOC_ALO) not SMARTS.

    for( q = p->second.begin() ; q != p->second.end() ; ++q ) {

      r = oe_subs.find( *q );
      if( r == oe_subs.end() ) {
	string msg = "No SMARTS definition for " + *q +
	  " needed for point type " + p->first;
	throw( msg );
      }
      subs = r->second;
      for( match = subs->Match( *this , true ) ; match ; ++match ) {
	OEIter<OEMatchPair<OEAtomBase> > mp = match->GetAtoms();
	next_ats.clear();
	for( ; mp ; ++mp )
	  next_ats.push_back( mp->target->GetIdx() );
	pphore_site_atoms_.push_back( next_ats );
	pphore_site_labels_.push_back( p->first );
      }
    }
  }

  //  report_pphore_sites( cout );

}

// ***********************************************************************
void SpivMolecule::report_pphore_sites( ostream &os ) {

  for( int i = 0 , is = pphore_site_atoms_.size() ; i < is ; ++i ) {
    os << pphore_site_labels_[i] << " : ";
    for( int j = 0 , js = pphore_site_atoms_[i].size() ; j < js ; ++j )
      os << pphore_site_atoms_[i][j] << " ";
    os << endl;
  }

}

// ***********************************************************************
// sites must be made before pairs, and sites need the SMARTS defs so can't
// be generated in this function.
void SpivMolecule::make_pphore_pairs() {

  if( pphore_site_labels_.empty() )
    return; // need sites for the pairs

  if( !atom_atom_dists_ )
    make_atom_atom_dists_matrix();

  pphore_pairs_.clear();

  //  cout << "XXXXXXXXXXXXXXXXXXXXXXXX" << GetTitle() << endl;

  for( int i = 0 , is = pphore_site_labels_.size() - 1 ; i < is ; ++i ) {
    for( int j = i + 1 , js = pphore_site_labels_.size() ; j < js ; ++j ) {
      pphore_pairs_.push_back( make_spiv_pair( i , j ) );
    }
  }

  sort( pphore_pairs_.begin() , pphore_pairs_.end() , SpivPairIsLess() );
  pphore_pairs_.erase( unique( pphore_pairs_.begin() , pphore_pairs_.end() ,
			       SpivPairIsSame() ) , pphore_pairs_.end() );

#ifdef NOTYET
  for( int i = 0 , is = pphore_pairs_.size() ; i < is ; ++i )
    cout << pphore_pairs_[i].label_ << " : "
	 << pphore_pairs_[i].site1_ << " - " << pphore_pairs_[i].site2_
	 << " : " << pphore_pairs_[i].dist_ << endl;
#endif

}

// ***********************************************************************
// pairs must be made before triplets
void SpivMolecule::make_pphore_triplets() {

  if( pphore_pairs_.empty() )
    make_pphore_pairs();
  if( pphore_pairs_.size() < 3 )
    return; // need 3 pairs for the triplets

  pphore_triplets_.clear();

  // the triplets are encoded using the algorithm of Abrahamian et al.
  // (paper 273, JCICS, 43, 458-468). The three features are labelled
  // f1, f2 and f3.  f2 is the feature common to the longest and shortest
  // edges, f1 is at the other end of the longest edge, f3 the other end of
  // the shortest edge.  If two edges have the same distance, priority is
  // given to the one with the higher label (label1>label2).
  SPIV_PAIR triplet_pairs[3];
  SPIV_TRIPLET spiv_triplet;
  for( int i = 0 , is = pphore_site_labels_.size() - 2 ; i < is ; ++i ) {
    for( int j = i + 1 , js = pphore_site_labels_.size() - 1 ; j < js ; ++j ) {
      for( int k = j + 1 , ks = pphore_site_labels_.size() ; k < ks ; ++k ) {

	triplet_pairs[0] = make_spiv_pair( i , j );
	triplet_pairs[1] = make_spiv_pair( j , k );
	triplet_pairs[2] = make_spiv_pair( i , k );

	sort( triplet_pairs , triplet_pairs + 3 , SpivPairIsLonger() );
	// spiv_triplet.sites_[1] is site in common for sides 2 and 0 (longest
	// and shortest).  spiv_triplet.sites_[0] is other end of side 2,
	// spiv_triplet.sites_[2] is other end of site 0
	if( triplet_pairs[2].site1_ == triplet_pairs[0].site1_ ) {
	  spiv_triplet.sites_[1] = triplet_pairs[2].site1_;
	  spiv_triplet.sites_[0] = triplet_pairs[2].site2_;
	} else if( triplet_pairs[2].site1_ == triplet_pairs[0].site2_ ) {
	  spiv_triplet.sites_[1] = triplet_pairs[2].site1_;
	  spiv_triplet.sites_[0] = triplet_pairs[2].site2_;
	} else if( triplet_pairs[2].site2_ == triplet_pairs[0].site1_ ) {
	  spiv_triplet.sites_[1] = triplet_pairs[2].site2_;
	  spiv_triplet.sites_[0] = triplet_pairs[2].site1_;
	} else if( triplet_pairs[2].site2_ == triplet_pairs[0].site2_ ) {
	  spiv_triplet.sites_[1] = triplet_pairs[2].site2_;
	  spiv_triplet.sites_[0] = triplet_pairs[2].site1_;
	}
	if( spiv_triplet.sites_[1] == triplet_pairs[0].site1_ )
	  spiv_triplet.sites_[2] = triplet_pairs[0].site2_;
	else
	  spiv_triplet.sites_[2] = triplet_pairs[0].site1_;

	spiv_triplet.dists_[0] = triplet_pairs[2].dist_;
	spiv_triplet.dists_[1] = triplet_pairs[1].dist_;
	spiv_triplet.dists_[2] = triplet_pairs[0].dist_;

	spiv_triplet.site_labels_[0] =
	  pphore_site_labels_[spiv_triplet.sites_[0]];
	spiv_triplet.site_labels_[1] =
	  pphore_site_labels_[spiv_triplet.sites_[1]];
	spiv_triplet.site_labels_[2] =
	  pphore_site_labels_[spiv_triplet.sites_[2]];

	spiv_triplet.label_ = triplet_pairs[2].label_ + "-" +
	  triplet_pairs[1].label_ + "-" + triplet_pairs[0].label_;

	pphore_triplets_.push_back( spiv_triplet );
      }
    }    
  }

  sort( pphore_triplets_.begin() , pphore_triplets_.end() ,
	SpivTripletIsLess() );
  pphore_triplets_.erase( unique( pphore_triplets_.begin() ,
				  pphore_triplets_.end() ,
				  SpivTripletIsSame() ) ,
			  pphore_triplets_.end() );

}

// ***********************************************************************
// get the atoms that define the named feature. Empty vectors will be returned
// if not relevant, e.g. if it's a Pairs feature, atoms3 will be empty. 
void SpivMolecule::get_feature_atoms( const string &feature_type ,
				      const string &feature_name ,
				      vector<unsigned int> &atoms1 ,
				      vector<unsigned int> &atoms2 ,
				      vector<unsigned int> &atoms3 ) const {

  atoms1.clear();
  atoms2.clear();
  atoms3.clear();

  if( string( "Sites" ) == feature_type ) {
    get_sites_atoms( feature_name , atoms1 );
  } else if( string( "Pairs" ) == feature_type ) {
    get_pairs_atoms( feature_name , atoms1 , atoms2 );
  } else if( string( "Triplets" ) == feature_type ) {
    get_triplets_atoms( feature_name , atoms1 , atoms2 , atoms3 );
  }

}

// ***********************************************************************
void SpivMolecule::make_atom_atom_dists_matrix() {

  if( atom_atom_dists_ )
    return; // only want to do it once

  const int max_idx = GetMaxAtomIdx();
  DACLIB::make_square_matrix( atom_atom_dists_ , max_idx );

  int half_max = numeric_limits<int>::max() / 2;
  fill( atom_atom_dists_[0] , atom_atom_dists_[0] + max_idx * max_idx ,
	half_max );
  for( int i = 0 ; i < max_idx ; ++i )
    atom_atom_dists_[i][i] = 0;

  // put the bonds in
  OEIter<OEAtomBase> atom , conns;
  for( atom = GetAtoms() ; atom ; ++atom ) {
    int this_idx = atom->GetIdx();
    for( conns = atom->GetAtoms() ; conns ; ++conns )
      atom_atom_dists_[this_idx][conns->GetIdx()] =
	atom_atom_dists_[conns->GetIdx()][this_idx] = 1;
  }

  // now use Floyd's algorithm to find the minimum through-bond distance
  // between each pair of atoms. This is taken from 'Practical Algorithms
  // in C++' by Flamig p 386.
  for( int k = 0 ; k < max_idx ; ++k ) {
    for( int i = 0 ; i < max_idx ; ++i ) {
      for( int j = 0 ; j < max_idx ; ++j ) {
	int old_dist = atom_atom_dists_[i][j];
	int new_dist = atom_atom_dists_[i][k] + atom_atom_dists_[k][j];
	if( new_dist < old_dist )
	  atom_atom_dists_[i][j] = new_dist;
      }
    }
  }

}

// *******************************************************************
// find the shortest distance between an atom in site1 and an atom in
// site2.
int SpivMolecule::shortest_site_site_dist( int site1 , int site2 ) {

  if( !atom_atom_dists_ )
    make_atom_atom_dists_matrix();

  int shortest_dist = numeric_limits<int>::max();
  vector<unsigned int> &i_atoms = pphore_site_atoms_[site1];
  vector<unsigned int> &j_atoms = pphore_site_atoms_[site2];
  for( int k = 0 , ks = i_atoms.size() ; k < ks ; ++k ) {
    for( int l = 0 , ls = j_atoms.size() ; l < ls ; ++l ) {
      if( atom_atom_dists_[i_atoms[k]][j_atoms[l]] < shortest_dist )
	shortest_dist = atom_atom_dists_[i_atoms[k]][j_atoms[l]];
    }
  }

  return shortest_dist;

}

// **************************************************************************
bool spiv_triplet_matches_criteria( const SPIV_TRIPLET &spiv_triplet ,
				    string *site_labels , int *min_dists ,
				    int *max_dists ) {

  // the distances are stored thus: site_labels[0] to site_labels[1]
  // dists min_dists[0] to max_dists[0], 1 to 2, min_dists[1], 2 to 0,
  // min_dists[2] but the sites won't necessarily be in the same order as
  // in the triplet, where dists_[0] is the longest dist, between 1 and 0,
  // dists_[2] is the shortest dist, between 1 and 2, and dists_[1] is the
  // other distance, between 0 and 2.  Need to deal with the 3 cyclic
  // permutations of this.
  // if the corners correspond, then the spiv_triplet.site_labels_ correspond
  // to min_dists and max_dists in order 0 => 0, 1 => 2, 2 => 1.

  if( ( site_labels[0] == "*" ||
	spiv_triplet.site_labels_[0] == site_labels[0] ) &&
      ( site_labels[1] == "*" ||
	spiv_triplet.site_labels_[1] == site_labels[1] ) &&
      ( site_labels[2] == "*" ||
	spiv_triplet.site_labels_[2] == site_labels[2] ) &&
      spiv_triplet.dists_[0] >= min_dists[0] &&
      spiv_triplet.dists_[0] <= max_dists[0] &&
      spiv_triplet.dists_[2] >= min_dists[1] &&
      spiv_triplet.dists_[2] <= max_dists[1] &&
      spiv_triplet.dists_[1] >= min_dists[2] &&
      spiv_triplet.dists_[1] <= max_dists[2] )
    return true;

  if( ( site_labels[2] == "*" ||
	spiv_triplet.site_labels_[0] == site_labels[2] ) &&
      ( site_labels[0] == "*" ||
	spiv_triplet.site_labels_[1] == site_labels[0] ) &&
      ( site_labels[1] == "*" ||
	spiv_triplet.site_labels_[2] == site_labels[1] ) &&
      spiv_triplet.dists_[0] >= min_dists[2] &&
      spiv_triplet.dists_[0] <= max_dists[2] &&
      spiv_triplet.dists_[2] >= min_dists[0] &&
      spiv_triplet.dists_[2] <= max_dists[0] &&
      spiv_triplet.dists_[1] >= min_dists[1] &&
      spiv_triplet.dists_[1] <= max_dists[1] )
    return true;

  if( ( site_labels[1] == "*" ||
	spiv_triplet.site_labels_[0] == site_labels[1] ) &&
      ( site_labels[2] == "*" ||
	spiv_triplet.site_labels_[1] == site_labels[2] ) &&
      ( site_labels[0] == "*" ||
	spiv_triplet.site_labels_[2] == site_labels[0] ) &&
      spiv_triplet.dists_[0] >= min_dists[1] &&
      spiv_triplet.dists_[0] <= max_dists[1] &&
      spiv_triplet.dists_[2] >= min_dists[2] &&
      spiv_triplet.dists_[2] <= max_dists[2] &&
      spiv_triplet.dists_[1] >= min_dists[0] &&
      spiv_triplet.dists_[1] <= max_dists[0] )
    return true;

  if( ( site_labels[0] == "*" ||
	spiv_triplet.site_labels_[0] == site_labels[0] ) &&
      ( site_labels[2] == "*" ||
	spiv_triplet.site_labels_[1] == site_labels[2] ) &&
      ( site_labels[1] == "*" ||
	spiv_triplet.site_labels_[2] == site_labels[1] ) &&
      spiv_triplet.dists_[0] >= min_dists[2] &&
      spiv_triplet.dists_[0] <= max_dists[2] &&
      spiv_triplet.dists_[2] >= min_dists[1] &&
      spiv_triplet.dists_[2] <= max_dists[1] &&
      spiv_triplet.dists_[1] >= min_dists[0] &&
      spiv_triplet.dists_[1] <= max_dists[0] )
    return true;

  if( ( site_labels[2] == "*" ||
	spiv_triplet.site_labels_[0] == site_labels[2] ) &&
      ( site_labels[1] == "*" ||
	spiv_triplet.site_labels_[1] == site_labels[1] ) &&
      ( site_labels[0] == "*" ||
	spiv_triplet.site_labels_[2] == site_labels[0] ) &&
      spiv_triplet.dists_[0] >= min_dists[1] &&
      spiv_triplet.dists_[0] <= max_dists[1] &&
      spiv_triplet.dists_[2] >= min_dists[0] &&
      spiv_triplet.dists_[2] <= max_dists[0] &&
      spiv_triplet.dists_[1] >= min_dists[2] &&
      spiv_triplet.dists_[1] <= max_dists[2] )
    return true;

  if( ( site_labels[1] == "*" ||
	spiv_triplet.site_labels_[0] == site_labels[1] ) &&
      ( site_labels[0] == "*" ||
	spiv_triplet.site_labels_[1] == site_labels[0] ) &&
      ( site_labels[2] == "*" ||
	spiv_triplet.site_labels_[2] == site_labels[2] ) &&
      spiv_triplet.dists_[0] >= min_dists[0] &&
      spiv_triplet.dists_[0] <= max_dists[0] &&
      spiv_triplet.dists_[2] >= min_dists[2] &&
      spiv_triplet.dists_[2] <= max_dists[2] &&
      spiv_triplet.dists_[1] >= min_dists[1] &&
      spiv_triplet.dists_[1] <= max_dists[1] )
    return true;

  return false;

}

// ***********************************************************************
SPIV_PAIR SpivMolecule::make_spiv_pair( int site1 , int site2 ) {

  SPIV_PAIR spiv_pair;

  static ostringstream oss;
  // want the shortest distance between an atom in site i and another atom
  // in site j
  int shortest_dist = shortest_site_site_dist( site1 , site2 );

  // the site names are created in alphabetical order, so this should always
  // create consistent labels, with the first site name <= second site name.
  if( site1 > site2 )
    std::swap( site1 , site2 );

  spiv_pair.site1_ = site1;
  spiv_pair.site2_ = site2;
  spiv_pair.site_label1_ = pphore_site_labels_[spiv_pair.site1_];
  spiv_pair.site_label2_ = pphore_site_labels_[spiv_pair.site2_];
  spiv_pair.dist_ = shortest_dist;
  oss.str( "" );
  oss << spiv_pair.site_label1_ << ":" << shortest_dist
      << ":" << spiv_pair.site_label2_;
  spiv_pair.label_ = oss.str();

  return spiv_pair;

}

// ***********************************************************************
void SpivMolecule::get_sites_atoms( const string &feature_name ,
				    vector<unsigned int> &atoms1 ) const {

  for( int i = 0 , is = pphore_site_atoms_.size() ; i < is ; ++i ) {
    if( feature_name == pphore_site_labels_[i] ) {
      atoms1 = pphore_site_atoms_[i];
      return;
    }
  }

}

// ***********************************************************************
void SpivMolecule::get_pairs_atoms( const string &feature_name ,
				    vector<unsigned int> &atoms1 ,
				    vector<unsigned int> &atoms2 ) const {

  for( int i = 0 , is = pphore_pairs_.size() ; i < is ; ++i ) {
    if( feature_name == pphore_pairs_[i].label_ ) {
      atoms1 = pphore_site_atoms_[pphore_pairs_[i].site1_];
      atoms2 = pphore_site_atoms_[pphore_pairs_[i].site2_];
      return;
    }
  }

}

// ***********************************************************************
void SpivMolecule::get_triplets_atoms( const string &feature_name ,
				       vector<unsigned int> &atoms1 ,
				       vector<unsigned int> &atoms2 ,
				       vector<unsigned int> &atoms3 ) const {
  
  for( int i = 0 , is = pphore_triplets_.size() ; i < is ; ++i ) {
    if( feature_name == pphore_triplets_[i].label_ ) {
      atoms1 = pphore_site_atoms_[pphore_triplets_[i].sites_[0]];
      atoms2 = pphore_site_atoms_[pphore_triplets_[i].sites_[1]];
      atoms3 = pphore_site_atoms_[pphore_triplets_[i].sites_[2]];
      return;
    }
  }

}


