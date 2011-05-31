// comparer functions for std::sort

struct cmp_grid {
  bool operator()(particle_data const& a, particle_data const& b) const {
    if (a.grid < b.grid) return true;
    else return false;
  }
};

struct cmp_proc {
  bool operator()(particle_data const& a, particle_data const& b) const {
    if (a.proc < b.proc) return true;
    else return false;
  }
};

struct cmp_star_grid {
  bool operator()(star_data const& a, star_data const& b) const {
    if (a.grid < b.grid) return true;
    else return false;
  }
};

struct cmp_star_proc {
  bool operator()(star_data const& a, star_data const& b) const {
    if (a.proc < b.proc) return true;
    else return false;
  }
};

struct cmp_hkey {
  bool operator()(hilbert_data const& a, hilbert_data const& b) const {
    if (a.hkey < b.hkey) return true;
    else return false;
  }
};

struct cmp_grid1 {
  bool operator()(two_int const& a, two_int const& b) const {
    if (a.grid < b.grid) return true;
    else return false;
  }
};

struct cmp_proc1 {
  bool operator()(two_int const& a, two_int const& b) const {
    if (a.proc < b.proc) return true;
    else return false;
  }
};

#ifdef TRANSFER
#include "PhotonPackage.h"
struct cmp_ss {
  bool operator()(PhotonPackageEntry const& a, 
		  PhotonPackageEntry const& b) const {
    if (a.CurrentSource < b.CurrentSource)
      return true;
    else if (a.CurrentSource > b.CurrentSource)
      return false;
    else {
      if ( a.level < b.level)
	return true;
      else if ( a.level > b.level)
	return false;
      else {
	if ( a.ipix < b.ipix)
	  return true;
	else if ( a.ipix > b.ipix)
	  return false;
	else {
	  if ( a.Type < b.Type)
	    return true;
	  else if ( a.Type > b.Type)
	    return false;
	}
      }
    } // ENDELSE (top)
    return true;
  } // END bool operator()
};
#endif /* TRANSFER */
