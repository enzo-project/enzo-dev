//----------------------------------------------------------------------
//
// Usage: jb-derive <operation> <new-metric> <argument1> [<argument2> ...]
//
// Description: Derive a new metric in terms of existing counters or
// Description: other previously-derived metrics
//
// Arguments:
//
//    Arguments depend on the operation; below we list the available operations
//    and their arguments.  We use the following notation for argument names:
//
//    "x" is the new derived counter
//    "y" and "z" are existing counters
//    "z[]" is an attribute, indexed by i
//    "a" is a (real) scalar
//    "i" in the descriptions refers to the index of the independent metric
//
//
// LOCAL UNARY OPERATORS
// ---------------------
// 
// Arguments:
//
//    Counter x
//    Real    a
//    Counter y
//
// Operators: 
//
//    scale      x(i) <-- a * y(i)
//    offset     x(i) <-- a + y(i)
//    invscale   x(i) <-- a / y(i)
//
//
// LOCAL BINARY OPERATORS
// ----------------------
//
// Arguments:
//
//    Counter x
//    Counter y
//    Counter z
//
// Operators: 
//
//    [add | +]   x(i) <-- y(i) + z(i)
//    [sub | -]   x(i) <-- y(i) - z(i)
//    [div | /]   x(i) <-- y(i) / z(i)
//    [mul | *]   x(i) <-- y(i) * z(i)
//
//
// REDUCTION OPERATORS
// -------------------
//
// Arguments:
//
//    Counter x
//    Counter y
//    Attribute z
//
// Operators: 
//
//    sum         x(i) <-- sum of y(i) over different z 
//    min         x(i) <-- minimum of y(i) over different z
//    max         x(i) <-- maximum of y(i) over different z 
//    avg         x(i) <-- average of y(i) over different z 
//   
//
// FINITE DIFFERENCE OPERATORS
// ---------------------------
//
// (For I and D operators, x(i) are sorted such that x(i-1) < x(i))
//
// Arguments:
//
//    Counter x
//    Counter y
//    Attribute z
//
// Operators: 
//
//    I           x(n) <-- y(0) + (z(1)-z(0))*y(1) + ... + (z(n)-z(n-1))*y(n)
//                
//    D (i = 0)   x(i) <-- y(0)
//    D (i > 0)   x(i) <-- (y(i) - y(i-1)) / (z(i) - z(i-1)) 
//
//
// INTERPOLATION OPERATOR
// ----------------------
//
// Arguments:
//
//    Counter x
//    Counter y
//    Attribute z
//
// Operators: 
//
//    interpolate (z(i) != "*") x(i) <-- z(i) 
//    interpolate (z(i) != "*") x(i) <-- z(i1) + (z(i2)-z(i1)) / (y(i2)-y(i1))
//             where i1 s.t. y(i1) < y(i) and z(k) < z(i1) ==> @@@@@
//             and   i2 = minarg_k y(k) > y(i)
//
//----------------------------------------------------------------------
//
// Copyright 2004 James Bordner
// Copyright 2004 Laboratory for Computational Astrophysics
// Copyright 2004 Regents of the University of California
//
//----------------------------------------------------------------------

#include "jb.h"
#include "sort.h"

main(int argc, char **argv)
{

  string operation;
  string new_metric;
  vector<string> arguments;
  map<string,string> Globals;
  vector<string> Attribute_names;
  vector<string> Counter_names;
  vector<string> Counter_types;
  string augregion;
  vector<string> counters;
  vector<string> augregion_split;
  map<string,vecstr> Records;

  vector<struct reg_struct *> aug_regions;

  char svalue[20];

  // Log start

  jb_log("start", argv[0]);

  // Parse command-line arguments

  operation  = argv[1];
  new_metric = argv[2];
  for (int i=3; i<argc; i++) {
    arguments.push_back (argv[i]);
  }

  bool valid_operation = false;

  FILE *fp = stdin;

  // -----------------------------------------------------------------------
  // Read and write globals
  // -----------------------------------------------------------------------

  jb_read_globals(fp, Globals);
  jb_print_globals (Globals);

  // -----------------------------------------------------------------------
  // Read, modify, and write data header
  // -----------------------------------------------------------------------

  jb_read_header(fp, Attribute_names, Counter_names, Counter_types);

  Counter_names.push_back (new_metric);
  Counter_types.push_back ("derived");
  
  jb_print_header (Attribute_names, Counter_names, Counter_types);

  // -----------------------------------------------------------------------
  // Local unary operators
  // -----------------------------------------------------------------------
  
  if (operation == "scale" ||
      operation == "offset" ||
      operation == "invscale") {

    valid_operation = true;

    // Get arguments

    double a;
    int iy;

    a = atof (arguments[0].c_str());
    for (iy=0; iy < Counter_names.size() && arguments[1] != Counter_names[iy]; iy++)
      ;
  
    // Check arguments

    if (iy >= Counter_names.size()) {
      fprintf (stderr, "%s Error: no Counter named %s!\n",
	       argv[0],arguments[1].c_str());
    }

    // Perform operation

    while (jb_read_next_record(fp,Attribute_names.size(),Counter_names.size()-1,
			       augregion, counters)) {

      char value[20];
      if (counters[iy] == "*") {
	sprintf (value,"*");
      } else {
	double y = atof(counters[iy].c_str());
	if (operation == "scale")    sprintf (value,"%g",a*y);
	if (operation == "offset")   sprintf (value,"%g",a+y);
	if (operation == "invscale") sprintf (value,"%g",a/y);
      }

      counters.push_back (value);

      jb_print_record (augregion, counters);
      
    }
  }

  // -----------------------------------------------------------------------
  // Local binary operators
  // -----------------------------------------------------------------------
  
  if (operation == "add" ||
      operation == "sub" ||
      operation == "div" ||
      operation == "mul" ||
      operation == "+" ||
      operation == "-" ||
      operation == "/" ||
      operation == "*") {

    valid_operation = true;

    // Get arguments

    int ix,iy;

    for (iy=0; 
	 iy < Counter_names.size() && arguments[0] != Counter_names[iy]; 
	 iy++)
      ;
    for (ix=0; 
	 ix < Counter_names.size() && arguments[1] != Counter_names[ix]; 
	 ix++)
      ;
  
    // Check arguments

    if (iy >= Counter_names.size()) {
      fprintf (stderr, "%s Error: no Counter named %s!\n",
	       argv[0],arguments[0].c_str());
    }
    if (ix >= Counter_names.size()) {
      fprintf (stderr, "%s Error: no Counter named %s!\n",
	       argv[0],arguments[1].c_str());
    }

    // Perform operation

    while (jb_read_next_record(fp,Attribute_names.size(),Counter_names.size()-1,
			       augregion, counters)) {

      if (counters[ix] == "*" || counters[iy] == "*") {
	sprintf (svalue,"*");
      } else {
	double x = atof(counters[ix].c_str());
	double y = atof(counters[iy].c_str());
	if (operation == "add" || operation == "+") sprintf (svalue,"%g",y+x);
	if (operation == "sub" || operation == "-") sprintf (svalue,"%g",y-x);
	if (operation == "mul" || operation == "*") sprintf (svalue,"%g",y*x);
	if (operation == "div" || operation == "/") sprintf (svalue,"%g",y/x);
      }

      counters.push_back (svalue);

      jb_print_record (augregion, counters);
      
    }
  }

// =======================================================================
// OPERATIONS THAT REQUIRE SORTING
// =======================================================================

  if (operation == "sum" ||
      operation == "min" ||
      operation == "max" ||
      operation == "avg" ||
      operation == "I" ||
      operation == "D" ||
      operation == "interpolate") {

    int ix,iy;

    for (iy=0; 
	 iy < Counter_names.size() && arguments[0] != Counter_names[iy]; 
	 iy++)
      ;
    for (ix=0; 
	 ix < Attribute_names.size() && arguments[1] != Attribute_names[ix]; 
	 ix++)
      ;
  
    // Check arguments

    if (iy >= Counter_names.size()) {
      fprintf (stderr, "%s Error: no Counter named %s!\n",
	       argv[0],arguments[0].c_str());
    }
    if (ix >= Attribute_names.size()) {
      fprintf (stderr, "%s Error: no Attribute named %s!\n",
	       argv[0],arguments[1].c_str());
    }

    // Adjust attribute index
    
    ++ix;

    // Generate list aug_regions of
    // [sorted-aug-regions, original aug-regions, indices, values-to-reduce]

    int index = 0;
    
    jb_read_records(fp,Attribute_names.size(),Counter_names.size()-1,Records);

    map<string,vecstr>::iterator it_augregion;
    for (it_augregion = Records.begin();
	 it_augregion != Records.end();
	 ++it_augregion) {

      double x,y;

      string augregion = it_augregion->first;

      if (Records[augregion][iy] == "*") {
	y = STAR_FLOAT;
      } else {
	y = atof (Records[augregion][iy].c_str());
      }

      string augregion_r;
      jb_augregion_rotate_(augregion,ix,augregion_r);
      
      jb_split (augregion,augregion_split);

      if (augregion_split[ix] == "*") {
	x = STAR_FLOAT;
      } else {
	x = atof (augregion_split[ix].c_str());
      }

      struct reg_struct *pregyx = new reg_struct;
      pregyx->sorted   = augregion_r;
      pregyx->original = augregion;
      pregyx->y = y;
      pregyx->x = x;

      aug_regions.push_back (pregyx);
      
      ++ index;
    }

    // Sort depending on operation type

    if (operation == "sum" ||
	operation == "min" ||
	operation == "max" ||
	operation == "avg") {
      jb_sort_sorted(aug_regions,0,aug_regions.size()-1);
    } else if (operation == "D" ||
	       operation == "I" ||
	       operation == "interpolate") {
      jb_sort_x(aug_regions,0,aug_regions.size()-1);
    }

    int istart, iend;
    string root_start,root_end;
    double value;
    istart = 0;

    // Get aug_region's root

    bool isFinished = false;
    
    while (! isFinished) {

        // Find ending index for reduction

      iend = istart;
      root_start = jb_augregion_root_(aug_regions[istart]->sorted);
      root_end   = root_start;
      while (! isFinished  && root_start == root_end) {
	++ iend;
	if (iend >= aug_regions.size()) {
	  isFinished = true;
	} else {
	  root_end = jb_augregion_root_(aug_regions[iend]->sorted);
	}
      }


      // ----------------------------------------------------------------------
      // Reduction operations
      // ----------------------------------------------------------------------

      if (operation == "sum" ||
	  operation == "min" ||
	  operation == "max" ||
	  operation == "avg") {
    
	value = aug_regions[istart]->y;

	double x,y;

	for (int i=istart+1; i<iend; i++) {
	  
	  y = aug_regions[i]->y;
	  if (value == STAR_FLOAT or y == STAR_FLOAT) {
            value = STAR_FLOAT;
	  } else if (operation == "sum") {                
	    value = value + y;
	  } else if (operation == "min" and value > y) {  
	    value = y;
	  } else if (operation == "max" and value < y) {  
	    value = y;
	  } else if (operation == "avg") {
	    value = value + y;
	  }
	}
                            
	if (operation == "avg" and value != STAR_FLOAT) {
	  value = value / (iend - istart);
	}

	// Store the values

	sprintf (svalue,"%g",value);
	for (int i=istart; i<iend; i++) {
	  augregion = aug_regions[i]->original;
	  Records[augregion].push_back(svalue);
	}
      }

// ----------------------------------------------------------------------
// Difference/quadrature operations
// ----------------------------------------------------------------------

      if (operation == "D" ||
	  operation == "I") {

	double x,y,xprev,yprev;

	y = aug_regions[istart]->y;
	x = aug_regions[istart]->x;
	yprev = y;
	xprev = x;
	value = y;

	// First value z(0) is y(0)

	augregion = aug_regions[istart]->original;
	sprintf (svalue,"%g",value);
	Records[augregion].push_back(svalue);
	for (int i=istart+1; i<iend; i++) {
	  y = aug_regions[i]->y;
	  x = aug_regions[i]->x;
	  // compute value
	  if (value == STAR_FLOAT || y == STAR_FLOAT || x == STAR_FLOAT) {
	    value = STAR_FLOAT;
	  } else {
	    double xv = x;
	    double xvprev = xprev;
	    if (operation == "I") { 
	      value = value + (xv - xvprev) * y;
	    } else if (operation == "D") { 
	      value = (y - yprev) / (xv - xvprev);
	    }
	  }

	  // store value
	  augregion = aug_regions[i]->original;
	  sprintf (svalue,"%g",value);
	  Records[augregion].push_back(svalue);
	  yprev = y;
	  xprev = x;
	}
      }

// ----------------------------------------------------------------------
// Interpolation operation
// ----------------------------------------------------------------------

      if (operation == "interpolate") {

	double x,y,x_first,y_first,x_second,y_second;
	int i_first,i_second,i_last;

	x = aug_regions[istart]->x;
	y = aug_regions[istart]->y;

	// Find i_first
	i_first = istart;
	while (i_first < iend) {
	  if (aug_regions[i_first]->y != STAR_FLOAT) break;
	  ++i_first;
	}

	x_first = aug_regions[i_first]->x;
	y_first = aug_regions[i_first]->y;

	// Find i_last
	i_last = iend - 1;
	while (i_last >= i_first) {
	  if (aug_regions[i_last]->y != STAR_FLOAT) break;
	  --i_last;
	}

	// If at most one y is not "*"...

	if (i_first == iend || i_first == i_last) {
	  // nothing to interpolate
	  if (i_first == iend) value = STAR_FLOAT;
	  // one value to "interpolate"--assume constant
	  if (i_first == i_last) value = aug_regions[i_first]->y;
	  for (int i=istart; i<iend; i++) {
	    augregion = aug_regions[i]->original;
	    sprintf (svalue,"%g",value);
	    Records[augregion].push_back(svalue);
	  }
	} else {

	  // If at least two y's are not "*"...

	  // Find i_second
	  i_second = i_first + 1;
	  while (i_second < iend) {
	    if (aug_regions[i_second]->y != STAR_FLOAT) break;
	    ++ i_second;
	  }

	  x_second = aug_regions[i_second]->x;
	  y_second = aug_regions[i_second]->y;

	  for (int i=istart; i<iend; i++) {
	    x = aug_regions[i]->x;
	    y = aug_regions[i]->y;
	    if (i < i_second && y == STAR_FLOAT) {
	      value = y_first + (x - x_first) 
		/     (x_second - x_first) * (y_second - y_first);
	    }	else if (i < i_second && y != STAR_FLOAT) {
	      value = y;
	    }	else {
	      // assert (i == i_second), so we know y != "*"
	      value = aug_regions[i]->y;
	      // also move to next (i_first, i_second) if needed
	      if (i_second != i_last) {
		i_first = i_second;
		x_first = aug_regions[i_first]->x;
		y_first = aug_regions[i_first]->y;
		  
		// Find i_second
		i_second = i_first+1;
		while (i_second < iend) {
		  if (aug_regions[i_second]->y != STAR_FLOAT) break;
		  ++i_second;
		}
		x_second = aug_regions[i_second]->x;
		y_second = aug_regions[i_second]->y;

	      }
	    }
	    augregion = aug_regions[i]->original;
	    sprintf (svalue,"%g",value);
	    Records[augregion].push_back(svalue);
	  }
	}
      }                 
      // Prepare for next reduction

      istart = iend;

    }
    jb_print_records (Records);
    
  }
  jb_log("stop", argv[0]);
}
