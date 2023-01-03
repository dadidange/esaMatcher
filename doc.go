
// The esaMatcher module is a library to find matches between two (or more) strings using the enhanced suffix array (ESA). 
// This library is also a go wrapper for libdivsufsort and libsais (https://github.com/y-256/libdivsufsort and https://github.com/IlyaGrebnov/libsais) to construct the suffix array. 
// The wrapping functions are adapted from https://github.com/evolbioinf/esa. 
// The matching algorithms are adapted from Ohlebusch - Bioinformatics Algorithms (2013), phylonium (https://github.com/evolbioinf/phylonium) and my own masters thesis, the literate program par (https://github.com/dadidange/par_lp). 
//
// Introduction
//
// This package provides either the esaMatcher.Esa type or independent functions to work on byte slices directly.  
//
// To initialize the esa type, create a new instance of a esamatcher.Esa type, providing a slice of bytes that contain the text to create the ESA from. 
// Also, the library to use for creating the suffix array must be provided which is either "SaSais", 'SaDivSufSort' or 'SaNaive' for a naive suffix array construction.  
// For an ESA that also includes the reverse complement, a different construction method can be used.  
//	
//	e := esaMatcher.NewEsa(data, "sais")
//	//including reverse
//	eRev := esaMatcher.NewRevEsa(data, "sais")
//
// For details about the returned struct see the documentation below.
// Most parts of the documentation are adopted from the documentation in par_lp.
package esaMatcher