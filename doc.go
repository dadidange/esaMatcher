
// The esaMatcher module is a library to find matches between two (or more) strings using the enhanced suffix array (ESA). 
// This library is also a go wrapper for libdivsufsort and libsais (https://github.com/y-256/libdivsufsort and https://github.com/IlyaGrebnov/libsais) to construct the suffix array. 
// The wrapping functions are adapted from https://github.com/evolbioinf/esa. 
// The matching algorithms are adapted from Ohlebusch - Bioinformatics Algorithms (2013), phylonium (https://github.com/evolbioinf/phylonium) and my own masters thesis, the literate program par (https://github.com/dadidange/par_lp). 
//
// Prerequisites
//
// To run this library both, *libdivsufsort* and *libsais* should be installed, also GO. 
//
// Benchmarks
//
// Since both, *libsivsufsort* and *libsais* are optimized for certain hardware they should be benchmarked for better results. 
// For comparison, go's own [suffix-array library](https://pkg.go.dev/index/suffixarray) is tested as well.
// But it is currently not included for building the SA. 
//
// To benchmark run
// 
// 	$ go test -bench=. -benchtime=20s 
// 
// this returned for me
// 
// 	goos: linux
// 	goarch: amd64
// 	pkg: github.com/dadidange/esaMatcher
// 	cpu: AMD Ryzen 5 5600X 6-Core Processor             
// 	Benchmark_Sa_LibSais_5MBP-12                 176         133370608 ns/op
// 	Benchmark_Sa_LibDivSufSort_5MBP-12           100         202835960 ns/op
// 	Benchmark_Sa_GoSa_5MP-12                     123         193818537 ns/op
// 	Benchmark_Sa_LibSais_50MBP-12                 14        1579467532 ns/op
// 	Benchmark_Sa_LibDivSufSort_50MBP-12            8        2734851842 ns/op
// 	Benchmark_Sa_GoSa_50MP-12                      7        3328684829 ns/op
// 	Benchmark_Sa_LibSais_Eco-12                  162         148501200 ns/op
// 	Benchmark_Sa_LibDivSufSort_Eco-12            100         223299174 ns/op
// 	Benchmark_Sa_GoSa_Eco-12                     100         225492325 ns/op
// 	Benchmark_Esa_LibSais_Eco-12                  69         339440262 ns/op
// 	Benchmark_Esa_LibDivSufSort_Eco-12            52         438802555 ns/op
// 	PASS
// 	ok      github.com/dadidange/esaMatcher 327.585s
// 
//
// Quickstart
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