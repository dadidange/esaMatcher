# esaMatcher
Library to effectively find matches between two (or more) strings using the enhanced suffix array (ESA). 
This library is also a go wrapper for [libdivsufsort](https://github.com/y-256/libdivsufsort) and (https://github.com/IlyaGrebnov/libsais) since both are used to construct the suffix array structure. 
The wrapping functions are adapted from [esa](https://github.com/evolbioinf/esa). 
The matching algorithms are adapted from [Ohlebusch (2013)](https://www.uni-ulm.de/en/in/theo/m/ohlebusch/book-bioinformatics-algorithms/), [phylonium](https://github.com/evolbioinf/phylonium) and my own masters thesis. 

## Prerequisites
To run this library both, *libdivsufsort* and *libsais* should be installed.  

## Benchmarks
Since both, *libsivsufsort* and **libsais* are optimized for certain hardware please benchmark for better results. 

To benchmark run
```
$ go test -bench=. -benchtime=20s
```
returned for me
```
goos: linux
goarch: amd64
pkg: github.com/dadidange/esaMatcher
cpu: AMD Ryzen 5 5600X 6-Core Processor             
BenchmarkLibSais5MBP-12           	     183	 130556902 ns/op
BenchmarkLibDivSufSort5MBP-12     	     100	 200920172 ns/op
BenchmarkLibSais50MBP-12          	      14	1598579853 ns/op
BenchmarkLibDivSufSort50MBP-12    	       8	2706890674 ns/op
BenchmarkSaisEco-12               	     160	 148516043 ns/op
BenchmarkLibDivSufSortEco-12      	     100	 224297015 ns/op
BenchmarkEsaSaisEco-12            	      63	 349417870 ns/op
BenchmarkEsaDivSufSortEco-12      	      54	 424569795 ns/op
PASS
ok  	github.com/dadidange/esaMatcher	243.965s
``` 