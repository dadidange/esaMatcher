# esaMatcher
Library to effectively find matches between two (or more) strings using the enhanced suffix array (ESA). 
This library is also a go wrapper for [libdivsufsort](https://github.com/y-256/libdivsufsort) and (https://github.com/IlyaGrebnov/libsais) since both are used to construct the suffix array structure. 
The wrapping functions are adapted from [esa](https://github.com/evolbioinf/esa). 
The matching algorithms are adapted from [Ohlebusch (2013)](https://www.uni-ulm.de/en/in/theo/m/ohlebusch/book-bioinformatics-algorithms/), [phylonium](https://github.com/evolbioinf/phylonium) and my own masters thesis. 

## Prerequisites
To run this library both, *libdivsufsort* and *libsais* should be installed.  

## Benchmarks
Since both, *libsivsufsort* and *libsais* are optimized for certain hardware please benchmark for better results. 
For comparison, go's own [suffix-array library](https://pkg.go.dev/index/suffixarray) is tested as well but not included in the package yet. 

The memory benchmark for DivSufSort seems to not work so far.

To benchmark run
```
$ go test -bench=. -benchtime=20s -benchmem
```
returned for me
```
goos: linux
goarch: amd64
pkg: github.com/dadidange/esaMatcher
cpu: AMD Ryzen 5 5600X 6-Core Processor             
Benchmark_Sa_LibSais_5MBP-12                 176         133370608 ns/op        40001538 B/op          1 allocs/op
Benchmark_Sa_LibDivSufSort_5MBP-12           100         202835960 ns/op               0 B/op          0 allocs/op
Benchmark_Sa_GoSa_5MP-12                     123         193818537 ns/op        20004944 B/op          2 allocs/op
Benchmark_Sa_LibSais_50MBP-12                 14        1579467532 ns/op        400007168 B/op         1 allocs/op
Benchmark_Sa_LibDivSufSort_50MBP-12            8        2734851842 ns/op               0 B/op          0 allocs/op
Benchmark_Sa_GoSa_50MP-12                      7        3328684829 ns/op        200007760 B/op         2 allocs/op
Benchmark_Sa_LibSais_Eco-12                  162         148501200 ns/op        44228609 B/op          1 allocs/op
Benchmark_Sa_LibDivSufSort_Eco-12            100         223299174 ns/op               0 B/op          0 allocs/op
Benchmark_Sa_GoSa_Eco-12                     100         225492325 ns/op        22118483 B/op          2 allocs/op
Benchmark_Esa_LibSais_Eco-12                  69         339440262 ns/op        232318246 B/op        21 allocs/op
Benchmark_Esa_LibDivSufSort_Eco-12            52         438802555 ns/op        188089607 B/op        20 allocs/op
PASS
ok      github.com/dadidange/esaMatcher 327.585s
``` 

