

# esaMatcher
`import "github.com/dadidange/esaMatcher"`

* [Overview](#pkg-overview)
* [Index](#pkg-index)

## <a name="pkg-overview">Overview</a>
The esaMatcher module is a library to find matches between two (or more) strings using the enhanced suffix array (ESA).
This library is also a go wrapper for libdivsufsort and libsais (<a href="https://github.com/y-256/libdivsufsort">https://github.com/y-256/libdivsufsort</a> and <a href="https://github.com/IlyaGrebnov/libsais">https://github.com/IlyaGrebnov/libsais</a>) to construct the suffix array.
The wrapping functions are adapted from <a href="https://github.com/evolbioinf/esa">https://github.com/evolbioinf/esa</a>.
The matching algorithms are adapted from Ohlebusch - Bioinformatics Algorithms (2013), phylonium (<a href="https://github.com/evolbioinf/phylonium">https://github.com/evolbioinf/phylonium</a>) and my own masters thesis, the literate program par (<a href="https://github.com/dadidange/par_lp">https://github.com/dadidange/par_lp</a>).

### Prerequisites
To run this library both, *libdivsufsort* and *libsais* should be installed, also GO.

### Benchmarks
Since both, *libsivsufsort* and *libsais* are optimized for certain hardware they should be benchmarked for better results.
For comparison, go's own [suffix-array library](<a href="https://pkg.go.dev/index/suffixarray">https://pkg.go.dev/index/suffixarray</a>) is tested as well.
But it is currently not included for building the SA.

To benchmark run


	$ go test -bench=. -benchtime=20s

this returned for me

goos: linux
goarch: amd64
pkg: github.com/dadidange/esaMatcher


	cpu: AMD Ryzen 5 5600X 6-Core Processor
	Benchmark_Sa_LibSais_5MBP-12                 176         133370608 ns/op
	Benchmark_Sa_LibDivSufSort_5MBP-12           100         202835960 ns/op
	Benchmark_Sa_GoSa_5MP-12                     123         193818537 ns/op
	Benchmark_Sa_LibSais_50MBP-12                 14        1579467532 ns/op
	Benchmark_Sa_LibDivSufSort_50MBP-12            8        2734851842 ns/op
	Benchmark_Sa_GoSa_50MP-12                      7        3328684829 ns/op
	Benchmark_Sa_LibSais_Eco-12                  162         148501200 ns/op
	Benchmark_Sa_LibDivSufSort_Eco-12            100         223299174 ns/op
	Benchmark_Sa_GoSa_Eco-12                     100         225492325 ns/op
	Benchmark_Esa_LibSais_Eco-12                  69         339440262 ns/op
	Benchmark_Esa_LibDivSufSort_Eco-12            52         438802555 ns/op
	PASS
	ok      github.com/dadidange/esaMatcher 327.585s

### Quickstart
This package provides either the esaMatcher.Esa type or independent functions to work on byte slices directly.

To initialize the esa type, create a new instance of a esamatcher.Esa type, providing a slice of bytes that contain the text to create the ESA from.
Also, the library to use for creating the suffix array must be provided which is either "SaSais", 'SaDivSufSort' or 'SaNaive' for a naive suffix array construction.
For an ESA that also includes the reverse complement, a different construction method can be used.


	e := esaMatcher.NewEsa(data, "sais")
	//including reverse
	eRev := esaMatcher.NewRevEsa(data, "sais")

For details about the returned struct see the documentation below.
Most parts of the documentation are adopted from the documentation in par_lp.




## <a name="pkg-index">Index</a>
* [func Cld(lcp []int) []int](#Cld)
* [func Lcp(t []byte, sa []int) []int](#Lcp)
* [func RevComp(seq []byte) []byte](#RevComp)
* [func RevCompObs(s []byte) []byte](#RevCompObs)
* [func Sa(t []byte, method string) []int](#Sa)
* [type Esa](#Esa)
  * [func NewEsa(s []byte, saLib string) Esa](#NewEsa)
  * [func NewRevEsa(s []byte, saLib string) Esa](#NewRevEsa)
  * [func (e *Esa) Cld() []int](#Esa.Cld)
  * [func (e *Esa) GetInterval(i EsaInterval, c byte) EsaInterval](#Esa.GetInterval)
  * [func (e *Esa) GetMatch(query []byte) EsaInterval](#Esa.GetMatch)
  * [func (e *Esa) Lcp() []int](#Esa.Lcp)
  * [func (e *Esa) Print(numSeq int)](#Esa.Print)
  * [func (e *Esa) Sa() []int](#Esa.Sa)
  * [func (e *Esa) Sequence() []byte](#Esa.Sequence)
  * [func (e *Esa) StrandSize() int](#Esa.StrandSize)
* [type EsaInterval](#EsaInterval)
  * [func EmptyEsaInterval() EsaInterval](#EmptyEsaInterval)
  * [func NewEsaInterval(start, end int, e Esa) EsaInterval](#NewEsaInterval)
  * [func (i *EsaInterval) End() int](#EsaInterval.End)
  * [func (i *EsaInterval) L() int](#EsaInterval.L)
  * [func (i *EsaInterval) Mid() int](#EsaInterval.Mid)
  * [func (i *EsaInterval) Start() int](#EsaInterval.Start)


#### <a name="pkg-files">Package files</a>
[doc.go](/src/github.com/dadidange/esaMatcher/doc.go) [esaMatcher.go](/src/github.com/dadidange/esaMatcher/esaMatcher.go) [naive.go](/src/github.com/dadidange/esaMatcher/naive.go) 





## <a name="Cld">func</a> [Cld](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=170, int=4262, int=4287))
``` go
func Cld(lcp []int) []int
```
Cld returns the child array from a LCP array.

### Concept
A decrease in the LCP-array indicates a local minimum, that is, a node in our tree.
Therefore, a right child, left child or both are added at these positions.
Minima inside the lcp array indicate boundaries between lcp intervals and
hence different paths through our suffix tree.
We use these boundaries to navigate through the tree structure
which can be described as “guided binary search” (Frith and Shrestha, 2018).

The left and right child pointers CLD.L and CLD.R , respectively,
can be merged together to reduce memory requirements.



## <a name="Lcp">func</a> [Lcp](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=129, int=3205, int=3239))
``` go
func Lcp(t []byte, sa []int) []int
```
Lcp returns the LCP-array of a given text t and the corresponding suffic array sa.

The LCP-array contains the length of the common prefix of an element with its
predecessor in the alphabetically sorted suffix array.



## <a name="RevComp">func</a> [RevComp](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=213, int=5176, int=5207))
``` go
func RevComp(seq []byte) []byte
```
RevComp returns the reverse complement of a given DNA string.

The DNA nucleotide characters 'A','C','G' and 'T' will be translated to their counterparts 'T','G','C' and 'A'.
Any other character will be converted to a 'N'.



## <a name="RevCompObs">func</a> [RevCompObs](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=245, int=5746, int=5778))
``` go
func RevCompObs(s []byte) []byte
```
This function is obsolete and will be removed in later versions



## <a name="Sa">func</a> [Sa](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=107, int=2561, int=2599))
``` go
func Sa(t []byte, method string) []int
```
Calculate the suffix array for a text t using the given method.
Options are empyt ("") for default, SaDivSufSort, SaSais and SaNaive.




## <a name="Esa">type</a> [Esa](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=27, int=502, int=610))
``` go
type Esa struct {
    // contains filtered or unexported fields
}

```
The Esa type holds relevant properties of the ESA.
All properties are initialized when the Esa is constructed and
can be acessed by calling their corresponding getters.







### <a name="NewEsa">func</a> [NewEsa](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=48, int=1209, int=1248))
``` go
func NewEsa(s []byte, saLib string) Esa
```
initialize new ESA of text t without the reverse complement


### <a name="NewRevEsa">func</a> [NewRevEsa](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=64, int=1553, int=1595))
``` go
func NewRevEsa(s []byte, saLib string) Esa
```
make new ESA that includes of text t the reverse complement





### <a name="Esa.Cld">func</a> (\*Esa) [Cld](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=40, int=828, int=853))
``` go
func (e *Esa) Cld() []int
```
Return the child array of Esa.




### <a name="Esa.GetInterval">func</a> (\*Esa) [GetInterval](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=401, int=10886, int=10947))
``` go
func (e *Esa) GetInterval(i EsaInterval, c byte) EsaInterval
```
Given an interval i on the ESA and a character c GetInterval returns the subinterval
of i that starts with c.

### Implementation
GetInterval starts by looking for singletons, i.e. intervals that only contain one element.
If i is not a singleton interval we want to loop through the subintervals one level below
the given interval, the child intervals. We initialize our interval by setting the
upper and lower bounds. The first child interval starts where i starts and ends mid, the
first local minimum.
We loop through the child intervals and check if any interval starts with c.




### <a name="Esa.GetMatch">func</a> (\*Esa) [GetMatch](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=441, int=11887, int=11934))
``` go
func (e *Esa) GetMatch(query []byte) EsaInterval
```
GetMatch returns the longest prefix of the query that matches the ESA,
that is any suffix of the reference.

### Implementation
GetMatch calls GetInterval once per character at most. For every character in the
query we can call GetInterval with the child interval returned by the previous character.




### <a name="Esa.Lcp">func</a> (\*Esa) [Lcp](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=38, int=745, int=770))
``` go
func (e *Esa) Lcp() []int
```
Return the longest common prefix array of Esa.




### <a name="Esa.Print">func</a> (\*Esa) [Print](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=73, int=1754, int=1785))
``` go
func (e *Esa) Print(numSeq int)
```
Print the ESA to stdout.




### <a name="Esa.Sa">func</a> (\*Esa) [Sa](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=36, int=647, int=671))
``` go
func (e *Esa) Sa() []int
```
Return the suffix array of Esa.




### <a name="Esa.Sequence">func</a> (\*Esa) [Sequence](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=42, int=913, int=944))
``` go
func (e *Esa) Sequence() []byte
```
Return the sequence for the Esa.




### <a name="Esa.StrandSize">func</a> (\*Esa) [StrandSize](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=45, int=1090, int=1120))
``` go
func (e *Esa) StrandSize() int
```
Return the single strand size hold by the esa.
Equals len(Sequence) if the ESA was initialized w/o the reverse complement.




## <a name="EsaInterval">type</a> [EsaInterval](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=331, int=7996, int=8067))
``` go
type EsaInterval struct {
    // contains filtered or unexported fields
}

```
The type EsaInterval represents an interval inside our ESA.

It contains an index for its starting and ending position.
Also, it has its middle mid which is the end of its first child
interval and the length l defined.
The length l represents the length of the lcp at mid during GetInterval and GetMatch.
When the actual match is returned we will write the length of the match to l.







### <a name="EmptyEsaInterval">func</a> [EmptyEsaInterval](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=386, int=10198, int=10233))
``` go
func EmptyEsaInterval() EsaInterval
```
Returns a new, empty inerval.


### <a name="NewEsaInterval">func</a> [NewEsaInterval](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=366, int=9761, int=9815))
``` go
func NewEsaInterval(start, end int, e Esa) EsaInterval
```
Initialise new EsaInterval using the starting and ending indeces and the esa on which the interval lies.

### Concept
For most cases we can find the values for l and mid with the help of the child array.
Let’s take a look again at the left pointer CLD.L. CLD.L points to the first local minimum
of its “left” (above) interval. If an interval i ends at position j, CLD[J+1].L points to the
end of the first child interval of i normally. We only have to be careful for singletons
and the last child interval. If the given interval {i,j} is the last child interval of some
parent interval, CLD.L[j+1] might point to an minimum that starts before i since it
always points to the first minimum of the interval {h,j} that ends there with h < i ≤ j.
For those two intervals that end at j, CLD[j+1].L points to the minimum of the larger
interval, that is {h,j}. Hence we can not take this pointer to find the minimum of {i,j}
In this case we can make use of CLD[h].R that points to the next minimum of {h,j}.
Since {i,j} must be a child of {h,j} we can follow the right pointers until we will
eventually arrive at the start index of {i,j}.

Both, the left and the right pointer are stored in a single list
with CLD[i].L = CLD[i-1] and CLD[i].R = CLD[i] to save space.





### <a name="EsaInterval.End">func</a> (\*EsaInterval) [End](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=341, int=8210, int=8240))
``` go
func (i *EsaInterval) End() int
```
Ending index of the interval in the ESA.




### <a name="EsaInterval.L">func</a> (\*EsaInterval) [L](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=345, int=8397, int=8425))
``` go
func (i *EsaInterval) L() int
```
length of the lcp at mid or match length




### <a name="EsaInterval.Mid">func</a> (\*EsaInterval) [Mid](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=343, int=8308, int=8338))
``` go
func (i *EsaInterval) Mid() int
```
middle index which is the end of its first child.




### <a name="EsaInterval.Start">func</a> (\*EsaInterval) [Start](https://github.com/dadidange/esaMatcher/blob/main/esaMatcher.go%!(EXTRA string=/target/esaMatcher.go, int=339, int=8116, int=8148))
``` go
func (i *EsaInterval) Start() int
```
Starting index of the interval in the ESA.








- - -
Generated by [godoc2md](http://godoc.org/github.com/davecheney/godoc2md)
