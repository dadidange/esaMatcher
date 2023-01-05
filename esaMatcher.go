package esaMatcher

/*
#cgo CFLAGS: -I/usr/local/include
#cgo LDFLAGS: -ldivsufsort64 -L/usr/local/include/ -lsais
#include <divsufsort64.h>
#include <libsais.h>
#include <stdlib.h>
*/
import "C"
import (
	"fmt"
	"log"
	"reflect"
	"unsafe"
	"bytes"
)

const (
	// Default SA Construction Method
	defaultSa = "SaSais"
)

// The Esa type holds relevant properties of the ESA. 
// All properties are initialized when the Esa is constructed and  
// can be acessed by calling their corresponding getters. 
type Esa struct {
	s          []byte
	sa         []int
	lcp        []int
	cld        []int
	strandSize int
}

// Return the suffix array of Esa.
func (e *Esa) Sa() []int        { return e.sa }
// Return the longest common prefix array of Esa.
func (e *Esa) Lcp() []int       { return e.lcp }
// Return the child array of Esa.
func (e *Esa) Cld() []int       { return e.cld }
// Return the sequence for the Esa.
func (e *Esa) Sequence() []byte { return e.s }
// Return the single strand size hold by the esa. 
// Equals len(Sequence) if the ESA was initialized w/o the reverse complement.
func (e *Esa) StrandSize() int  { return e.strandSize }

// Initialize new ESA with default values
func NewBaseEsa(s []byte) Esa{
	return NewEsa(s, defaultSa)
}

// Initialize new ESA of text t with given suffix array library.
// Does not include the reverse complement.
func NewEsa(s []byte, saLib string) Esa {
	strandSize := len(s)
	s = append(s, '$')

	sa := Sa(s, saLib)
	lcp := Lcp(s, sa)
	// Add last element to lcp if necessary
	if lcp[len(lcp)-1] != -1 {
		lcp = append(lcp, -1)
	}
	cld := Cld(lcp)

	return Esa{s, sa, lcp, cld, strandSize}
}

//Initialize a new ESA of text t that also includes the reverse complement.
func NewRevEsa(s []byte, saLib string) Esa {
	l := len(s)
	s = append(s, append([]byte{'#'}, RevComp(s)...)...)
	esa := NewEsa(s, saLib)
	esa.strandSize = l
	return esa
}

// Print the ESA to stdout. 
func (e *Esa) Print(numSeq int) {
	sa := e.Sa()
	lcp := e.Lcp()

	fmt.Print("i\t")
	fmt.Print("SA\t")
	fmt.Print("LCP\t")
	fmt.Print("CLD\t")
	fmt.Print("S[SA[i]..]\n")
	fmt.Print("-----------------------------------\n")
	m := len(sa)
	if numSeq == 0 || numSeq > len(lcp) {
		numSeq = len(lcp)
	}
	for i := 0; i < numSeq; i++ {
		if i >= m {
			fmt.Printf("%d\t", i)
			fmt.Printf("%s\t", "-")
			fmt.Printf("%d\t", lcp[i])
			fmt.Printf("%d\t", e.Cld()[i])
			fmt.Printf("%s\n", "-")
			continue
		}
		fmt.Printf("%d\t", i)
		fmt.Printf("%d\t", sa[i])
		fmt.Printf("%d\t", lcp[i])
		fmt.Printf("%d\t", e.Cld()[i])
		fmt.Printf("%s\n", e.Sequence()[sa[i]:])
	}
}

// Calculate the suffix array for a text t using the given method. 
// Options are empyt ("") for default, SaDivSufSort, SaSais and SaNaive. 
func Sa(t []byte, method string) []int {
	if method == ""{
		method = defaultSa
	}
	var sa []int
	if method == "SaDivSufSort" {
		sa = saDivSufSort(t)
	} else if method == "SaSais" {
		sa = saSais(t)
	}  else if method == "SaNaive" {
		sa = saNaive(t)
	} else {
		s := "Current options are:\n\t-SaDivSufSort\n\t-SaSais"
		log.Fatalf("library saLib = %s not defined to compute SA\n%s\n", method, s)
	}
	return sa
}

// Lcp returns the LCP-array of a given text t and the corresponding suffic array sa.
//
// The LCP-array contains the length of the common prefix of an element with its 
// predecessor in the alphabetically sorted suffix array.
func Lcp(t []byte, sa []int) []int {
	//from https://github.com/EvolBioInf/esa/
	n := len(t)
	lcp := make([]int, n)
	isa := make([]int, n)
	for i := 0; i < n; i++ {
		isa[sa[i]] = i
	}
	lcp[0] = -1
	l := 0
	for i := 0; i < n; i++ {
		j := isa[i]
		if j == 0 {
			continue
		}
		k := sa[j-1]
		for k+l < n && i+l < n && t[k+l] == t[i+l] {
			l++
		}
		lcp[j] = l
		l -= 1
		if l < 0 {
			l = 0
		}
	}
	return lcp
}

// Cld returns the child array from a LCP array.
//
// Concept
//
// A decrease in the LCP-array indicates a local minimum, that is, a node in our tree. 
// Therefore, a right child, left child or both are added at these positions.
// Minima inside the lcp array indicate boundaries between lcp intervals and 
// hence different paths through our suffix tree. 
// We use these boundaries to navigate through the tree structure 
// which can be described as “guided binary search” (Frith and Shrestha, 2018). 
//
// The left and right child pointers CLD.L and CLD.R , respectively,
// can be merged together to reduce memory requirements.
func Cld(lcp []int) []int {
	// initialize stack
	stack := []int{}
	top := func() int {
		return stack[len(stack)-1]
	}
	pop := func() int {
		t := top()
		stack = stack[:len(stack)-1]
		return t
	}
	push := func(i int) {
		stack = append(stack, i)
	}

	n := len(lcp) - 1
	cld := make([]int, n+1)
	cld[0] = n
	push(0)
	var last int

	for k := 1; k <= n; k++ {
		for lcp[k] < lcp[top()] {
			last = pop()
			for lcp[top()] == lcp[last] {
				cld[top()] = last // CLD[k].R = CLD[k]
				last = pop()
			}
			if lcp[k] < lcp[top()] {
				cld[top()] = last // CLD[k].R = CLD[k]
			} else {
				cld[k-1] = last // CLD[k].L = CLD[k-1].R = CLD[k-1]
			}
		}
		push(k)
	}
	return cld
}

// RevComp returns the reverse complement of a given DNA string. 
//
// The DNA nucleotide characters 'A','C','G' and 'T' will be translated to their counterparts 'T','G','C' and 'A'. 
// Any other character will be converted to a 'N'.
func RevComp(seq []byte) []byte {
	n := len(seq)
	revSeq := make([]byte, n)

	//Reverse
	//for i, j := 0, n-1; i < j; i, j = i+1, j-1 {
	//revSeq[i], revSeq[j] = seq[j], seq[i]
	//}
	//slower but more secure
	for i:= 0; i<n; i++ {
		revSeq[(n-i)-1] = seq[i]
	}

	//Complement
	f := func (r rune) rune  {
		switch{
		case r == 'A':
			return 'T'
		case r == 'T':
			return 'A'
		case r == 'G':
			return 'C'
		case r == 'C':
			return 'G'
		default:
			return 'N'
		}
	}
	return  bytes.Map(f, revSeq)
}

// This function is obsolete and will be removed in later versions
func RevCompObs(s []byte) []byte {
	n := len(s)

	f := []byte("ACGTN")
	r := []byte("TGCAN")

	rev := make([]byte, len(s))
	var dic [256]byte
	for i := range dic {
		dic[i] = byte(i)
	}
	for i, v := range f {
		dic[v] = r[i]
	}
	for i, v := range rev {
		rev[i] = dic[v]
	}

	//Reverse
	for i, j := 0, n-1; i < j; i, j = i+1, j-1 {
		rev[i], rev[j] = dic[s[j]], dic[s[i]]
	}
	return rev
}

// Wrapper for the C-Library LibDivSufSort, adopded from https://github.com/EvolBioInf/esa/. 
// This function takes a text t and returns its suffix array SA. 
func saDivSufSort(t []byte) []int {
	//from https://github.com/EvolBioInf/esa/
	var sa []int
	header := (*reflect.SliceHeader)(unsafe.Pointer(&t))
	ct := (*C.sauchar_t)(unsafe.Pointer(header.Data))
	n := len(t)
	csa := (*C.saidx64_t)(C.malloc(C.size_t(n * C.sizeof_saidx64_t)))
	cn := C.saidx64_t(n)
	err := int(C.divsufsort64(ct, csa, cn))
	if err != 0 {
		log.Fatalf("divsufsort failed with code %d\n", err)
	}
	header = (*reflect.SliceHeader)((unsafe.Pointer(&sa)))
	header.Cap = n
	header.Len = n
	header.Data = uintptr(unsafe.Pointer(csa))
	return sa
}

// Wrapper for the C-Library Libsais.  
// This function takes a text t and returns its suffix array SA. 
func saSais(t []byte) []int {
	// var sa []int
	n := len(t)
	sh := (*reflect.SliceHeader)(unsafe.Pointer(&t))
	ct := (*C.uint8_t)(unsafe.Pointer(sh.Data))
	csa := (*C.int32_t)(C.malloc(C.size_t(n * C.sizeof_int32_t)))
	cn := C.int32_t(n)
	cz := C.int32_t(0)
	czp := (*C.int32_t)(nil)
	err := int(C.libsais(ct, csa, cn, cz, czp))
	if err != 0 {
		log.Fatalf("divsufsort failed with code %d\n", err)
	}

	var sa []int32

	sh = (*reflect.SliceHeader)((unsafe.Pointer(&sa)))
	sh.Cap = n
	sh.Len = n
	sh.Data = uintptr(unsafe.Pointer(csa))

	sa2 := make([]int, len(sa))

	for i, v := range sa {
		sa2[i] = int(v)
	}

	return sa2
}

// The type EsaInterval represents an interval inside our ESA. 
//
// It contains an index for its starting and ending position. 
// Also, it has its middle mid which is the end of its first child
// interval and the length l defined. 
// The length l represents the length of the lcp at mid during GetInterval and GetMatch. 
// When the actual match is returned we will write the length of the match to l.
type EsaInterval struct {
	start int
	end   int
	mid   int
	l     int
}

// Starting index of the interval in the ESA. 
func (i *EsaInterval)Start() int{return i.start}
// Ending index of the interval in the ESA. 
func (i *EsaInterval)End() int{return i.end}
// middle index which is the end of its first child.
func (i *EsaInterval)Mid() int{return i.mid}
// length of the lcp at mid or match length
func (i *EsaInterval)L() int{return i.l}

// Initialise new EsaInterval using the starting and ending indices and the esa on which the interval lies. 
//
// Concept
//
// For most cases we can find the values for l and mid with the help of the child array.
// Let’s take a look again at the left pointer CLD.L. CLD.L points to the first local minimum
// of its “left” (above) interval. If an interval i ends at position j, CLD[J+1].L points to the
// end of the first child interval of i normally. We only have to be careful for singletons
// and the last child interval. If the given interval {i,j} is the last child interval of some
// parent interval, CLD.L[j+1] might point to an minimum that starts before i since it
// always points to the first minimum of the interval {h,j} that ends there with h < i ≤ j.
// For those two intervals that end at j, CLD[j+1].L points to the minimum of the larger
// interval, that is {h,j}. Hence we can not take this pointer to find the minimum of {i,j}
// In this case we can make use of CLD[h].R that points to the next minimum of {h,j}.
// Since {i,j} must be a child of {h,j} we can follow the right pointers until we will 
// eventually arrive at the start index of {i,j}.
func NewEsaInterval(start, end int, e Esa) EsaInterval {
	//Check for empty, invalid or singleton interval
	if start >= end {
		//singleton
		if start >= 0 {
			return EsaInterval{start, end, start, e.lcp[end]}
		} else {
			//empty or invalid
			return EmptyEsaInterval()
		}
	}

	m := e.cld[end] //CLD.L(m+1) = cld(m)
	for m <= start {
		m = e.cld[m]
	}
	return EsaInterval{start, end, m, e.lcp[m]}
}

// Returns a new, empty interval.
func EmptyEsaInterval() EsaInterval {
	return EsaInterval{-1, -1, -1, -1}
}

// Given an interval i on the ESA and a character c GetInterval returns the subinterval
// of i that starts with c.
//
// Implementation
//
// GetInterval starts by looking for singletons, i.e. intervals that only contain one element.
// If i is not a singleton interval we want to loop through the subintervals one level below
// the given interval, the child intervals. We initialize our interval by setting the
// upper and lower bounds. The first child interval starts where i starts and ends mid, the
// first local minimum.
// We loop through the child intervals and check if any interval starts with c. 
func (e *Esa)GetInterval(i EsaInterval, c byte) (EsaInterval){
	// Check Singleton Interval
	if i.start == i.end{
	  if(e.s[e.sa[i.start]] == c){
		return i
	  } else {
		//Return empty interval
		return EmptyEsaInterval()
	  }
	}
	lower := i.start
	upper := i.mid
	l := i.l
	for e.lcp[upper] == l {
	  if (e.s[e.sa[lower]+l] == c){
		//match found
		return NewEsaInterval(lower, upper-1, *e)
	  }
	  //increment interval boundaries
	  lower = upper
	  //check for singleton
	  if (lower == i.end){
		break
	  }
	  upper = e.cld[upper] //CLD.R(m) = cld(m)
	}
	if (e.s[e.sa[lower] + l] == c){
	  return NewEsaInterval(lower, i.end, *e)
	} else {
	  return EmptyEsaInterval()
	}
}
   
// GetMatch returns the longest prefix of the query that matches the ESA, 
// that is any suffix of the reference. 
//
// Implementation
//
// GetMatch calls GetInterval once per character at most. For every character in the
// query we can call GetInterval with the child interval returned by the previous character.
func (e *Esa)GetMatch(query []byte) EsaInterval{
	in := NewEsaInterval(0, len(e.s)-1, *e)
	cld := EmptyEsaInterval()
	k := 0
	m := len(query)
	for k < m{
		cld = e.GetInterval(in, query[k])
		//check if empty interval was returned -> no match this round
		if (cld.start == -1 && cld.end == -1){
			// if k > 0, there where matches in previous rounds
			if (k == 0){
				return cld
			}
			in.l = k
			return in
		}

		k++
		//extend the match and skip common prefixes
		in = cld
		l := in.l
		if(in.start == in.end || l > m){
			l = m
		}		
		for saIdx:=e.sa[in.start]; k < l; k++ {
			if(e.s[saIdx+k] != query[k]){
				in.l = k
				return in
			}
		}
	}
	in.l = m
	return in
}