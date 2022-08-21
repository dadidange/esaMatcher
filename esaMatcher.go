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
)

const (
	// Default SA Construction Method
	defaultSa = "SaSais"
)

type Esa struct {
	s          []byte
	sa         []int
	lcp        []int
	cld        []int
	strandSize int
}

func (e *Esa) Sa() []int        { return e.sa }
func (e *Esa) Lcp() []int       { return e.lcp }
func (e *Esa) Cld() []int       { return e.cld }
func (e *Esa) Sequence() []byte { return e.s }
func (e *Esa) StrandSize() int  { return e.strandSize }

//make new ESA that includes the reverse complement
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

//make new ESA that includes the reverse complement
func NewRevEsa(s []byte, saLib string) Esa {
	l := len(s)
	s = append(s, append([]byte{'#'}, RevComp(s)...)...)
	esa := NewEsa(s, saLib)
	esa.strandSize = l
	return esa
}

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

func Cld(lcp []int) []int {
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

func RevComp(s []byte) []byte {
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

type EsaInterval struct {
	start int
	end   int
	mid   int
	l     int
}

func (i *EsaInterval)Start() int{return i.start}
func (i *EsaInterval)End() int{return i.end}
func (i *EsaInterval)Mid() int{return i.mid}
func (i *EsaInterval)L() int{return i.l}

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

func EmptyEsaInterval() EsaInterval {
	return EsaInterval{-1, -1, -1, -1}
}


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
  
  
  func (e *Esa)GetMatch(query []byte) EsaInterval{
	in := NewEsaInterval(0, len(e.s)-1, *e)
	cld := EmptyEsaInterval()
	k := 0
	m := len(query)
	for k < m{
	  cld = e.GetInterval(in, query[k])
	  if (cld.start == -1 && cld.end == -1){
		if (k == 0){
		  return cld
		}
		in.l = k
		return in
	  }
  
	  k++ //the k-th character was matched in
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