package esaMatcher

/*
#cgo LDFLAGS: -ldivsufsort64 
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

func NewEsa(s []byte, saLib string) Esa {
	strandSize := len(s)
	s = append(s, '$')

	var sa []int
	if saLib == "SaDivSufSort" {
		sa = SaDivSufSort(s)
	} else if saLib == "SaLibSais" {

	} else {
		s := "Current options are:\n\t-SaDivSufSort\n\t-SaLibSais"
		log.Fatalf("library saLib = %s not defined to compute SA\n%s\n", saLib, s)
	}
	lcp := Lcp(s, sa)
	// Add last element to lcp if necessary
	if lcp[len(lcp)-1] != -1 {
		lcp = append(lcp, -1)
	}
	cld := Cld(lcp)

	return Esa{s, sa, lcp, cld, strandSize}
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



func SaDivSufSort(t []byte) []int {
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

func Lcp(t []byte, sa []int) []int {
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