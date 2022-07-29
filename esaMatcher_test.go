package esaMatcher

import (
	"fmt"
	"testing"
)

func TestDebug(t *testing.T){
	s1 := []byte("ACTTGACAA")//ranseq
	s2 := []byte("ACAAACATAT")//OhleBusch Book
	s3 := []byte("AAGTAAGG")
	
	s := [][]byte{s1,s2,s3}
	 
	var e Esa

	for _,seq := range s {
		e = NewEsa(seq, "SaSais")
		fmt.Printf("\nESA for %s\n", seq)
		e.Print(0)
	}
}

func TestDebug2(t *testing.T){
	s1 := []byte("ACTTGACAA")//ranseq
	
	 
	s := Sa(s1, "SaSais")
	fmt.Println(s)

	s = Sa(s1, "SaDivSufSort")
	fmt.Println(s)
}



