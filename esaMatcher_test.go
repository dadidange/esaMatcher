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
		e = NewEsa(seq, "SaDivSufSort")
		fmt.Printf("\nESA for %s\n", seq)
		e.Print(0)
	}
}



