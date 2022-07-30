package esaMatcher

import (
	"bufio"
	"math/rand"
	// "fmt"
	"os"
	"testing"
)


var path="testSeqs/AE005174.fa"
var ecoSeq []byte
var ranseq50MBP []byte
var ranseq5MBP []byte
var result []int

func TestMain(m *testing.M) {
	// call flag.Parse() here if TestMain uses flags
	ranseq50MBP = ranseq(50000000, "ACGT")
	ranseq5MBP = ranseq(5000000, "ACGT")
	ecoSeq = readFile(path)
	os.Exit(m.Run())
}

// func TestDebug(t *testing.T){
// 	s1 := []byte("ACTTGACAA")//ranseq
// 	s2 := []byte("ACAAACATAT")//OhleBusch Book
// 	s3 := []byte("AAGTAAGG")
	
// 	s := [][]byte{s1,s2,s3}
	 
// 	var e Esa

// 	for _,seq := range s {
// 		e = NewEsa(seq, "SaSais")
// 		fmt.Printf("\nESA for %s\n", seq)
// 		e.Print(0)
// 	}
// }

// func TestDebug2(t *testing.T){
// 	s1 := []byte("ACTTGACAA")//ranseq
	
	 
// 	s := Sa(s1, "SaSais")
// 	fmt.Println(s)

// 	s = Sa(s1, "SaDivSufSort")
// 	fmt.Println(s)
// }


func BenchmarkLibSais5MBP(b *testing.B){
	s := Sa(ranseq5MBP, "SaSais")
	result=s
}


func BenchmarkLibDivSufSort5MBP(b *testing.B){
	s := Sa(ranseq5MBP, "SaDivSufSort")
	result=s
}

func BenchmarkLibSais50MBP(b *testing.B){
	s := Sa(ranseq50MBP, "SaSais")
	result=s
}


func BenchmarkLibDivSufSort50MBP(b *testing.B){
	s := Sa(ranseq50MBP, "SaDivSufSort")
	result=s
}

func BenchmarkSaisEco(b *testing.B){
	s := Sa(ecoSeq, "SaSais")
	result=s
}


func BenchmarkLibDivSufSortEco(b *testing.B){
	s := Sa(ecoSeq, "SaDivSufSort")
	result=s
}

//helpers

//generate sequence with length len and the given alphabet nuc
func ranseq(len int, nuc string) []byte{
	seq := make([]byte, len)
	for i := 0; i < len; i++ {
		rd := rand.Int() % 4
		//  nuc at random position
		seq[i] = nuc[rd]
	}
	return seq
}

func readFile(path string) []byte{
	fileHandle, _ := os.Open(path)
	defer fileHandle.Close()
	sc  := bufio.NewScanner(fileHandle)
	var t []byte
	var seq []byte

	for sc.Scan() {
		t = sc.Bytes()
		//Skip enmpty lines or fasta header
		if len(t) > 0 && t[0] != '>' {
			seq = append(seq, t...)
		}
	}
	return seq
}