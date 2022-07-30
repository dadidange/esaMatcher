package esaMatcher

import (
	"bufio"
	"math/rand"
	"os"
	"testing"
)


var path="testSeqs/AE005174.fa"
var ecoSeq []byte
var ranseq50MBP []byte
var ranseq5MBP []byte
var result []int
var esaRes Esa

func TestMain(m *testing.M) {
	// call flag.Parse() here if TestMain uses flags
	ranseq50MBP = ranseq(50000000, "ACGT")
	ranseq5MBP = ranseq(5000000, "ACGT")
	ecoSeq = readFile(path)
	os.Exit(m.Run())
}

//-----------------------------
//Benchmarks
//-----------------------------

func BenchmarkLibSais5MBP(b *testing.B){
	for n:=0; n<b.N; n++{
		Sa(ranseq5MBP, "SaSais")
	}
}

func BenchmarkLibDivSufSort5MBP(b *testing.B){
	for n:=0; n<b.N; n++{
		Sa(ranseq5MBP, "SaDivSufSort")
	}
}

func BenchmarkLibSais50MBP(b *testing.B){
	for n:=0; n<b.N; n++{
		Sa(ranseq50MBP, "SaSais")
	}
}

func BenchmarkLibDivSufSort50MBP(b *testing.B){
	for n:=0; n<b.N; n++{
		Sa(ranseq50MBP, "SaDivSufSort")
	}
}

func BenchmarkSaisEco(b *testing.B){
	for n:=0; n<b.N; n++{
		Sa(ecoSeq, "SaSais")
	}
}

func BenchmarkLibDivSufSortEco(b *testing.B){
	for n:=0; n<b.N; n++{
		Sa(ecoSeq, "SaDivSufSort")
	}
}

func BenchmarkEsaSaisEco(b *testing.B) {
	for n:=0; n<b.N; n++{
		NewEsa(ecoSeq, "SaSais")
	}
}

func BenchmarkEsaDivSufSortEco(b *testing.B) {
	for n:=0; n<b.N; n++{
		NewEsa(ecoSeq, "SaDivSufSort")
	}
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

//read fasta file from input string
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