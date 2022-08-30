package esaMatcher

import (
	"bufio"
	"fmt"
	"index/suffixarray"
	"math/rand"
	"os"
	"strings"
	"testing"
)

var path = "testSeqs/AE005174.fa"
var ecoSeq []byte
var ranseq50KBP []byte
var ranseq50MBP []byte
var ranseq5MBP []byte
var aseq5MBP []byte

func TestMain(m *testing.M) {
	// call flag.Parse() here if TestMain uses flags
	ranseq50MBP = ranseq(50000000, "ACGT")
	ranseq5MBP = ranseq(5000000, "ACGT")
	ranseq50KBP = ranseq(50000, "ACGT")
	aseq5MBP = ranseq(5000, "T")
	ecoSeq = readFile(path)
	os.Exit(m.Run())
}

//-----------------------------
// Functionality Tests
//-----------------------------

// func TestDummy(t *testing.T) {
// 	x := []byte("ACAAACATAT")
// 	e := NewEsa(x, "SaNaive")
// 	e.Print(0)
// 	e = NewEsa(x, "SaSais")
// 	e.Print(0)
// }

func TestSasComp(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACAAACATAT"), //OhleBusch Book
		[]byte("ACTTCACAAA"), //ranseq
		[]byte("AAGTAAGG"),   //phylonium & andi
		[]byte("TCTAATGAATATGTAGGATACGAATCGGATGATTCTCGAGCATCAGCATGGA " +
			"GGACTCCATGTGTTAAGTCACCGACTGCGTGCCACCTGCGTCTTCGAAACGG"),
		ecoSeq,
	}

	for _, seq := range seqs {
		saSais := Sa(seq, "SaSais")
		saDiv := Sa(seq, "SaDivSufSort")
		saNaive := Sa(seq, "SaNaive")
		saGo := goSa(seq)

		if len(saSais) != len(saDiv) ||
			len(saNaive) != len(saGo) ||
			len(saSais) != len(saNaive) {
			t.Error("Size Mismatch for Suffix Arrays")

		}
		n := len(saSais)
		for i := 0; i < n; i++ {
			if saSais[i] != saDiv[i] ||
				saNaive[i] != saSais[i] {
				t.Errorf("Idx Mismatch for Suffix Arrays at %d: %d vs. %d vs. %d", i, saSais[i], saDiv[i], saNaive[i])
				return
			}
		}

	}

}
func TestChildArray(t *testing.T) {
	s1 := []byte("ACTTCACAAA") //ranseq
	s2 := []byte("ACAAACATAT") //OhleBusch Book
	s3 := []byte("AAGTAAGG")   //phylonium & andi
	s4 := []byte("TCTAATGAATATGTAGGATACGAATCGGATGATT" +
		"CTCGAGCATCAGCATGGAGGACTCCATGTGTTAAG" +
		"TCACCGACTGCGTGCCACCTGCGTCTTCGAAACGG")

	s := [][]byte{s1, s2, s3, s4}

	var e Esa

	for _, seq := range s {
		e = NewEsa(seq, "")
		lcp := e.lcp
		if len(e.s) != (len(lcp) - 1) {
			t.Errorf("Size mismatch, len(e.s) must be len(lcp)-1, we got %d & %d",
				len(e.s), len(lcp))
		}

		if len(e.cld) != len(lcp) {
			t.Errorf("Size mismatch, len(e.cld) must be len(lcp)-1, we got %d & %d",
				len(e.cld), len(lcp))
		}
		for i := 0; i < len(seq)-1; i++ {
			if lcp[i] > lcp[i+1] { //cld(i) contains cld(i+1).L
				l := lcp[e.cld[i]]
				for k := e.cld[i]; k < i; k++ {
					if l > lcp[k] {
						t.Errorf("l=%d should be min in interval[%d..%d] \nbut %d is@ %d",
							l, e.cld[i], i, lcp[k], k)
					}
				}
			} else { //cld(i) contains cld(i).R
				l := lcp[i]
				for k := i; k < e.cld[i]; k++ {
					if l > lcp[k] {
						t.Errorf("l=%d should be min in interval[%d..%d] \nbut %d is@ %d",
							l, i, e.cld[i], lcp[k], k)
					}
				}
			}
		}
	}
}

func TestGetInterval_Single(t *testing.T) {
	seq := []byte("ACTTCACAAA") //ranseq
	e := NewEsa(seq, "")

	nucleotides := []byte("ATGC")
	m := len(e.Sequence()) - 1
	i := NewEsaInterval(0, m, e)

	for _, c := range nucleotides {
		cld := e.GetInterval(i, c)
		for k := cld.start; k < cld.end; k++ {
			if e.s[e.sa[k]] != c {
				t.Errorf("Found child that does not start with %s in Interval{%d,%d} @%d\n",
					string(c), cld.start, cld.end, k)
			}
		}
	}
}

func TestGetInterval_OwnSuffixes(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACAAACATAT"), //OhleBusch Book
		[]byte("ACTTCACAAA"), //ranseq
		[]byte("AAGTAAGG"),   //phylonium & andi
		[]byte("TCTAATGAATATGTAGGATACGAATCGGATGATTCTCGAGCATCAGCATGGA" +
			"GGACTCCATGTGTTAAGTCACCGACTGCGTGCCACCTGCGTCTTCGAAACGG"),
	}

	var e Esa
	var inter EsaInterval
	var subI EsaInterval

	verbose := false

	for seqIdx, seq := range seqs {
		e = NewEsa(seq, "")
		n := len(seq)
		for subS := 0; subS < n; subS++ {
			//remove first character to get next substring/suffix
			if subS > 0 {
				seq = seq[1:]
			}
			if verbose {
				fmt.Printf("Testing %d for %s\n", seqIdx, string(seq))
			}

			m := len(e.Sequence()) - 1
			inter = NewEsaInterval(0, m, e)
			i := 0
			for i < len(seq) {
				// fmt.Printf("Checking %d=%s\n", i, string(seq[i]))
				subI = e.GetInterval(inter, seq[i])
				if verbose {
					fmt.Printf("{%d,%d, %d} for %s returned {%d,%d,%d}\n",
						inter.start, inter.end, inter.mid, string(seq[i]), subI.start, subI.end, subI.mid)
				}
				if subI.start == subI.end {
					// Two options -> only the matching sequence left
					if subI.start < 0 {
						t.Errorf("\nempty interval returned for %s, index %d at Interval {%d,%d}\n",
							string(seq[i]), i, inter.start, inter.end)
						break
					}
					if e.s[e.sa[subI.start]+i] == seq[i] {
						for k := i; k < len(seq); k++ {
							if e.s[e.sa[subI.start]+k] != seq[k] {
								t.Errorf("\n Mismatch found at %d in %s\n", i+k, string(seq))
							}
						}
						//No mismatching elements found in only remaining suffix
						break
					}
				}
				i++
				//Skipping characters that are equal in lcp
				l := min(len(seq), subI.l)
				for ; i < l; i++ {
					if e.s[e.sa[subI.start]+i] != seq[i] {
						err := fmt.Sprintf("\nMismatch found at %d", i)
						t.Errorf("%s\nS:%s\n", err, string(seq))
					}
				}
				inter = subI
			}
		}
	}
}

func TestGetMatch(t *testing.T) {
	seqs := [][]byte{
		[]byte("ACAAACATAT"), //OhleBusch Book
		[]byte("ACTTCACAAA"), //ranseq
		[]byte("AAGTAAGG"),   //phylonium & andi
		[]byte("TCTAATGAATATGTAGGATACGAATCGGATGATTCTCGAGCATCAGCATGGA " +
			"GGACTCCATGTGTTAAGTCACCGACTGCGTGCCACCTGCGTCTTCGAAACGG"),
	}
	verbose := false

	var e Esa
	if verbose {
		t.Logf("\nTesting own Subsequences\n")
	}
	for _, seq := range seqs {
		e = NewEsa(seq, "")
		n := len(seq)
		for s := 0; s < n; s++ {
			seq = seq[1:]
			if verbose {
				t.Logf("\nMatching %s to %s\n", string(seq), string(e.s))
			}
			if string(seq) == "TAT" {
				verbose = true
			}
			m := e.GetMatch(seq)
			if m.start == -1 && m.end == -1 {
				t.Errorf("Mismatch found for %s in %s\n", string(seq), string(e.s))
			}

			for i := m.start; i <= m.end; i++ {
				refStr := string(e.s[e.sa[i]:])
				qStr := string(seq)
				if !strings.Contains(refStr, qStr) {
					t.Errorf("\nFound %s to contain %s\n", refStr, qStr)
				}
			}
		}
	}

	if verbose {
		t.Logf("\nTesting Mismatches\n")
	}
	seq := []byte("TCTAATGAATATGTAGGATACGAATCGGATGATTCTCGAGCATCAGCATGGA " +
		"GGACTCCATGTGTTAAGTCACCGACTGCGTGCCACCTGCGTCTTCGAAACGG")
	e = NewEsa(seq, "")
	q := []byte("M")
	m := e.GetMatch(q)
	if m.start != -1 || m.end != -1 {
		t.Errorf("Match found for %s in %s", string(q), string(e.s))
	}

}

//-----------------------------
//Benchmarks
//-----------------------------

func Benchmark_Sa_LibSais_50KBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq50KBP, "SaSais")
	}
}

func Benchmark_Sa_LibDivSufSort_50KBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq50KBP, "SaDivSufSort")
	}
}

func Benchmark_Sa_GoSa_50KBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		goSa(ranseq50KBP)
	}
}

func Benchmark_Sa_Naive_50KBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq50KBP, "SaNaive")
	}
}

func Benchmark_Sa_LibSais_5MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq5MBP, "SaSais")
	}
}

func Benchmark_Sa_LibDivSufSort_5MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq5MBP, "SaDivSufSort")
	}
}

func Benchmark_Sa_GoSa_5MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		goSa(ranseq5MBP)
	}
}

func Benchmark_Sa_Naive_5MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq5MBP, "SaNaive")
	}
}

// func Benchmark_Sa_LibSais_T_5MBP(b *testing.B) {
// 	for n := 0; n < b.N; n++ {
// 		Sa(aseq5MBP, "SaSais")
// 	}
// }

// func Benchmark_Sa_LibDivSufSort_T_5MBP(b *testing.B) {
// 	for n := 0; n < b.N; n++ {
// 		Sa(aseq5MBP, "SaDivSufSort")
// 	}
// }

// func Benchmark_Sa_Naive_T_5MBP(b *testing.B) {
// 	for n := 0; n < b.N; n++ {
// 		Sa(aseq5MBP, "SaNaive")
// 	}
// }

// func Benchmark_Sa_GoSa_T_5MBP(b *testing.B) {
// 	for n := 0; n < b.N; n++ {
// 		goSa(aseq5MBP)
// 	}
// }

func Benchmark_Sa_LibSais_50MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq50MBP, "SaSais")
	}
}

func Benchmark_Sa_LibDivSufSort_50MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq50MBP, "SaDivSufSort")
	}
}

func Benchmark_Sa_GoSa_50MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		goSa(ranseq50MBP)
	}
}

func Benchmark_Sa_Naive_50MBP(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ranseq50MBP, "SaNaive")
	}
}

func Benchmark_Sa_LibSais_Eco(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ecoSeq, "SaSais")
	}
}

func Benchmark_Sa_LibDivSufSort_Eco(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ecoSeq, "SaDivSufSort")
	}
}

func Benchmark_Sa_GoSa_Eco(b *testing.B) {
	for n := 0; n < b.N; n++ {
		goSa(ecoSeq)
	}
}

func Benchmark_Sa_Naive_Eco(b *testing.B) {
	for n := 0; n < b.N; n++ {
		Sa(ecoSeq, "SaNaive")
	}
}

func Benchmark_Esa_LibSais_Eco(b *testing.B) {
	for n := 0; n < b.N; n++ {
		NewEsa(ecoSeq, "SaSais")
	}
}

func Benchmark_Esa_LibDivSufSort_Eco(b *testing.B) {
	for n := 0; n < b.N; n++ {
		NewEsa(ecoSeq, "SaDivSufSort")
	}
}

//-----------------------------
// Helper Functions
//-----------------------------

//generate sequence with length len and the given alphabet nuc
func ranseq(seqLen int, nuc string) []byte {
	n := len(nuc)
	seq := make([]byte, seqLen)
	for i := 0; i < seqLen; i++ {
		rd := rand.Int() % n
		//  nuc at random position
		seq[i] = nuc[rd]
	}
	return seq
}

//read fasta file from input string
func readFile(path string) []byte {
	fileHandle, _ := os.Open(path)
	defer fileHandle.Close()
	sc := bufio.NewScanner(fileHandle)
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

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func goSa(data []byte) []byte {
	index := suffixarray.New(data)
	return index.Bytes()
}
