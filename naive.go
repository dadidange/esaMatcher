package esaMatcher

import "sort"


type saSorter struct{
	sa []int
	str []byte
	size int
}

func (m saSorter) Len() int{
	return m.size
}

func (m saSorter) Less(i,j int) bool{
	return less(m.sa[i], m.sa[j],m.size, m.str)
}

func less(i,j,n int, str[]byte) bool{
	if i == n || j == n {
		return i == n
	} else if str[i] != str[j]{
		return str[i] < str[j]
	} else { 
		return less(i+1, j+1, n, str) 
	}
}

func (m saSorter) Swap(i,j int){
	m.sa[i], m.sa[j] = m.sa[j], m.sa[i]
}

func saNaive(t []byte) []int {
	n := len(t)
	indeces := make([]int, n)

	for idx := range t{
		indeces[idx] = idx
	}
	sa := saSorter{indeces, t, n}
	
	sort.Sort(sa)
	return sa.sa
}
