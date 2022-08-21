package esaMatcher

import "sort"


type saSorter struct{
	sa []int
	idx []int
	str []byte
}

func (m saSorter) Len() int{
	return len(m.idx) //is the same as len(m.Names)
}

func (m saSorter) Less(i,j int) bool{
	return less(m.sa[i], m.sa[j], m.str)
}

func less(i,j int, str[]byte) bool{
	if i == len(str) || j == len(str) {
		return i == len(str)
	} else if str[i] != str[j]{
		return str[i] < str[j]
	} else { 
		return less(i+1, j+1, str) 
	}
}

func (m saSorter) Swap(i,j int){
	m.sa[i], m.sa[j] = m.sa[j], m.sa[i]
}

func saNaive(t []byte) []int {
	n := len(t)
	// if t[n-1] != '$'{
	// 	n += 1
	// 	t =append(t, '$')
	// } 
	indeces := make([]int, n)

	for idx := range t{
		indeces[idx] = idx
	}
	sa := saSorter{indeces, indeces, t}
	
	sort.Sort(sa)
	return sa.idx
}
