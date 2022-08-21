package esaMatcher

import "sort"


type saSorter struct{
	idx []int
	str []byte
}

func (m saSorter) Len() int{
	return len(m.idx) //is the same as len(m.Names)
}

func (m saSorter) Less(i,j int) bool{
	return less(i,j, m.str)
}

func less(i,j int, str[]byte) bool{
	if str[i] != str[j]{
		return str[i] < str[j]
	} else {
		return less(i+1, j+1, str)
	}

}

func (m saSorter) Swap(i,j int){
	m.idx[i], m.idx[j] = m.idx[j], m.idx[i]
}

func saNaive(t []byte) []int {
	n := len(t)
	indeces := make([]int, n + 1)
	t = append(t, '$')

	for idx := range t{
		indeces[idx] = idx
	}
	sa := saSorter{indeces, t}
	
	sort.Sort(sa)
	return sa.idx
}
