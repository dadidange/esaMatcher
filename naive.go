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
	// return less(m.sa[i], m.sa[j],m.size, m.str)
	istr := m.str[m.sa[i]:]
	jstr := m.str[m.sa[j]:]
	for n := 0; n < len(istr) && n < len(jstr); n++ {
		if istr[n] == jstr[n]{
			continue
		} else {
			return istr[n] < jstr[n]			
		}
	}
	return len(istr) < len(jstr)
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
