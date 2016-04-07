#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
IntegerVector binary_search(IntegerVector breaks, IntegerVector search, bool left_index) {
	
	int first, last, middle;
	int n = breaks.size();

	int n2 = search.size();
	IntegerVector index(n2);

	for(int i = 0; i < n2; i ++) {
		first = 0;
		last = n - 1;
		middle = (first+last)/2;
		// printf("first = %i, middle = %i, last = %i\n", first, middle, last);

		if(search[i] < breaks[0]) {
			if(!left_index) {
				index[i] = 0;
			} else {
				index[i] = NA_INTEGER;
			}
			continue;
		}
		if(search[i] > breaks[last]) {
			if(left_index) {
				index[i] = n - 1;
			} else {
				index[i] = NA_INTEGER;
			}
			continue;
		}
		if(n <= 1) {
			search[i] = NA_INTEGER;
			continue;
		}

		while (first <= last) {
			if(breaks[first] == search[i]) {
				index[i] = first;
				break;
			} else if(breaks[last] == search[i]) {
				index[i] = last;
				break;
			} else if(last - first == 1) {
				if(left_index) {
					index[i] = first;
				} else {
					index[i] = last;
				}
				break;
			} else if (breaks[middle] == search[i]) {
				index[i] = middle;
				break;
			} else if (breaks[middle] < search[i]) {
				first = middle;    
			} else {
				last = middle;
			}
			middle = (first + last)/2;
			// printf("first = %i, middle = %i, last = %i\n", first, middle, last);
		}
	}

	return index;
}

// [[Rcpp::export]]
IntegerVector extract_sites(IntegerVector start, IntegerVector end, IntegerVector site, bool return_index, int min_sites) {

	IntegerVector index1 = binary_search(site, start, false);
	IntegerVector index2 = binary_search(site, end, true);

	int len = 0;
	int n = index1.size();
	IntegerVector passed_index(n);
	int l;

	for(int i = 0; i < n; i ++) {
		if(IntegerVector::is_na(index1[i]) || IntegerVector::is_na(index2[i])) {
			continue;
		}
		if(index2[i] < index1[i]) {
			continue;
		}

		l = index2[i] - index1[i] + 1;

		if(l >= min_sites) {
			len += l;
			passed_index[i] = i;
		}
	}

	IntegerVector out(len);
	if(len == 0) {
		return out;
	}

	int k = 0;
	if(return_index) {
		for(int i = 0; i < n; i ++) {
			if(!IntegerVector::is_na(passed_index[i])) {
				for(int j = index1[i]; j <= index2[i]; j ++) {
					out[k] = j;
					k ++;
				}
			}
		}
	} else {
		for(int i = 0; i < n; i ++) {
			if(!IntegerVector::is_na(passed_index[i])) {
				for(int j = index1[i]; j <= index2[i]; j ++) {
					out[k] = site[j];
					k ++;
				}
			}
		}
	}

	return out;
}
