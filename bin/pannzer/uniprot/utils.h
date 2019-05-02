#include <string>
#include <cstring>
#include <vector>
#include <unordered_map>

using namespace std;

/*
 * UTILITY FUNCTIONS
 */
 
 /*
  * Splits string at occurences of one or more delimiters.
  * Two or more delimiters in a sequence are interpreted as
  * a single splitting point.
  *
  * NOTE: NOT suited for splitting data-lines with empty data-fields
  *
  * based on: strtok()
  */
vector<string> split(const string str, const char* del){
	
	char str_copy[1000];
	strcpy(str_copy,str.c_str());
	
	vector<string> sp;
	char* tok= strtok(str_copy,del);
	while(tok != NULL){
		sp.push_back(string(tok));
		tok= strtok(NULL,del);
	}
	return sp;
}

/*
 * Splits string at every occurance of delimiters del.
 * Two delimiters in a sequence are interpreted as defining an empty string.
 *
 * NOTE: suited for splitting data-lines with empty entries
 *
 * based on: strpbrk()
 */
vector<string> split2(char* str, const char* del){

	char* str_start= str;
	
	vector<string> sp;
	char* pbrk = strpbrk(str_start,del);
	while(pbrk != NULL){
		*pbrk= '\0';
		sp.push_back( string(str_start) );
		
		str_start= pbrk+1;
		pbrk= strpbrk(str_start,del);
	}
	sp.push_back( string(str_start));

	return sp;
}

/* 
 * split2() for <string> input
 */
vector<string> split2(const string str, const char* del){
	char str_c[1000];
	strcpy(str_c,str.c_str());
	return split2(str_c,del);
}

/*
 * UTILS FOR VECTORS, VECTOR BASED MATRICES AND HASH MAPS
 *
 */
 
/*
 * Constructs and returns a prefix from the first k elements of V.
 * If V has fewer than k elements, constructs and returns a copy of V.
 */
 template <typename T>
 vector<T> get_prefix(const vector<T> &V, unsigned int k){
 	unsigned int prefix_length;
	if(V.size()<k)
		prefix_length= V.size();
	else
		prefix_length= k;
		
	vector<T> prefix(V.begin(),V.begin()+prefix_length);
	return prefix;
}

/*
 * Returns keys of a hash map in a vector.
 * The hash is unchanged.
 *
 */
template <typename T1,typename T2>
vector<T1> get_keys(unordered_map<T1,T2> &mymap){

	vector<T1> keys;
	
	for (auto it = mymap.begin(); it != mymap.end(); ++it ){
		keys.push_back(it->first);
	}
	sort(keys.begin(),keys.end());

	return keys;
}
	

/*
 * Returns values of a hash map in a vector.
 * The hash is unchanged.
 *
 */
template <typename T1,typename T2>
vector<T2> get_values(unordered_map<T1,T2> &mymap){

	vector<T1> keys;
	vector<T2> values;
	
	for (auto it = mymap.begin(); it != mymap.end(); ++it ){
		keys.push_back(it->first);
	}
	sort(keys.begin(),keys.end());
	
	for(unsigned int i=0; i<keys.size(); i++){
		values.push_back(mymap[keys[i]]);
	}
	
	
	return values;
}
