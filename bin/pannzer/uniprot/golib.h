#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <cstring>
#include <algorithm>
#include <queue>

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>


// other files
#include <stats.h>
#include <utils.h>

using namespace std;

/*
 * THIS IS A COLLECTION OF METHODS FOR MANIPULATION OF GO ANNOTATION DATA
 *
 * THE DATA IS OF THE FORM <geneid,goid,rank,score>
 */

struct go_annotation{
	string gene;
	unsigned int go;
	unsigned int rank;
	double score;
};


// comparison operators for quick sorting of data
bool compare_gene(const go_annotation &a, const go_annotation &b){
	return (a.gene < b.gene);
}
bool compare_go(const go_annotation &a, const go_annotation &b){
	return (a.go<b.go);
}
bool compare_rank(const go_annotation &a, const go_annotation &b){
	return (a.rank>b.rank);
}
bool compare_score(const go_annotation &a, const go_annotation &b){
	return (a.score>b.score);
}
bool compare_gene_score(const go_annotation &a,const go_annotation &b){
	if(a.gene != b.gene){
		return (a.gene < b.gene);}
	else{
		return (a.score > b.score);}
}
bool compare_go_score (const go_annotation &a,const go_annotation &b){
	if(a.go != b.go){ 
		return (a.go<b.go); }
	else{ 
		return (a.score>b.score);}
}
bool compare_gene_go_score(const go_annotation &a,const go_annotation &b){
	if(a.gene != b.gene){
		return (a.gene < b.gene);}
	else if(a.go != b.go){
		return (a.go < b.go);}
	else{
		return (a.score > b.score);}
}

/*
 * Returns a deep copy
 */
vector<go_annotation> copy(const vector<go_annotation> &data){
	vector<go_annotation> data_new= data;
	return data_new;
}

/*
 * Returns a deep copy
 */
vector<vector<go_annotation> > copy(const vector<vector<go_annotation> > &data){
	vector<vector<go_annotation> > data_new;
	for(unsigned int i=0; i<data.size(); i++){
		vector<go_annotation> gene_data= data[i];
		data_new.push_back(gene_data);
	}
	return data_new;
}

/*
 * Return minimum score
 */
double get_min_score(const vector<go_annotation> &data){
	if(data.size() == 0)
		return 0;
	double min_score	= data[0].score;
	for(unsigned int i=0; i<data.size(); i++){
		if(data[i].score < min_score)
			min_score = data[i].score;
	}
	return min_score;
}
/*
 * Return maximum score
 */
double get_max_score(const vector<go_annotation> &data){
	if(data.size() == 0)
		return 0;
	double max_score	= data[0].score;
	for(unsigned int i=0; i<data.size(); i++){
		if(data[i].score > max_score)
			max_score = data[i].score;
	}
	return max_score;
}

/*
 * Catenates two go_annotation vectors
 */
vector<go_annotation> catenate(const vector<go_annotation> &dat1, const vector<go_annotation> &dat2){
	vector<go_annotation> cat= dat1; // deep copy
	cat.insert(cat.end(),dat2.begin(),dat2.end());
	return cat;
}

/*
 * Reads text file containing go annotation into a vector of <gene,go,rank,score> structures.
 * 
 * Go annotations are ranked by their order number in the input file.
 *
 * After reading data is sorted gene>go>score
 *
 * in:
 * fin		file in the format "gene_id\t[GO:|go:|]7u\tscore\n"
 *
 * out:
 * data		vector of go_annotations <gene,go,rank,score>
 *
 */
vector<go_annotation> read_go_annotation_data(const char* fin){

	ifstream in(fin, ios::in);
	if(!in.is_open()){
		cerr << "failed to open file \'"<<fin<<"\': exiting\n";
		exit(1);
	}
	char buffer[1000];
	string f1;
	string f2;
	string f3;
	unsigned int go;
	double sc;
	
	vector<go_annotation> data;
	
	int line_num=0;
	while(in.getline(buffer,1000)){
		line_num++;
		if(buffer[0] == '#' || buffer[0] == '!')
			continue;
		string line(buffer);
		if(line == "")
			continue;
		istringstream iss (line,istringstream::in);
		iss >> f1;
		iss >> f2;
		iss >> f3;
		if( f2.substr(0,3) == "GO:" || f2.substr(0,3)=="go:"){ // f2.substr(0,3) == "go:"
			f2= f2.substr(3);
		}
		go = atoi(f2.c_str());
		sc = atof(f3.c_str());

		// invalid format
		if(iss.fail() || go==0){
			cerr << "ERROR: invalid format: "<<line_num<<":"<<fin<<"\n";
			exit(1);
		}
		
		go_annotation temp= {f1,go,0,sc};
		data.push_back(temp);
	}
	in.close();
	
	// ranking
	for(unsigned int i=0; i<data.size(); i++)
		data[i].rank= i;
	
	
	// sorting gene>score
	sort(data.begin(),data.end(),compare_gene_score);
	
	return data;
}


/*
 * Restructure vector of go_annotations so that there is one vector of annotations for each gene/go-term.
 * Sort your data by gene/go-term before calling this function.
 *
 * data_in	vector of go_annotations (sorted by gene/go-term)
 * factor	"gene" or "go"
 * data_out	vector containing one vector of go_annotations for each unique gene/go-term
 *
 */
vector<vector<go_annotation> > group_go_annotation_data(const vector<go_annotation> &data_in, const string factor){

	// assumed to be sorted by
	// sort(data_in.begin(),data_in.end(),compare_gene);
	// sort(data_in.begin(),data_in.end(),compare_go);

	vector<vector<go_annotation> > data_out;
	
	if(data_in.size()==0){
		return data_out;
	}
	
	if(factor == "gene"){
		string gene0 = data_in[0].gene;
		vector<go_annotation> data_temp;
		data_temp.push_back(data_in[0]);
		for(unsigned int i=1;i<data_in.size(); i++){	
			if(data_in[i].gene != gene0){
				vector<go_annotation> data_temp_copy= data_temp; // create a copy
				data_out.push_back(data_temp_copy);
				data_temp.clear();
				gene0 = data_in[i].gene;
			}
			data_temp.push_back(data_in[i]);
		}
		if(data_temp.size() > 0){
			vector<go_annotation> data_temp_copy= data_temp; // create a copy
			data_out.push_back(data_temp_copy);
		}
	}
	else if(factor == "go"){
		unsigned int go0 = data_in[0].go;
		vector<go_annotation> data_temp;
		data_temp.push_back(data_in[0]);
		for(unsigned int i=1;i<data_in.size(); i++){	
			if(data_in[i].go != go0){
				vector<go_annotation> data_temp_copy= data_temp; // create a copy
				data_out.push_back(data_temp_copy);
				data_temp.clear();
				go0 = data_in[i].go;
			}
			data_temp.push_back(data_in[i]);
		}
		if(data_temp.size() > 0){
			vector<go_annotation> data_temp_copy= data_temp; // create a copy
			data_out.push_back(data_temp_copy);
		}
	}
	else{
		cout << "# ERROR: group_go_annotation_data:invalid group factor:"<<factor<<"\n";
		exit(1);
	}
	
	return data_out;		
}

void print_go_annotation_data(const vector<go_annotation> &data){
	for(unsigned int i=0; i<data.size(); i++){
		printf("%s\tGO:%07u\t%u\t%.6f\n",data[i].gene.c_str(),data[i].go,data[i].rank,data[i].score);
	}
}

void print_go_annotation_data(const vector<vector<go_annotation> > &data){
	
	for(unsigned int i=0; i<data.size(); i++){
		for(unsigned int j=0; j<data[i].size(); j++){
			printf("%s\tGO:%07u\t%u\t%.6f\n",data[i][j].gene.c_str(),data[i][j].go,data[i][j].rank,data[i][j].score);
		}
		printf("<>\n");
	}
}


/*
 * Reads information content for each go-term into a hash "pmap"
 *
 */
void read_information_content(const char* fin, unordered_map <unsigned int,double>& icmap){

	ifstream in(fin, ios::in);
	if(!in.is_open()){
		cerr << "failed to open file \'"<<fin<<"\': exiting\n";
		exit(1);
	}
	char buffer[1000];
	string f1;
	string f2;
	unsigned int go;
	double ic;
	
	int line_num=0;
	while(in.getline(buffer,1000)){
		line_num++;
		if(buffer[0] == '#' || buffer[0] == '!')
			continue;
		string line(buffer);
		if(line == "")
			continue;
		istringstream iss (line,istringstream::in);
		iss >> f1;
		iss >> f2;
		if( f1.substr(0,3) == "GO:" || f1.substr(0,3)=="go:"){
			f1= f1.substr(3);
		}
		go = atoi(f1.c_str());
		ic = atof(f2.c_str());

		// invalid format
		if(iss.fail() || go==0 || ic<0){
			cerr << "ERROR: invalid format: "<<line_num<<":"<<fin<<"\n";
			exit(1);
		}
		
		icmap[go]= ic;
	}
	in.close();
}

/* 
 * Print information content to file
 *
 */
void print_information_content(const char* fout, unordered_map <unsigned int,double>& icmap){

	FILE * OUT = fopen (fout,"w");
	if(OUT == NULL){
		cerr<< "ERROR: failed to open file: "<<fout<<"\n";
		exit(1);
	}
	
	vector<unsigned int> goid_list= get_keys(icmap);
	fprintf(OUT,"#GO\tIC\n");
	for(unsigned int i=0; i<goid_list.size(); i++){
	
		fprintf(OUT,"%07u\t%1.6f\n",goid_list[i],icmap[goid_list[i]]);
	}
	fclose(OUT);
}


/*
 * Reads go parent nodes to a hash "pmap" given by reference
 * 
 * Parent nodes of node i can be accessed as pmap[i]
 *
 * Parent nodes for each child are sorted by ascending order (for set operations)
 */
void read_goparents(const char* fin, unordered_map <unsigned int,vector<unsigned int> >& pmap){
	ifstream in(fin, ios::in);
	if(!in.is_open()){
		cerr << "failed to open file \'"<<fin<<"\': exiting\n";
		exit(1);
	}
	char buffer[10000];
	unsigned int child;
	unsigned int parent;
	int n=0;
	char* pch;
	
	while(in.getline(buffer,10000)){
		n++;
		if(buffer[0] == '#' || buffer[0] == '!')
			continue;
		string line(buffer);
		if(line == "")
			continue;
		
		
		// parsing child
		pch= strtok(buffer,"\t");
		if(pch==NULL){
			cerr<< "# ERROR: "<<n<<":"<<fin<<"\n";
			//exit(1);
			continue;
			}
		child= atoi(pch);
		if(child == 0){
			cerr << "# ERROR: "<<n<<":"<<fin<<"\n";
			//exit(1);
			continue;
		}
		
		// parsing parents
		vector<unsigned int> parents;
		pch= strtok(NULL,",;");
		if(pch==NULL){
			cerr<< "# ERROR: "<<n<<":"<<fin<<"\n";
			//exit(1);
			continue;
		}
  		while (pch != NULL){
			parent= atoi(pch);
			if(parent == 0){
				cerr << "# ERROR: "<<n<<":"<<fin<<"\n";
				//exit(1);
				continue;
			}
			
			parents.push_back(parent);
    			pch = strtok (NULL, ",;");
		}
		sort(parents.begin(),parents.end());
		pmap[child]= parents;
	}
	in.close();
}

/*
 * Prints parent nodes for first n keys in pmap
 */
void print_goparents(unordered_map <unsigned int,vector<unsigned int> >& pmap,unsigned int n){
	vector<unsigned int> keys;
	for (auto it = pmap.begin(); it != pmap.end(); ++it ){
		keys.push_back(it->first);
	}
	sort(keys.begin(),keys.end());
	
	for(unsigned int i=0; i<keys.size() && i<n; i++){
		vector<unsigned int> gos = pmap[ keys[i] ];
		
		printf("%07u\t",keys[i]);
		if(gos.size()>0)
			printf("%07u",gos[0]);
		for(unsigned int j=1; j<gos.size(); j++){
			printf(",%07u",gos[j]);
		}
		printf("\n");
	}
}






/*
 * Adds parent GO terms for each go_annotation in input data.
 * Added go_annotations are set with the same rank/score as their children.
 * If several child nodes have the same parent, only the one with the highest scoring child is added.
 *
 * data		vector of go_annotations
 *
 * pmap		map from go-term to vector of parents
 * 
 */
void propagate_parents(vector<go_annotation> &data, unordered_map<unsigned int,vector<unsigned int> > &pmap){
	
		unsigned int data_size= data.size();
		
		for(unsigned int i=0; i<data_size; i++){
			
			if(pmap.count(data[i].go)==0){
				fprintf(stderr,"ERROR: propagate_parents: no pmap[GO:%07u]\n",data[i].go);
				exit(1);
			}
			
			vector<unsigned int> parent_gos= pmap[data[i].go];
			
			for(unsigned int k=0; k<parent_gos.size(); k++){
				// parent gets its childs rank and score
				go_annotation temp= {data[i].gene,parent_gos[k],data[i].rank,data[i].score}; 
				data.push_back(temp);
			}
		}
		
		// select the top scoring term for each <gene,go> tuple
		//    sort gene > go > score
		sort(data.begin(), data.end(), compare_gene_go_score);
		
		//    select best score for each unique <gene,go> tuple
		vector<go_annotation> data_new;
		data_new.push_back(data[0]);

		string prev_gene= data[0].gene;
		unsigned int prev_go= data[0].go;
		for(unsigned int j=1; j<data.size(); j++){
			if(prev_gene == data[j].gene && prev_go == data[j].go){
				continue;}
			else{
				data_new.push_back(data[j]);
				prev_gene= data[j].gene;
				prev_go= data[j].go;
			}
		}
		
		// sort gene>score
		sort(data_new.begin(), data_new.end(), compare_gene_score);
		
		data= data_new;
}

vector<vector<double> > get_sim_matrix(const vector<unsigned int> &pred,
				const vector<unsigned int> &pos,
				unordered_map<unsigned int,vector<unsigned int> > &pmap,
				unordered_map<unsigned int,double> &icmap,
				const string simf){
				
	// SANITY CHECK
	if( !(simf=="pJacc" || simf=="Resnik" || simf=="Lin" || simf=="PathLin" || simf=="JC" || simf=="id") ){
		cerr << "ERROR: get_sim_matrix(): invalid simf:"<<simf<<"\n";
		exit(1);
	}
	
	vector<vector<double> > sim_matrix(pred.size());
	for(unsigned int i=0; i<pred.size(); i++){
		vector<double> sim_row(pos.size(),0);
		sim_matrix[i]= sim_row;
	}
			
	for(unsigned int i=0; i<pred.size(); i++){
	
		for(unsigned int j=0; j<pos.size(); j++){
			
			if(pmap.count(pred[i])==0){
				fprintf(stderr,"ERROR: get_sim_matrix(): no parents for pred GO:%07u]\n",pred[i]);
				exit(1);
			}
			if(pmap.count(pos[j])==0){
				fprintf(stderr,"ERROR: get_sim_matrix(): no parents for pos GO:%07u]\n",pos[j]);
				exit(1);
			}						
			
			
			// Resnik,Lin,PathLin & JC
			double MICA_IC = 0;
			double N1_IC= 0;
			double N2_IC= 0;
			if(simf=="Resnik" || simf=="Lin" || simf=="PathLin" || simf=="JC"){
				vector<unsigned int> CA= set_intersection(pmap[pred[i]],pmap[pos[j]]); // Common Ancestors
				vector<double> CA_IC;
				for(unsigned int k=0; k<CA.size(); k++){
					if( icmap.count(CA[k])==0 ){
						fprintf(stderr,"ERROR: get_sim_matrix(): no IC for GO:%07u\n",CA[k]);
						exit(1);
					}				
					CA_IC.push_back( icmap[CA[k]]);
				}
				if( CA.size() == 0 )
					MICA_IC = 0;
				else{
					MICA_IC = max(CA_IC);
				}
			}
			
			
			// Lin,PathLin & JC
			if(simf=="Lin" || simf=="PathLin" || simf=="JC"){
				if( icmap.count(pred[i]) == 0){
					fprintf(stderr,"ERROR: get_sim_matrix(): no IC for GO:%07u\n",pred[i]);
					exit(1);
				}
				if( icmap.count(pos[j]) == 0){
					fprintf(stderr,"ERROR: get_sim_matrix(): no IC for GO:%07u\n",pos[j]);
					exit(1);
				}
					
				N1_IC= icmap[pred[i]];
				N2_IC= icmap[pos[j]];
			}			
			

			if(simf=="pJacc"){
				sim_matrix[i][j]= get_Jaccard(pmap[pred[i]],pmap[pos[j]]);
			}
			else if(simf=="Resnik"){
				sim_matrix[i][j]= MICA_IC;
			}
			else if(simf=="Lin"){
				double Lin= 0;
				if( (N1_IC+N2_IC)>0 )
					Lin= 2.0*MICA_IC/(N1_IC+N2_IC);	
				sim_matrix[i][j]= Lin;				
			}
			else if(simf=="JC"){
				double JC= 1.0+2.0*MICA_IC - N1_IC - N2_IC; // TODO: check this
				
				sim_matrix[i][j]= JC;
			}
			else if(simf=="PathLin"){
				// IC(N1_ONLY_PARENTS UNION N2_ONLY_PARENTS) = IC(PATH_SET)
				vector<unsigned int> N1_only_parents= set_difference(pmap[pred[i]],pmap[pos[j]]);
				vector<unsigned int> N2_only_parents= set_difference(pmap[pos[j]],pmap[pred[i]]);
				vector<unsigned int> path_set= N1_only_parents;
				path_set.insert(path_set.end(),N2_only_parents.begin(),N2_only_parents.end());
				double PATH_SET_IC= 0;
				for(unsigned int k=0; k<path_set.size(); k++){
					if( icmap.count(path_set[k]) == 0 ){
						fprintf(stderr,"ERROR: get_sim_matrix(): no IC for GO:%07u\n",path_set[k]);
						exit(1);
					}
					PATH_SET_IC += icmap[path_set[k]];
				}
				
				double PathLin= 0;
				if( (PATH_SET_IC) > 0){
					PathLin = 2.0*MICA_IC/(PATH_SET_IC);
				}
				
				sim_matrix[i][j]= PathLin;
			}
			else if(simf=="id"){
				if( pred[i] == pos[j])
					sim_matrix[i][j]= 1;
				else
					sim_matrix[i][j]= 0;	
			}
			else{
				cerr << "ERROR: get_sim_matrix: invalid simf:"<<simf<<"\n";
				exit(1);
			}
		}
	}
	
	return sim_matrix;
}

/*
 * Matches gene labels in pred and pos data.
 * 
 * Saves matching data to data_pred_bygene_match and data_pos_bygene_match.
 * 
 * If gene labels mismatch significantly: prints warning
 *
 * If all gene labels mismatch: prints error and exits.
 *
 */
void match_gene_labels(const vector<vector<go_annotation> > &data_pred_bygene,
			const vector<vector<go_annotation> > &data_pos_bygene,
			vector<vector<go_annotation> > &data_pred_bygene_match,
			vector<vector<go_annotation> > &data_pos_bygene_match,
			bool verbal=false){
			
	unsigned int pred_size= data_pred_bygene.size();
	unsigned int pos_size= data_pos_bygene.size();
	unsigned int intersect_size= 0;
	double gene_rc= 0;
	double gene_pr= 0;
	unsigned int gi_pred= 0; // gene indices
	unsigned int gi_pos= 0;
	
	data_pred_bygene_match.clear();
	data_pos_bygene_match.clear();
	
	while(gi_pred<pred_size && gi_pos<pos_size){
		if(data_pred_bygene[gi_pred][0].gene < data_pos_bygene[gi_pos][0].gene)
			gi_pred++;
		else if(data_pred_bygene[gi_pred][0].gene > data_pos_bygene[gi_pos][0].gene)
			gi_pos++;
		else{
			vector<go_annotation> pred = data_pred_bygene[gi_pred];
			vector<go_annotation> pos  = data_pos_bygene[gi_pos];
			
			data_pred_bygene_match.push_back( pred );
			data_pos_bygene_match.push_back( pos );
			
			intersect_size++;
			gi_pred++;
			gi_pos++;
		}
	}
	
	gene_rc = double(intersect_size)/double(pos_size);
	gene_pr = double(intersect_size)/double(pred_size);
	
	if(intersect_size == 0){
		cerr << "ERROR: no common genes in pred and pos\n";
		exit(1);
	}
	if( (verbal && gene_rc<1.0)  ||  (gene_rc < 0.5) ){
		fprintf(stderr,"WARNING: Mismatch in gene labels?: recall of genes is less then 100%%: gene_rc=%0.3f\n",gene_rc);
	}
	if( (verbal && gene_pr<1.0)  ||  (gene_pr < 0.5) ){
		fprintf(stderr,"WARNING: Mismatch in gene labels?: precision of genes is less then 100%%: gene_pr=%0.3f\n",gene_pr);
	}
}


/*
 * Matches GO-terms in pred and pos data.
 * 
 * Saves matching data to data_pred_bygo_match and data_pos_bygo_match.
 * 
 * If GO-terms mismatch significantly: prints warning
 *
 * If all GO-terms mismatch: prints error and exits.
 *
 * PT modified this: This code excludes small GO classes from analysis
 */
void match_go_labels_filt_small(const vector<vector<go_annotation> > &data_pred_bygo,
			const vector<vector<go_annotation> > &data_pos_bygo,
			vector<vector<go_annotation> > &data_pred_bygo_match,
			vector<vector<go_annotation> > &data_pos_bygo_match,
			unsigned int SizeThreshold= 0,
			bool verbal=false){
			
	unsigned int pred_size= data_pred_bygo.size();
	unsigned int pos_size= data_pos_bygo.size();
	unsigned int intersect_size= 0;
	double go_rc= 0;
	double go_pr= 0;
	unsigned int gi_pred= 0; // GO indices
	unsigned int gi_pos= 0;
	
	data_pred_bygo_match.clear();
	data_pos_bygo_match.clear();
	
	while(gi_pred<pred_size && gi_pos<pos_size){
		if(data_pred_bygo[gi_pred][0].go < data_pos_bygo[gi_pos][0].go)
			gi_pred++;
		else if(data_pred_bygo[gi_pred][0].go > data_pos_bygo[gi_pos][0].go)
			gi_pos++;
		else{
			if( data_pos_bygo[gi_pos].size() > SizeThreshold){
				vector<go_annotation> pred = data_pred_bygo[gi_pred];
				vector<go_annotation> pos  = data_pos_bygo[gi_pos];
			
				data_pred_bygo_match.push_back( pred );
				data_pos_bygo_match.push_back( pos );
			}
			
			intersect_size++;
			gi_pred++;
			gi_pos++;
		}
	}
	
	go_rc = double(intersect_size)/double(pos_size);
	go_pr = double(intersect_size)/double(pred_size);
	
	if(intersect_size == 0){
		cerr << "ERROR: no common GO-terms in pred and true sets\n";
		exit(1);
	}
	if( (verbal && go_rc < 1.0) || (go_rc < 0.5) ){
		fprintf(stderr,"WARNING: recall of GO terms is less then 100%%: GO_rc=%0.3f\n",go_rc);
	}
	if( (verbal && go_pr < 1.0) || (go_pr < 0.5) ){
		fprintf(stderr,"WARNING: precision of GO terms is less then 100%%: GO_pr=%0.3f\n",go_pr);
	}	
}

	
/*
 * Returns a hash map of scores for each gene.
 * These can also be called gene-centric score, f.e. gene-centric ROC AUC values.
 * 
 * This function does NOT apply thresholds to the predicted data. Instead for each gene
 * the full list of predictions is included.
 *
 * data_pred		vector of go_annotations predicted by a classifier
 * data_pos		vector of go_annotations defining the positives
 * pmap			go -> list of parents
 * icmap		go -> information content
 *
 * scoref		function used to score similarity between two GO-lists
 *	AUC		ROC AUC
 * 	Jacc		Jaccard Index
 *	score_A		average of pairwise scores
 *	score_B		average of col maxima
 *	score_C		average of row maxima
 * 	score_D		average of score_B and score_C
 *	score_E		minimum of score_B and score_C
 *	score_F		averages of pooled col and row maxima
 *
 * simf			function used to score similarity between two GO-terms
 * go_count		number of uniq GO terms (used for AUC3)
 * 
 */
unordered_map<string,double> get_score_bygene(vector<go_annotation> &data_pred,
				vector<go_annotation> &data_pos,
				unordered_map<unsigned int,vector<unsigned int> > &pmap,
				unordered_map<unsigned int,double> &icmap,
				const string scoref,
				const string simf,
				unsigned int go_count_ontology = 0,
				bool verbal=false){
			
	if(scoref == "AUC3b" && go_count_ontology == 0){
		cerr << "ERROR: get_score_bygene(scoref=AUC3b): go_count_ontology==0\n";
		exit(1);
	}

	// sort gene>score
	sort(data_pred.begin(),data_pred.end(),compare_gene_score);
	sort(data_pos.begin(),data_pos.end(),compare_gene_score);
	
	// group by gene
	vector<vector<go_annotation> > data_pred_bygene_all = group_go_annotation_data(data_pred,"gene");
	vector<vector<go_annotation> > data_pos_bygene_all  = group_go_annotation_data(data_pos,"gene");
	
	// filter pred and pos data with matching gene-labels
	vector<vector<go_annotation> > data_pred_bygene;
	vector<vector<go_annotation> > data_pos_bygene;
	match_gene_labels(data_pred_bygene_all,data_pos_bygene_all, data_pred_bygene,data_pos_bygene,verbal);
	
	// number of uniq GOs in the database
	unordered_map<unsigned int,unsigned int> go_set_database;
	for(unsigned int i=0; i<data_pos.size(); i++)
		go_set_database[data_pos[i].go]= 1;
	for(unsigned int i=0; i<data_pred.size(); i++)
		go_set_database[data_pred[i].go]= 1;
	unsigned int go_count_database = go_set_database.size();
	
	
	// gene-wise scores
	unordered_map<string,double> score_map;
	
	for(unsigned int gi=0; gi<data_pred_bygene.size(); gi++){
	
		vector<unsigned int> pred(data_pred_bygene[gi].size(),0);
		vector<double> scores(data_pred_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pred_bygene[gi].size(); k++){
			pred[k]= data_pred_bygene[gi][k].go;
			scores[k]= data_pred_bygene[gi][k].score;
		}
				
		vector<unsigned int> pos(data_pos_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pos_bygene[gi].size(); k++)
			pos[k]=  data_pos_bygene[gi][k].go;
	
		// calling <stats.h> methods
		if(scoref == "U"){
			//score_map[data_pred_bygene[gi][0].gene] = get_U(pred,pos);
			score_map[data_pred_bygene[gi][0].gene] = get_U_with_ranks(pred,pos, scores);
		}
		else if(scoref == "AUC1"){
			//score_map[data_pred_bygene[gi][0].gene] = get_auc(pred,pos,1);
			score_map[data_pred_bygene[gi][0].gene] = get_U_with_ranks(pred,pos, scores, 1);
		}
		else if(scoref == "AUC2"){
			//score_map[data_pred_bygene[gi][0].gene] = get_auc(pred,pos,2);
			score_map[data_pred_bygene[gi][0].gene] = get_U_with_ranks(pred,pos, scores, 2);
		}
		else if(scoref == "AUC3a"){
			//score_map[data_pred_bygene[gi][0].gene] = get_auc(pred,pos,3,(go_count_database - pos.size()));
			score_map[data_pred_bygene[gi][0].gene] = get_U_with_ranks(pred,pos, scores, 3, (go_count_database - pos.size()));
		}
		else if(scoref == "AUC3b"){
			//score_map[data_pred_bygene[gi][0].gene] = get_auc(pred,pos,3,(go_count_ontology - pos.size()));
			score_map[data_pred_bygene[gi][0].gene] = get_U_with_ranks(pred,pos, scores, 3, (go_count_ontology - pos.size()));
		}
		else if(scoref == "AUCPR"){
			score_map[data_pred_bygene[gi][0].gene] = get_PR_AUC_ranks(pred,pos, scores);
		}		

		else if(scoref == "Jacc"){
			score_map[data_pred_bygene[gi][0].gene] = get_Jaccard(pred,pos);
		}
		else if(scoref == "wJacc"){
			score_map[data_pred_bygene[gi][0].gene] = get_WJaccard(pred,pos,icmap);	
		}	
		else if(scoref == "score_A" 
				|| scoref == "score_B" 
				|| scoref == "score_C" 
				|| scoref == "score_D"
				|| scoref == "score_E"
				|| scoref == "score_F"){
			vector<vector<double> > sim_mat= get_sim_matrix(pred,pos,pmap,icmap,simf);			
			if(scoref == "score_A")
				score_map[data_pred_bygene[gi][0].gene] = mean(as_vector(sim_mat));
			else if(scoref == "score_B")
				score_map[data_pred_bygene[gi][0].gene] = mean(col_max(sim_mat));
			else if(scoref == "score_C")
				score_map[data_pred_bygene[gi][0].gene] = mean(row_max(sim_mat));
			else if(scoref == "score_D"){
				vector<double> temp;
				temp.push_back( mean(col_max(sim_mat)));
				temp.push_back( mean(row_max(sim_mat)));
				score_map[data_pred_bygene[gi][0].gene] = mean(temp);
			}
			else if(scoref == "score_E"){
				vector<double> temp;
				temp.push_back( mean(col_max(sim_mat)));
				temp.push_back( mean(row_max(sim_mat)));
				score_map[data_pred_bygene[gi][0].gene] = min(temp);
			}
			else if(scoref == "score_F"){
				vector<double> tmp1= col_max(sim_mat);
				vector<double> tmp2= row_max(sim_mat);
				tmp1.insert(tmp1.end(),tmp2.begin(),tmp2.end());
				score_map[data_pred_bygene[gi][0].gene] = mean(tmp1);
			}
		}
		else{
			cerr << "ERROR: get_score_bygene: invalid SF:"<<scoref<<"\n";
			exit(1);
		}
	}	
	
	return score_map;
}

/*
 * Returns a hash map of scores for each gene.
 * These can also be called gene-centric scores, f.e. gene-centric ROC AUC values.
 *
 * Predictions are thresholded seperately for each gene in a ranked fashion:
 * for gene A:
 *     the best predictions
 *     the best and the second best
 *     ...
 *     all predictions
 *
 * For each gene scores across different thresholds are converted to a single value using sumf2.
 *
 * data_pred		vector of go_annotations predicted by a classifier
 * data_pos		vector of go_annotations defining the positives
 * pmap			go -> list of parents
 * icmap		go -> information content
 *
 * scoref		function used to score similarity between two GO-lists
 *	AUC		ROC AUC
 * 	Jacc		Jaccard Index
 *	score_A		average of pairwise scores
 *	score_B		average of col maxima
 *	score_C		average of row maxima
 * 	score_D		average of score_B and score_C
 *	score_E		minimum of score_B and score_C
 *	score_F		averages of pooled col and row maxima
 *
 * simf			function used to score similarity between two GO-terms
 *
 * sumf2	 	function used to convert scores for diff thresholds to a single value
 * 
 */
unordered_map<string,double> get_score_bygene_geneth(
				vector<go_annotation> &data_pred,
				vector<go_annotation> &data_pos,
				unordered_map<unsigned int,vector<unsigned int> > &pmap,
				unordered_map<unsigned int,double> &icmap,
				const string scoref,
				const string simf,
				const string sumf2){
	// sort gene>score
	sort(data_pred.begin(),data_pred.end(),compare_gene_score);
	sort(data_pos.begin(),data_pos.end(),compare_gene_score);
	
	// group by gene
	vector<vector<go_annotation> > data_pred_bygene_all = group_go_annotation_data(data_pred,"gene");
	vector<vector<go_annotation> > data_pos_bygene_all  = group_go_annotation_data(data_pos,"gene");
	
	// filter pred and pos data with matching gene-labels
	vector<vector<go_annotation> > data_pred_bygene;
	vector<vector<go_annotation> > data_pos_bygene;
	match_gene_labels(data_pred_bygene_all,data_pos_bygene_all, data_pred_bygene,data_pos_bygene);
	
	
	// gene-wise scores
	unordered_map<string,double> score_map;
	
	for(unsigned int gi=0; gi<data_pred_bygene.size(); gi++){
		vector<unsigned int> pred(data_pred_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pred_bygene[gi].size(); k++)
			pred[k]= data_pred_bygene[gi][k].go;
				
		vector<unsigned int> pos(data_pos_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pos_bygene[gi].size(); k++)
			pos[k]=  data_pos_bygene[gi][k].go;
				
		// precalculate SIM matrix
		vector<vector<double> > sim_mat;
		if(scoref == "score_A" 
				|| scoref == "score_B" 
				|| scoref == "score_C" 
				|| scoref == "score_D"
				|| scoref == "score_E"
				|| scoref == "score_F"){
			sim_mat = get_sim_matrix(pred,pos,pmap,icmap,simf);
		}	
			
		// scoring pred at stepwise thresholds (equal to scoring prefixes of increasing length)
		vector<double> score_list_th(data_pred_bygene[gi].size(),0);
			
		for(unsigned int ti=0; ti<score_list_th.size(); ti++){
				
			if(scoref == "Jacc"){
				vector<unsigned int> pred_pref= get_prefix(pred,ti+1);
				score_list_th[ti] = get_Jaccard(pred_pref,pos);
			}
			if(scoref == "wJacc"){
				vector<unsigned int> pred_pref= get_prefix(pred,ti+1);
				score_list_th[ti] = get_WJaccard(pred_pref,pos,icmap);
			}				
			else if(scoref == "score_A")
				score_list_th[ti] = mean(as_vector( get_prefix(sim_mat,ti+1) ) );
			else if(scoref == "score_B")
				score_list_th[ti] = mean(col_max( get_prefix(sim_mat,ti+1)) );
			else if(scoref == "score_C")
				score_list_th[ti] = mean(row_max( get_prefix(sim_mat,ti+1)) );
			else if(scoref == "score_D"){
				vector<double> temp;
				temp.push_back( mean(col_max( get_prefix(sim_mat,ti+1) )));
				temp.push_back( mean(row_max( get_prefix(sim_mat,ti+1) )));
				score_list_th[ti] = mean(temp);
			}
			else if(scoref == "score_E"){
				vector<double> temp;
				temp.push_back( mean(col_max( get_prefix(sim_mat,ti+1) )));
				temp.push_back( mean(row_max( get_prefix(sim_mat,ti+1) )));
				score_list_th[ti]= min(temp);
			}
			else if(scoref == "score_F"){
				vector<double> tmp1= col_max( get_prefix(sim_mat,ti+1) );
				vector<double> tmp2= row_max( get_prefix(sim_mat,ti+1) );
				tmp1.insert(tmp1.end(),tmp2.begin(),tmp2.end());
				score_list_th[ti]= mean(tmp1);
			}
			else{
				cerr << "ERROR: get_score_bygene: unsupported SF:"<<scoref<<"\n";
				exit(1);
			}
		}
		
		// summarizing over thresholds
		if(sumf2 == "max")
			score_map[data_pred_bygene[gi][0].gene] = max(score_list_th);
		else if(sumf2 == "mean")
			score_map[data_pred_bygene[gi][0].gene] = mean(score_list_th);
		else if(sumf2 == "median")
			score_map[data_pred_bygene[gi][0].gene] = median(score_list_th);
		else{
			cerr<< "ERROR: get_score_bygene: unsupported SUMF2: "<<sumf2<<"\n";
			exit(1);
		}
	}
	
	return score_map;
}

/*
 * Returns a list of scores, that are generated by iterating thresholds from th_min to th_max
 * in 100 steps (th_min/th_max are looked up from the predicted data). For each th looks for
 * genes that have predictions over th and score them using the specified scoref.
 * At each th converts a list of scores across genes to a single value using sumf1 (max/mean/median).
 *
 * data_pred		vector of go_annotations predicted by a classifier
 * data_pos		vector of go_annotations defining the positives
 * pmap			go -> list of parents
 * icmap		go -> information content
 *
 * scoref		function used to score similarity between two GO-lists
 *	AUC		ROC AUC
 * 	Jacc		Jaccard Index
 *	score_A		average of pairwise scores
 *	score_B		average of col maxima
 *	score_C		average of row maxima
 * 	score_D		average of score_B and score_C
 *	score_E		minimum of score_B and score_C
 *	score_F		averages of pooled col and row maxima
 *
 * simf			function used to score similarity between two GO-terms
 *
 * sumf1		function used to summarize score across genes
 * 
 */
vector<double> get_score_bygene_th(vector<go_annotation> &data_pred,
				vector<go_annotation> &data_pos,
				unordered_map<unsigned int,vector<unsigned int> > &pmap,
				unordered_map<unsigned int,double> &icmap,
				const string scoref,
				const string simf,
				const string sumf1,
				bool verbal=false){
	// sort gene>score
	sort(data_pred.begin(),data_pred.end(),compare_gene_score);
	sort(data_pos.begin(),data_pos.end(),compare_gene_score);
	
	// group by gene
	vector<vector<go_annotation> > data_pred_bygene_all = group_go_annotation_data(data_pred,"gene");
	vector<vector<go_annotation> > data_pos_bygene_all  = group_go_annotation_data(data_pos,"gene");
	
	// filter pred and pos data with matching gene-labels
	vector<vector<go_annotation> > data_pred_bygene;
	vector<vector<go_annotation> > data_pos_bygene;
	match_gene_labels(data_pred_bygene_all,data_pos_bygene_all, data_pred_bygene,data_pos_bygene,verbal);
		
	// RESTRUCTURE DATA: PRED AND POS VECTOR OF GO-LISTS
	vector<vector<unsigned int> > pred_golist_list(data_pred_bygene.size());
	vector<vector<unsigned int> > pos_golist_list(data_pos_bygene.size());
	for(unsigned int gi=0; gi< data_pred_bygene.size(); gi++){
		vector<unsigned int> golist(data_pred_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pred_bygene[gi].size(); k++)
			golist[k]= data_pred_bygene[gi][k].go;
		pred_golist_list[gi]= golist;
	}
	for(unsigned int gi=0; gi< data_pos_bygene.size(); gi++){
		vector<unsigned int> golist(data_pos_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pos_bygene[gi].size(); k++)
			golist[k]= data_pos_bygene[gi][k].go;
		pos_golist_list[gi]= golist;
	}
	
	// PRECALCULATE SIM MATRICES
	vector< vector<vector<double> > > simmat_list;
	if(scoref == "score_A" 
		|| scoref == "score_B" 
		|| scoref == "score_C" 
		|| scoref == "score_D"
		|| scoref == "score_E"
		|| scoref == "score_F"){
		for(unsigned int gi=0; gi< data_pos_bygene.size(); gi++){
			vector<vector<double> > simmat = get_sim_matrix(pred_golist_list[gi],pos_golist_list[gi],
									pmap,icmap,simf);
			simmat_list.push_back(simmat);
		}
	}
	
	// ITERATE THRESHOLDS
	vector<unsigned int> pred;
	vector<unsigned int> pos;	
	vector<double> sc_over_th;	
	vector<double> sc_over_genes;
	
	// THRESHOLD VALUES
	vector<double> score_values;
	for(unsigned int i=0; i<data_pred.size(); i++)
		score_values.push_back(data_pred[i].score); 
	double th_min = min(score_values);
	double th_max= max(score_values);
	double th_step= (th_max - th_min)/double(100);
		
	for(double th= th_min; th<=th_max; th+= th_step){
		
		sc_over_genes.clear();
		
		for(unsigned int gi=0; gi<data_pred_bygene.size(); gi++){
			unsigned int k = 0;
			while(k<data_pred_bygene[gi].size() && data_pred_bygene[gi][k].score >= th){
				k++;
			}
		
			if(k > 0){				
				pred= get_prefix(pred_golist_list[gi],k);
				pos = pos_golist_list[gi];	
				if(scoref == "Jacc"){
					sc_over_genes.push_back( get_Jaccard(pred,pos) );
				}
				else if(scoref == "wJacc"){
					sc_over_genes.push_back( get_WJaccard(pred,pos,icmap) );
				}					
				else if(scoref == "score_A")
					sc_over_genes.push_back( mean(as_vector( get_prefix(simmat_list[gi],k))) );
				else if(scoref == "score_B")
					sc_over_genes.push_back( mean(col_max( get_prefix(simmat_list[gi],k))) );
				else if(scoref == "score_C")
					sc_over_genes.push_back( mean(row_max( get_prefix(simmat_list[gi],k))) );
				else if(scoref == "score_D"){
					vector<double> temp;
					temp.push_back( mean(col_max( get_prefix(simmat_list[gi],k) )));
					temp.push_back( mean(row_max( get_prefix(simmat_list[gi],k) )));
					sc_over_genes.push_back( mean(temp) );
				}
				else if(scoref == "score_E"){
					vector<double> temp;
					temp.push_back( mean(col_max( get_prefix(simmat_list[gi],k) )));
					temp.push_back( mean(row_max( get_prefix(simmat_list[gi],k) )));
					sc_over_genes.push_back( min(temp) );
				}
				else if(scoref == "score_F"){
					vector<double> temp= col_max( get_prefix(simmat_list[gi],k) );
					vector<double> temp2= row_max( get_prefix(simmat_list[gi],k) );
					temp.insert(temp.end(),temp2.begin(),temp2.end());
					sc_over_genes.push_back( mean(temp) );
				}
				else{
					cerr << "ERROR: get_score_bygene_th: unsupported SF:"<<scoref<<"\n";
					exit(1);
				}					
			}
		}
		
		// CONVERT SCORE LIST ACROSS GENES TO A SIGNLE VALUE
		if(sumf1 == "max")
			sc_over_th.push_back( max(sc_over_genes) );
		else if(sumf1 == "mean")
			sc_over_th.push_back( mean(sc_over_genes) );
		else if(sumf1 == "median")
			sc_over_th.push_back( median(sc_over_genes) );
		else{
			cerr << "ERROR: get_score_bygene: unsupported SUMF1:"<<sumf1<<"\n";
			exit(1);
		}		
	}
	
	// 
	return sc_over_th;
}



/*
 * CAFA 2: Fmax score based on pr(t) and rc(t)
 *
 */
double get_score_cafa_fmax(vector<go_annotation> &data_pred,
				vector<go_annotation> &data_pos){
	// sort gene>score
	sort(data_pred.begin(),data_pred.end(),compare_gene_score);
	sort(data_pos.begin(),data_pos.end(),compare_gene_score);
	
	// group by gene
	vector<vector<go_annotation> > data_pred_bygene_all = group_go_annotation_data(data_pred,"gene");
	vector<vector<go_annotation> > data_pos_bygene_all  = group_go_annotation_data(data_pos,"gene");
	
	// filter pred and pos data with matching gene-labels
	vector<vector<go_annotation> > data_pred_bygene;
	vector<vector<go_annotation> > data_pos_bygene;
	match_gene_labels(data_pred_bygene_all,data_pos_bygene_all, data_pred_bygene,data_pos_bygene);
	
	// RESTRUCTURE DATA: PRED AND POS VECTOR OF GO-LISTS
	vector<vector<unsigned int> > pred_golist_list(data_pred_bygene.size());
	vector<vector<unsigned int> > pos_golist_list(data_pos_bygene.size());
	for(unsigned int gi=0; gi< data_pred_bygene.size(); gi++){
		vector<unsigned int> golist(data_pred_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pred_bygene[gi].size(); k++)
			golist[k]= data_pred_bygene[gi][k].go;
		pred_golist_list[gi]= golist;
	}
	for(unsigned int gi=0; gi< data_pos_bygene.size(); gi++){
		vector<unsigned int> golist(data_pos_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pos_bygene[gi].size(); k++)
			golist[k]= data_pos_bygene[gi][k].go;
		pos_golist_list[gi]= golist;
	}	
	
	// ITERATE THRESHOLDS
	vector<unsigned int> pred;
	vector<unsigned int> pos;
	vector<unsigned int> intersect;	
	vector<double> F_over_th;	
	vector<double> pr_over_genes;
	vector<double> rc_over_genes;
	double pr,rc,F;
	
	// THRESHOLD VALUES
	vector<double> score_values;
	for(unsigned int i=0; i<data_pred.size(); i++)
		score_values.push_back(data_pred[i].score); 
	double th_min = min(score_values);
	double th_max= max(score_values);
	double th_step= (th_max - th_min)/double(100);
	
		
	for(double th= th_min; th<=th_max; th+= th_step){
		pr_over_genes.clear();
		rc_over_genes.clear();
		
		for(unsigned int gi=0; gi<data_pred_bygene.size(); gi++){
			unsigned int k = 0;
			while( k<data_pred_bygene[gi].size() && data_pred_bygene[gi][k].score >= th){
				k++;
			}
			if(k > 0){				
				pred= get_prefix(pred_golist_list[gi],k);
				pos = pos_golist_list[gi];
				intersect= set_intersection(pred,pos);
					
				pr_over_genes.push_back( double(intersect.size())/double(pred.size()) );
				rc_over_genes.push_back( double(intersect.size())/double(pos.size()) );
			}
		}
		
		// NOTE THAT TAKING MEAN SLIGTLY DIFFERS FROM CAFA I DEFINITION
		// pr= mean(pr_over_genes);
		// rc= mean(rc_over_genes);
		pr= sum(pr_over_genes)/double(pr_over_genes.size());
		rc= sum(rc_over_genes)/double(data_pred_bygene.size());
		F = (2.0*pr*rc) / (pr + rc);
		F_over_th.push_back( F );		
	}
	
	// Fmax VALUE OVER THRESHOLDS
	return max(F_over_th);
}


/*
 * CAFA 2: Smin score based on ru(t) and mi(t)
 *
 * Assumes that icmap values are defined for all GO-ids in data_pred and data_pos.
 *
 */
double get_score_cafa_smin(vector<go_annotation> &data_pred,
				vector<go_annotation> &data_pos,
				unordered_map<unsigned int,double> &icmap){
	// sort gene>score
	sort(data_pred.begin(),data_pred.end(),compare_gene_score);
	sort(data_pos.begin(),data_pos.end(),compare_gene_score);
	
	// group by gene
	vector<vector<go_annotation> > data_pred_bygene_all = group_go_annotation_data(data_pred,"gene");
	vector<vector<go_annotation> > data_pos_bygene_all  = group_go_annotation_data(data_pos,"gene");
	
	// filter pred and pos data with matching gene-labels
	vector<vector<go_annotation> > data_pred_bygene;
	vector<vector<go_annotation> > data_pos_bygene;
	match_gene_labels(data_pred_bygene_all,data_pos_bygene_all, data_pred_bygene,data_pos_bygene);
	
	// RESTRUCTURE DATA: PRED AND POS VECTOR OF GO-LITS
	vector<vector<unsigned int> > pred_golist_list(data_pred_bygene.size());
	vector<vector<unsigned int> > pos_golist_list(data_pos_bygene.size());
	for(unsigned int gi=0; gi< data_pred_bygene.size(); gi++){
		vector<unsigned int> golist(data_pred_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pred_bygene[gi].size(); k++)
			golist[k]= data_pred_bygene[gi][k].go;
		pred_golist_list[gi]= golist;
	}
	for(unsigned int gi=0; gi< data_pos_bygene.size(); gi++){
		vector<unsigned int> golist(data_pos_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pos_bygene[gi].size(); k++)
			golist[k]= data_pos_bygene[gi][k].go;
		pos_golist_list[gi]= golist;
	}	
	
	// ITERATE THRESHOLDS
	vector<unsigned int> pred;
	vector<unsigned int> pos;
	vector<unsigned int> ru_set;
	vector<unsigned int> mi_set;	
	vector<double> S_over_th;	
	vector<double> ru_over_genes;
	vector<double> mi_over_genes;
	double ru,mi,S;
	
	// threshold values
	vector<double> score_values;
	for(unsigned int i=0; i<data_pred.size(); i++)
		score_values.push_back(data_pred[i].score); 
	double th_min = min(score_values);
	double th_max= max(score_values);
	double th_step= (th_max - th_min)/double(100);
		
	for(double th= th_min; th<=th_max; th+= th_step){	
		ru_over_genes.clear();
		mi_over_genes.clear();
		
		for(unsigned int gi=0; gi<data_pred_bygene.size(); gi++){
			unsigned int k = 0;
			while( k<data_pred_bygene[gi].size() && data_pred_bygene[gi][k].score >= th){
				k++;
			}
			if(k > 0){				
				pred= get_prefix(pred_golist_list[gi],k);
				pos = pos_golist_list[gi];
					
				ru_set= set_difference(pos,pred);
				mi_set= set_difference(pred,pos);
					
				ru = 0;
				mi = 0;
				for(unsigned int i=0; i<ru_set.size(); i++)
					ru+= icmap[ru_set[i]];
				for(unsigned int i=0; i<mi_set.size(); i++)
					mi+= icmap[mi_set[i]];
				
				ru_over_genes.push_back( ru );
				mi_over_genes.push_back( mi );
			}
		}
		
		// MEAN VALUES ACROSS GENES
		ru= mean(ru_over_genes);
		mi= mean(mi_over_genes);
		S = sqrt( ru*ru + mi*mi );
		S_over_th.push_back( S );		
	}
	
	// Smin VALUE OVER THRESHOLDS
	return min(S_over_th);
}

/* Petri modification. 
 * Testing different ordering of summing in Smin calculus (Smin per gene first; calculate average of these scores)
 * Testing Smin calculus with power of 1 (not power of 2)		
 * Added weighted Jaccard (SimGIC) score also. Because it is easily obtained with data we already have
 *
 * Minimum of each Smin scores is selected as output over thresholds.
 * Maximum of SimGIC is selected as output over thresholds. 
 */    
vector<double> get_score_cafa_smin_PT(vector<go_annotation> &data_pred,
				vector<go_annotation> &data_pos,
				unordered_map<unsigned int,double> &icmap){
	// sort gene>score
	sort(data_pred.begin(),data_pred.end(),compare_gene_score);
	sort(data_pos.begin(),data_pos.end(),compare_gene_score);
	
	// group by gene
	vector<vector<go_annotation> > data_pred_bygene_all = group_go_annotation_data(data_pred,"gene");
	vector<vector<go_annotation> > data_pos_bygene_all  = group_go_annotation_data(data_pos,"gene");
	
	// filter pred and pos data with matching gene-labels
	vector<vector<go_annotation> > data_pred_bygene;
	vector<vector<go_annotation> > data_pos_bygene;
	match_gene_labels(data_pred_bygene_all,data_pos_bygene_all, data_pred_bygene,data_pos_bygene);
	
	// RESTRUCTURE DATA: PRED AND POS VECTOR OF GO-LITS
	vector<vector<unsigned int> > pred_golist_list(data_pred_bygene.size());
	vector<vector<unsigned int> > pos_golist_list(data_pos_bygene.size());
	for(unsigned int gi=0; gi< data_pred_bygene.size(); gi++){
		vector<unsigned int> golist(data_pred_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pred_bygene[gi].size(); k++)
			golist[k]= data_pred_bygene[gi][k].go;
		pred_golist_list[gi]= golist;
	}
	for(unsigned int gi=0; gi< data_pos_bygene.size(); gi++){
		vector<unsigned int> golist(data_pos_bygene[gi].size(),0);
		for(unsigned int k=0; k<data_pos_bygene[gi].size(); k++)
			golist[k]= data_pos_bygene[gi][k].go;
		pos_golist_list[gi]= golist;
	}	
	
	// ITERATE THRESHOLDS
	vector<unsigned int> pred;
	vector<unsigned int> pos;
	vector<unsigned int> ru_set;
	vector<unsigned int> mi_set;	
	vector<unsigned int> inters_set;	

	vector<double> S_over_th;  // these store data over threshold positions
	vector<double> S2_over_th;
	vector<double> S3_over_th;
	vector<double> SimUI_over_th;
	vector<double> SimGIC_over_th;
	vector<double> SimGIC2_over_th;	
	vector<double> ru_over_genes;
	vector<double> mi_over_genes;
	vector<double> Smin_per_gene;
	vector<double> SimUI_per_gene;
	vector<double> SimGIC_per_gene;
	double ru,mi,S,S2, S3, SimUI, SimGIC, tmp, inters;  // PT added here many variables
	double SimGIC2, inters_sum;
	
	// threshold values
	vector<double> score_values;
	for(unsigned int i=0; i<data_pred.size(); i++)
		score_values.push_back(data_pred[i].score); 
	double th_min = min(score_values);
	double th_max= max(score_values);
	double th_step= (th_max - th_min)/double(100);
		
	// Loop over threshold values
	for(double th= th_min; th<=th_max; th+= th_step){	
		ru_over_genes.clear();
		mi_over_genes.clear();
		Smin_per_gene.clear();
		SimGIC_per_gene.clear();
		SimUI_per_gene.clear();
		inters_sum = 0;
		// loop over genes
		for(unsigned int gi=0; gi<data_pred_bygene.size(); gi++){
			unsigned int k = 0;
			while( k<data_pred_bygene[gi].size() && data_pred_bygene[gi][k].score >= th){
				k++;
			}
			if(k > 0){				
				pred= get_prefix(pred_golist_list[gi],k);
				pos = pos_golist_list[gi];
					
				ru_set= set_difference(pos,pred);
				mi_set= set_difference(pred,pos);
				inters_set= set_intersection(pred,pos);
					
				ru = 0;
				mi = 0;
				inters = 0;
				for(unsigned int i=0; i<ru_set.size(); i++)
					ru+= icmap[ru_set[i]];
				for(unsigned int i=0; i<mi_set.size(); i++)
					mi+= icmap[mi_set[i]];
				for(unsigned int i=0; i<inters_set.size(); i++)
					inters+= icmap[inters_set[i]];
				
				ru_over_genes.push_back( ru );
				mi_over_genes.push_back( mi );
				
				// PT adds calculations
				inters_sum = inters_sum + inters;
				//reordered Smin score below
				tmp = sqrt( ru*ru + mi*mi );
                                Smin_per_gene.push_back( tmp );	
				if( inters_set.size() > 0){
                                   tmp = (float) inters_set.size() / (ru_set.size() + mi_set.size() + inters_set.size() );
				}
				else{  
				   tmp = 0;
				}
				SimUI_per_gene.push_back( tmp );
				// Next is IC weighted Jaccard
				// Code assumes (ru + mi + inters) > 0			
				if( inters > 0){
				    tmp = (float) inters / (ru + mi + inters );
                                }
				else{  
				   tmp = 0;
				}
				SimGIC_per_gene.push_back( tmp );
			}
		}
		
		// MEAN VALUES ACROSS GENES
		ru= mean(ru_over_genes);
		mi= mean(mi_over_genes);

		S = sqrt( ru*ru + mi*mi );
		S_over_th.push_back( S );	
		// PT adds calculations
		// S2 and S3 are alternative ways to calc Smin score
		S2 = ru + mi;
		S2_over_th.push_back( S2 );
		S3 = mean(Smin_per_gene);
                S3_over_th.push_back( S3 );
		// Unweighted Jaccard (SimUI)
		SimUI = mean(SimUI_per_gene);
		SimUI_over_th.push_back( SimUI );
		// Weighted Jaccard (SimGIC)
		SimGIC = mean(SimGIC_per_gene);
		SimGIC_over_th.push_back( SimGIC );
		// Weighted Jaccard - Smin hybrid
		if(  inters_sum > 0){
			SimGIC2 = (float) inters_sum / (inters_sum + sum(ru_over_genes) + sum(mi_over_genes) )	;
		}
		else{	
			SimGIC2 = 0;	
		}
		SimGIC2_over_th.push_back( SimGIC2 );
	}

	// These would go to output	
	double S_out = min(S_over_th);
	double S2_out = min(S2_over_th);
	double S3_out = min(S3_over_th);
	double SimGIC_out = max(SimGIC_over_th);
	double SimGIC2_out = max(SimGIC2_over_th);
	double SimUI_out = max(SimUI_over_th);
        vector <double> output(6, 0);
	output[0] = S_out;
	output[1] = S2_out;
	output[2] = S3_out;
	output[3] = SimUI_out;
	output[4] = SimGIC_out;
	output[5] = SimGIC2_out;
	return output;
}


/*
 * Returns a vector of scores for each go-term.
 * These are also referred to as term-centric scores.
 * 
 * data_pred		vector of go_annotations as predicted by a classifier
 * data_pos		vector of go_annotations defining the positives
 * scoref		scoring function
 * MinClass_Sz		Minimum size for GO class in Positive Set. If this parameter is set to zero then all classes are taken 
 * 
 */
vector<double> get_score_bygo(vector<go_annotation> &data_pred, 
				vector<go_annotation> &data_pos,
				const string scoref,
				const string sumf="NULL",
				unsigned int MinClass_Sz = 10,
				bool verbal=false){

	// sort by go>scrore
	sort(data_pred.begin(),data_pred.end(),compare_go_score);
	sort(data_pos.begin(),data_pos.end(),compare_go_score);
	
	// group by go
	vector<vector<go_annotation> > data_pred_bygo_all = group_go_annotation_data(data_pred,"go");
	vector<vector<go_annotation> > data_pos_bygo_all = group_go_annotation_data(data_pos,"go");
	
	// filter pred and pos data with matching go-labels
	vector<vector<go_annotation> > data_pred_bygo;
	vector<vector<go_annotation> > data_pos_bygo;
	//match_go_labels(data_pred_bygo_all,data_pos_bygo_all, data_pred_bygo,data_pos_bygo, verbal);
	match_go_labels_filt_small(data_pred_bygo_all,data_pos_bygo_all, data_pred_bygo,data_pos_bygo, MinClass_Sz,verbal);

	// map gene names to unique ints
	unordered_map<string, unsigned int> name2id;
	unsigned int id=0;
	for(unsigned int i=0; i<data_pos.size(); i++){
		if(name2id.count(data_pos[i].gene) < 1){
			name2id[data_pos[i].gene]= id;
			id++;
		}
	}
	for(unsigned int i=0; i<data_pred.size(); i++){
		if(name2id.count(data_pred[i].gene) < 1){
			name2id[data_pred[i].gene]= id;
			id++;
		}
	}


	// Scores
	
	vector <double> S(data_pred_bygo.size(),0);
	
	for(unsigned int i=0; i<data_pred_bygo.size(); i++){ // iterating terms
		
		if(data_pred_bygo[i][0].go != data_pos_bygo[i][0].go){
			cerr << "ERROR: get_score_bygo: mismatch of go terms in pred and pos:\n";
			cerr << "pred.go="<< data_pred_bygo[i][0].go <<"\n";
			cerr << "pos.go= "<< data_pos_bygo[i][0].go <<"\n";
			exit(1);}
	
		vector<unsigned int> pred(data_pred_bygo[i].size(),0);
		vector<unsigned int> pos(data_pos_bygo[i].size(),0);
		vector<double> scores(data_pred_bygo[i].size(),0);
		for(unsigned int j=0; j<data_pred_bygo[i].size(); j++){
			pred[j]= name2id[data_pred_bygo[i][j].gene];
			scores[j]= data_pred_bygo[i][j].score;
		}
		for(unsigned int j=0; j<data_pos_bygo[i].size(); j++)
			pos[j]= name2id[data_pos_bygo[i][j].gene];
	
	
		// calling <stats.h> methods
		if(scoref=="U"){
			S[i]= get_U_with_ranks(pred,pos, scores);
			//S[i]= get_U_with_ranks(pred,pos)
		}
		else if(scoref=="AUC1"){
			//S[i]= get_auc(pred,pos,1);
			S[i]= get_U_with_ranks(pred,pos, scores, 1);
		}
		else if(scoref=="AUC2"){
			//S[i]= get_auc(pred,pos,2);
			S[i]= get_U_with_ranks(pred,pos, scores, 2);
		}
		else if(scoref=="AUC3a"){
			//S[i]= get_auc(pred,pos,3,(name2id.size() - pos.size()));
			if( name2id.size() == pos.size() ){
				cerr << " get_score_bygo N-size 0\n";
				cerr << pos.size() << " \n" ;
			}
			S[i]= get_U_with_ranks(pred,pos, scores, 3, (name2id.size() - pos.size()));
		}
		else if(scoref=="AUCPR"){
			//S[i]= get_auc(pred,pos,3,(name2id.size() - pos.size()));
			S[i]= get_PR_AUC_ranks(pred,pos, scores);
		}
		else if(scoref=="Jacc"){
			S[i]= get_Jaccard(pred,pos);
		}		
		else{
			cerr << "ERROR: get_score_bygo: invalid SF:"<<scoref<<"\n";
			exit(1);
		}
		if(  S[i] < 0 || NotFinite(S[i])){
			cerr << "ERROR: get_score_bygo: error with output value="<<S[i]<<" SF="<<scoref<<"\n";
			exit(1);
		
		} 
	}
	
	return S; 
}

double get_score_singlelist(vector<go_annotation> &data_pred,
				vector<go_annotation> &data_pos,
				const string scoref,
				const string thf,
				const string sumf,
				unsigned int go_count_ontology = 0){
				
	if(scoref == "AUC3b" && go_count_ontology == 0){
		cerr << "ERROR: get_score_singlelist: go_count_ontology==0\n";
		exit(1);
	}
	
	// sort data by score
	sort(data_pred.begin(),data_pred.end(),compare_score);
	sort(data_pos.begin(),data_pos.end(),compare_score);
	
	// map gene names to unique ints
	unordered_map<string, unsigned int> name2id;
	unsigned int id=0;
	for(unsigned int i=0; i<data_pos.size(); i++){
		if(name2id.count(data_pos[i].gene) < 1){
			name2id[data_pos[i].gene]= id;
			id++;
		}
	}
	for(unsigned int i=0; i<data_pred.size(); i++){
		if(name2id.count(data_pred[i].gene) < 1){
			name2id[data_pred[i].gene]= id;
			id++;
		}
	}
	
	// number of uniq genes and GOs in the database
	unordered_map<unsigned int,unsigned int> go_set_database;
	for(unsigned int i=0; i<data_pos.size(); i++)
		go_set_database[data_pos[i].go]= 1;
	for(unsigned int i=0; i<data_pred.size(); i++)
		go_set_database[data_pred[i].go]= 1;
	unsigned int go_count_database = go_set_database.size();
	unsigned int gene_count_database= name2id.size();	

	
	// collect <geneid,go> tuples and convert them to unique long int keys
	vector<long long unsigned int> keys_pred;
	vector<double> scores;
	for(unsigned int i=0; i<data_pred.size(); i++){
		keys_pred.push_back( encode_key( name2id[data_pred[i].gene],data_pred[i].go));
		scores.push_back( data_pred[i].score );
	}
	
	vector<long long unsigned int> keys_pos;
	for(unsigned int i=0; i<data_pos.size(); i++){
		keys_pos.push_back( encode_key( name2id[data_pos[i].gene],data_pos[i].go));
	}
		
	
	// calling <stats.h> methods
	if(scoref=="U"){
		//return get_U(keys_pred,keys_pos);
		return get_U_with_ranks(keys_pred,keys_pos, scores);
	}
	else if(scoref=="AUC1"){
		//return(get_auc(keys_pred,keys_pos,1));
		return get_U_with_ranks(keys_pred,keys_pos, scores, 1);
	}
	else if(scoref == "AUC2"){
		//return(get_auc(keys_pred,keys_pos,2));
		return get_U_with_ranks(keys_pred,keys_pos, scores, 2);
	}
	else if(scoref == "AUC3a"){
		//return get_auc(keys_pred,keys_pos,3, (gene_count_database*go_count_database - keys_pos.size()));
		return get_U_with_ranks(keys_pred,keys_pos, scores, 3, (gene_count_database*go_count_database - keys_pos.size()));
	}
	else if(scoref == "AUC3b"){
		//return get_auc(keys_pred,keys_pos,3, (gene_count_database*go_count_ontology - keys_pos.size()));
		return get_U_with_ranks(keys_pred,keys_pos, scores, 3, (gene_count_database*go_count_ontology - keys_pos.size()));
	}
	else if(scoref == "AUCPR"){
		return get_PR_AUC_ranks(keys_pred,keys_pos, scores);
	}	
	else if(scoref=="Jacc" && (thf=="all" || thf=="NULL")){
		return(get_Jaccard(keys_pred,keys_pos));
	}
	else if(scoref=="Jacc" && thf=="th"){
		// THRESHOLD VALUES
		vector<double> score_values;
		for(unsigned int i=0; i<data_pred.size(); i++)
			score_values.push_back(data_pred[i].score); 
		double th_min = min(score_values);
		double th_max= max(score_values);
		double th_step= (th_max - th_min)/double(100);
		
		vector<long long unsigned int> keys_pred_th;
		vector<double> sc_over_th;
		
		for(double th= th_min; th<=th_max; th+= th_step){
		
			keys_pred_th.clear();
			unsigned int i=0;
			while(i<data_pred.size() && data_pred[i].score>=th){
				keys_pred_th.push_back(keys_pred[i]);
				i++;
			}
			
			sc_over_th.push_back( get_Jaccard(keys_pred_th,keys_pos));		
		}
		
		if(sumf == "max")
			return( max(sc_over_th) );
		else if(sumf == "mean")
			return( mean(sc_over_th) );
		else if(sumf == "median")
			return( median(sc_over_th) );
		else{
			cerr << "ERROR: get_score_singlelist: unsupported SUMF:"<<sumf<<"\n";
			exit(1);
		}
	}	
	else{
		cerr << "ERROR: get_score_singlelist: invalid SCOREF:"<<scoref<<"\n";
		exit(1);
	}
}


 
 
 

