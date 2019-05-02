#include <iostream>
#include <vector>
#include <set>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <float.h>
#include <getopt.h>

// Local code
#include <utils.h>
#include <gograph.h>


using namespace std;


/*
 * QUICK INTERFACE TO <gograph.h> FOR READING PARENTS/CHILDREN/NEIGHBOURS FROM GO OBO 1.2 FILE
 *
 */
 



void print_usage(const char* name){
	cerr <<    "# gorelatives v1.2: listing neighbourhood for a given node in OBO 1.2 graph\n"
		<< "# copyright (C) 2017 Ilja Pljusnin, University of Helsinki\n"
		<< "\n"
		<< "USAGE: "<< name <<" -b obo_file -q qid/qlist -r isa[,partof,regulates,altid,consider,replacedby] -d parents|children|relatives [-k max_num] [-l] [-v]\n"
		<<"\n"
		<<"OPTIONS:\n"
		<<"-b\t: OBO 1.2 file\n"
		<<"-q\t: quiry id (without GO-qualifier)\n"
		<<"  \t: or a list of query ids in a file (one per line)\n"
		<<"  \t: or literal \"all\" to print relatives for all nodes in the onthology\n"
		<<"-r\t: list of literals that specify relations/edges to search through. Valid literals:\n"
		<<"  \t  isa|partof|regulates|altid|consider|replacedby\n"
		<<"-d\t: literal that specifies the direction of the search:\n"
		<<"  \t  parents|children|relatives\n"
		<<"-k\t: maximum number of relatives retrieved\n"
		<<"-l\t: print namespace and name for each query node\n"
		<<"-v\t: verbal mode: e.g. print warning if no relatives found for a given GO id\n"
		<<"\n";
}

int main(int argc, char** argv) {

	// paramenters
	const char* prog_name= argv[0];
	const char* obo_file= NULL;
	const char* qid_file= NULL;
	string relations= "null";
	string direction= "null";
	unsigned int max_num= 10000;
	unsigned int qid= 0;
	bool print_label= false;
	bool verbal= false;
	
	if(argc < 2){
		print_usage(prog_name);
		exit(1);
	}
	
	int option= 0;
	while ((option = getopt(argc, argv,"b:q:r:d:k:lv")) != -1) {
        switch (option) {
		case 'b':
			obo_file= optarg;
			break;
		case 'q':
			qid= atoi(optarg);
			if(qid==0)
				qid_file= optarg;
			break;
		case 'r':
			relations= string(optarg);
			break;
		case 'd':
			direction= string(optarg);
			break;
		case 'k':
			max_num= atoi(optarg);
			break;
		case 'l':
			print_label= true;
			break;
		case 'v':
			verbal= true;
			break;	
		
		case '?':
        		cerr << "option -"<<optopt<<" requires an argument\n\n";
        		exit(1);
             	default:
	     		print_usage(prog_name); 
                	exit(1);
        	}
    	}
	
	// sanity check
	if(obo_file== NULL){
		cerr<<"ERROR: invalid option: -b obo_file\n\n";
		print_usage(prog_name);
		exit(1);
	}
	if(qid==0 && qid_file==NULL){
		cerr<<"ERROR: invalid option: -q query_id/query_id_list.txt/all\n\n";
		print_usage(prog_name);
		exit(1);
	}
	if(relations == "null"){
		cerr<<"ERROR: invalid option: -r isa[,partof,regulates,altid,consider,replacedby]\n\n";
		print_usage(prog_name);
		exit(1);
	}
	if(direction == "null"){
		cerr<<"ERROR: invalid option: -d parents|children\n\n";
		print_usage(prog_name);
		exit(1);
	}
	
	
	//cout << "# reading obo file..\n"<<flush;
	vector<go_node> go_nodes= read_obo_file(obo_file); 
	
	//cout << "# creating graph..\n"<<flush;
	go_graph graph= create_go_graph(go_nodes);
	
	vector<unsigned int> qlist;
	if(qid>0){
		qlist.push_back(qid);
		
	}
	else if(qid_file != NULL && string(qid_file) == "all"){
		for(unsigned int i=0; i<go_nodes.size(); i++){
			qlist.push_back(go_nodes[i].id);
			// if node has alt-ids, we add those as well
			if( go_nodes[i].altid_list.size()> 0){
				qlist.insert(qlist.end(),go_nodes[i].altid_list.begin(),go_nodes[i].altid_list.end());
			}
		}
	}
	else if( qid_file != NULL){
	
		ifstream in(qid_file, ios::in);
		if(!in.is_open()){
			cerr << "ERROR: failed to open \'"<<qid_file<<"\'\n";
			exit(1);
		}
		char buffer[1000];
		while(in.getline(buffer,1000)){
			qid= atoi(buffer);
			if(qid>0)
				qlist.push_back(qid);
		}		
	}
	else{
		cerr << "ERROR: invalid query id (-q)\n";
		exit(1);
	}
	
		
	vector<string> rel_list= split(relations,",");
	bool include_altids= false;
	for(unsigned int i=0; i<rel_list.size(); i++){
		if(rel_list[i] == "altid")
			include_altids= true;
	}
	
	//cout << "# get_relatives..\n"<<flush;
	for(unsigned int i=0; i<qlist.size(); i++){

		unsigned int qnode_id= qlist[i];	// parameter qid or it's alt-id

		//mapping alt-ids: 
		// if our query id is and alt-id, then we replace it with the node, for which it is an alt-id
		// direct pointer for such node is in graph.altid_parents
		if(include_altids){
			if( graph.altid_parents.count(qlist[i])>0 ){
				qnode_id= graph.altid_parents[qlist[i]][0];
				
			}
		}

		vector<unsigned int> relatives = get_relatives(graph,qnode_id,relations,direction,max_num,verbal);
		
		
		// print results
		if(print_label){
			if(graph.nodes.count(qnode_id)>0){
				go_node qnode= graph.nodes[qnode_id];
				printf("%07u\t%s\t%s\t",qlist[i],qnode.name_space.c_str(),qnode.name.c_str());
			}
			else{
				printf("%07u\t\t\t",qlist[i]);
			}
		}
		else{
			printf("%07u\t",qlist[i]);	
		}
		
		if(relatives.size()>0)
			printf("%07u",relatives[0]);
			
		for(unsigned int j=1; j<relatives.size(); j++)
			printf(",%07u",relatives[j]);
		printf("\n");	
	}
}





