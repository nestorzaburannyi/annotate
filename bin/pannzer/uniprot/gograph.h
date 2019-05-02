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

using namespace std;

/*
 * SET OF METHODS FOR MANIPULATING GENE ONTOLOGY GRAPH
 *
 */

struct go_node{
	unsigned int id;
	string name;
	string name_space; // Note that "namespace" is a reserved word in C++
	bool is_obsolete;
	vector<unsigned int> isa_list;
	vector<unsigned int> partof_list;
	vector<unsigned int> regulates_list;
	vector<unsigned int> replacedby_list;
	vector<unsigned int> consider_list;
	vector<unsigned int> altid_list;
};

void print(const go_node node){
	printf("id\t%07u\n",node.id);
	cout <<"name\t"<<node.name<<"\n";
	cout <<"namespace\t"<<node.name_space<<"\n"; 
	printf("is_obsolete\t%s\n",(node.is_obsolete?"yes":"no"));
	if(node.isa_list.size()==0){
		printf("is_a\t[]\n");
	}else{
		printf("is_a\t[%07u",node.isa_list[0]);
		for(unsigned int i=1; i<node.isa_list.size(); i++) 
			printf(",%07u",node.isa_list[i]);
		printf("]\n");
	}
	if(node.partof_list.size()==0){
		printf("part_of\t[]\n");
	}else{
		printf("part_of\t[%07u",node.partof_list[0]);
		for(unsigned int i=1; i<node.partof_list.size(); i++) 
			printf(",%07u",node.partof_list[i]);
		printf("]\n");
	}
	if(node.regulates_list.size()==0){
		printf("regulates\t[]\n");
	}else{
		printf("regulates\t[%07u",node.regulates_list[0]);
		for(unsigned int i=1; i<node.regulates_list.size(); i++) 
			printf(",%07u",node.regulates_list[i]);
		printf("]\n");
	}
	if(node.replacedby_list.size()==0){
		printf("replacedby\t[]\n");
	}else{
		printf("replacedby\t[%07u",node.replacedby_list[0]);
		for(unsigned int i=1; i<node.replacedby_list.size(); i++) 
			printf(",%07u",node.replacedby_list[i]);
		printf("]\n");
	}
	if(node.consider_list.size()==0){
		printf("consider\t[]\n");
	}else{
		printf("consider\t[%07u",node.consider_list[0]);
		for(unsigned int i=1; i<node.consider_list.size(); i++) 
			printf(",%07u",node.consider_list[i]);
		printf("]\n");
	}
	if(node.altid_list.size()==0){
		printf("alt_id\t[]\n");
	}else{
		printf("alt_id\t[%07u",node.altid_list[0]);
		for(unsigned int i=1; i<node.altid_list.size(); i++)
			printf(",%07u",node.altid_list[i]);
		printf("]\n");
	}
}


/* 
 * Reads Gene Ontology node definitions from OBO file format.
 *
 * Based on strtok()
 *
 * fin		Input file in OBO 1.2 format.
 *
 * returns	vector of go_nodes
 *
 * NOTE: not all fields are parsed
 * 
 */

vector<go_node> read_obo_file(	const char* fin, 
				bool include_obsolete = true){
	ifstream in(fin, ios::in);
	if(!in.is_open()){
		cerr << "read_obo_file(): failed to open file \'"<<fin<<"\': exiting\n";
		exit(1);
	}
	char buffer[10000];
	unsigned int id = 0;
	string name = "null";
	string name_space = "null";
	vector<unsigned int> isa_list;
	vector<unsigned int> partof_list;
	vector<unsigned int> regulates_list;
	vector<unsigned int> replacedby_list;
	vector<unsigned int> consider_list;
	vector<unsigned int> altid_list;
	bool is_obsolete = false;
	
	int line_num=0;
	char* pch;
	
	vector<go_node> go_nodes;
	
	while(in.getline(buffer,10000)){
	
		line_num++;
		//cout<< "line "<<line_num<<"\n";
		
		string line(buffer);
		if(line == ""){
			if(id>0){
			vector<unsigned int> isa_list_copy= isa_list;
			vector<unsigned int> partof_list_copy= partof_list;
			vector<unsigned int> regulates_list_copy= regulates_list;
			vector<unsigned int> replacedby_list_copy= replacedby_list;
			vector<unsigned int> consider_list_copy= consider_list;
			vector<unsigned int> altid_list_copy= altid_list;
			
			// sorting
			sort(isa_list_copy.begin(),isa_list_copy.end());
			sort(partof_list_copy.begin(),partof_list_copy.end());
			sort(regulates_list_copy.begin(),regulates_list_copy.end());
			sort(replacedby_list_copy.begin(),replacedby_list_copy.end());
			sort(consider_list_copy.begin(),consider_list_copy.end());
			sort(altid_list_copy.begin(),altid_list_copy.end());
			
			/* removing duplicates
			iterator it= unique(isa_list_copy.begin(),isa_list_copy.end());
			isa_list_copy.resize( distance(isa_list_copy.begin(),it) );
			*/
			go_node temp= {id,
					name,
					name_space,
					is_obsolete,
					isa_list_copy,
					partof_list_copy,
					regulates_list_copy,
					replacedby_list_copy,
					consider_list_copy,
					altid_list_copy};
			if(!is_obsolete || include_obsolete)
				go_nodes.push_back(temp);
			
			id = 0;
			name = "null";
			name_space = "null";
			is_obsolete= false;
			isa_list.clear();
			partof_list.clear();
			regulates_list.clear();
			replacedby_list.clear();
			consider_list.clear();
			altid_list.clear();
			continue;
			}
		}
		else{
			pch= strtok(buffer,": ");
			
			if(strcmp(pch,"id")==0){
				pch= strtok(NULL," GO:");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				id= atoi(pch);
			}
			else if(strcmp(pch,"name")==0){
				pch= strtok(NULL,":");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				name= string(pch);
			}
			else if(strcmp(pch,"namespace")==0){
				pch= strtok(NULL,": ");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				name_space= string(pch);
			}
			else if(strcmp(pch,"is_a")==0){
				pch= strtok(NULL," GO:");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				isa_list.push_back( atoi(pch) );
			}
			else if(strcmp(pch,"alt_id")==0){
				pch= strtok(NULL," GO:");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				altid_list.push_back( atoi(pch) );
			}
			else if(strcmp(pch,"replaced_by")==0){
				pch= strtok(NULL," GO:");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				replacedby_list.push_back( atoi(pch) );
			}
			else if(strcmp(pch,"consider")==0){
				pch= strtok(NULL," GO:");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				consider_list.push_back( atoi(pch) );
			}			
			else if(strcmp(pch,"relationship")==0){
				pch= strtok(NULL,": ");
				if(pch == NULL){
					cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
					exit(1);
				}
				if(strcmp(pch,"part_of")==0){
					pch= strtok(NULL," GO:");
					if(pch == NULL){
						cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
						exit(1);
					}
					partof_list.push_back( atoi(pch));
				}
				else if( strcmp(pch,"regulates")==0 || strcmp(pch,"positively_regulates")==0 || strcmp(pch,"negatively_regulates")==0 ){
					pch= strtok(NULL," GO:");
					if(pch == NULL){
						cerr << "ERROR: "<<line_num<<": "<<fin<<"\n";
						exit(1);
					}
					regulates_list.push_back( atoi(pch));
				}
			}
			else if(strcmp(pch,"is_obsolete")==0){
				is_obsolete= true;
			}
		}
	}
	
	return go_nodes;
}

namespace gograph{
const unsigned int CELLULAR_COMPONENT= 5575;
const unsigned int BIOLOGICAL_PROCESS= 8150;
const unsigned int MOLECULAR_FUNCTION= 3674;
}

struct go_graph{
	unordered_map <unsigned int,go_node> nodes;
	
	unordered_map <unsigned int,vector<unsigned int> > isa_parents;
	unordered_map <unsigned int,vector<unsigned int> > isa_children;
	
	unordered_map <unsigned int,vector<unsigned int> > partof_parents;
	unordered_map <unsigned int,vector<unsigned int> > partof_children;
	
	unordered_map <unsigned int,vector<unsigned int> > regulates_parents; // inlcudes "regulates","positively_regulates" and "negatively_regulates"
	unordered_map <unsigned int,vector<unsigned int> > regulates_children;
	
	unordered_map <unsigned int,vector<unsigned int> > altid_parents;
	unordered_map <unsigned int,vector<unsigned int> > altid_children;
	
	unordered_map <unsigned int,vector<unsigned int> > replacedby_parents;
	unordered_map <unsigned int,vector<unsigned int> > consider_parents;
};

/*
 * Creates go_graph from a list of go_nodes.
 *
 * Nodes are stored as id->go_node hash.
 * Edges are stored as id->parents_ids/children_ids hashes.
 *
 */
 
go_graph create_go_graph(const vector<go_node>& go_nodes){
	
	go_graph graph;
	
	for(unsigned int i=0; i<go_nodes.size(); i++){
	
	
		go_node node= go_nodes[i];
		
		// NODE HASH:
		graph.nodes[node.id]= node;
		
		
		
		// EDGE HASHES: child->parent edges
		// we use is_a|part_of|[positively_|negatively_]regulates|replaced_by|consider|alt_id
		
		if(node.isa_list.size()>0){
			if(graph.isa_parents.count(node.id) == 0){
				vector<unsigned int> temp;
				graph.isa_parents[node.id]= temp;
			}
			graph.isa_parents[node.id].insert(graph.isa_parents[node.id].end(),
							node.isa_list.begin(),
							node.isa_list.end());
		}
		if(node.partof_list.size()>0){
			if(graph.partof_parents.count(node.id) == 0){
				vector<unsigned int> temp;
				graph.partof_parents[node.id]= temp;
			}
			graph.partof_parents[node.id].insert(graph.partof_parents[node.id].end(),
							node.partof_list.begin(),
							node.partof_list.end());
		}
		if(node.regulates_list.size()>0){
			if(graph.regulates_parents.count(node.id) == 0){
				vector<unsigned int> temp;
				graph.regulates_parents[node.id]= temp;
			}
			graph.regulates_parents[node.id].insert(graph.regulates_parents[node.id].end(),
							node.regulates_list.begin(),
							node.regulates_list.end());
		}
		if(node.replacedby_list.size()>0){
			if(graph.replacedby_parents.count(node.id) == 0){
				vector<unsigned int> temp;
				graph.replacedby_parents[node.id]= temp;
			}
			graph.replacedby_parents[node.id].insert(graph.replacedby_parents[node.id].end(),
							node.replacedby_list.begin(),
							node.replacedby_list.end());
		}
		if(node.consider_list.size()>0){
			if(graph.consider_parents.count(node.id) == 0){
				vector<unsigned int> temp;
				graph.consider_parents[node.id]= temp;
			}
			graph.consider_parents[node.id].insert(graph.consider_parents[node.id].end(),
							node.consider_list.begin(),
							node.consider_list.end());
		}
		// note that altid_list contains children not parents
		// thus we invert the edge befor adding to the altid_parents edge-hash
		for(unsigned int j=0; j<node.altid_list.size(); j++){
			if(graph.altid_parents.count(node.altid_list[j]) == 0){
				vector<unsigned int> temp;
				graph.altid_parents[node.altid_list[j]] = temp;
			}
			graph.altid_parents[node.altid_list[j]].push_back(node.id);
		}
		
		
		// EDGE HASHES: parent->child edges
		// we use is_a|part_of|regulates|alt_id
		
		for(unsigned int j=0; j<node.isa_list.size(); j++){
			if(graph.isa_children.count(node.isa_list[j]) == 0){
				vector<unsigned int> temp;
				graph.isa_children[node.isa_list[j]] = temp;
			}
			graph.isa_children[node.isa_list[j]].push_back(node.id);
		}
		for(unsigned int j=0; j<node.partof_list.size(); j++){
			if(graph.partof_children.count(node.partof_list[j]) == 0){
				vector<unsigned int> temp;
				graph.partof_children[node.partof_list[j]] = temp;
			}
			graph.partof_children[node.partof_list[j]].push_back(node.id);
		}
		for(unsigned int j=0; j<node.regulates_list.size(); j++){
			if(graph.regulates_children.count(node.regulates_list[j]) == 0){
				vector<unsigned int> temp;
				graph.regulates_children[node.regulates_list[j]] = temp;
			}
			graph.regulates_children[node.regulates_list[j]].push_back(node.id);
		}
		
		// note that node.altid_list contains nodes children not parents
		if(node.altid_list.size() > 0){
			if(graph.altid_children.count(node.id) == 0){
				vector<unsigned int> temp;
				graph.altid_children[node.id]= temp;
			}
			graph.altid_children[node.id].insert(graph.altid_children[node.id].end(),
								node.altid_list.begin(),
								node.altid_list.end());
		}		
	}
	
	// TODO: are edge duplicates possible? -> remove
	
	return graph;
}



/*
 * Return a list of parent, child or relative ids for a given query id.
 *
 * Query id is not included in this list (that is node is not considered to be it's own parent/child)
 *
 * grah
 * query_id
 * relations	comma-sep list of edge types to be included in the search:"isa","partof","regulates","consider","replacedby"
 * direction	literal "parents","children" or "relatives"
 * max_num	max number of relatives to return
 * [verbal]	check those cases when no relatives are found and prints likely reasons to STDERR
 * [filter]	only return nodes that have keys in this filter
 *		An empty filter will not be used.
 *
 */

vector<unsigned int> get_relatives(go_graph& graph,
				unsigned int query_id,
				string relations,
				string direction, 
				unsigned int max_num,
				bool verbal= false,
				unordered_map<unsigned int,double> filter = unordered_map<unsigned int,double>()){
	
	
	// NOTE: pointers essential here, otherwise slowdown
	vector<unordered_map<unsigned int,vector<unsigned int> > *> edge_hash_list;
	
	vector<string> rel_list= split(relations,",");
	for(unsigned int i=0; i<rel_list.size(); i++){
		if(rel_list[i] == "isa"){
			if(direction == "parents")
				edge_hash_list.push_back(&graph.isa_parents);
			else if(direction == "children")
				edge_hash_list.push_back(&graph.isa_children);
			else if(direction == "relatives"){
				edge_hash_list.push_back(&graph.isa_parents);
				edge_hash_list.push_back(&graph.isa_children);
			}
		}
		else if(rel_list[i] == "partof"){
			if(direction == "parents")
				edge_hash_list.push_back(&graph.partof_parents);
			else if(direction == "children")
				edge_hash_list.push_back(&graph.partof_children);
			else if(direction == "relatives"){
				edge_hash_list.push_back(&graph.partof_parents);
				edge_hash_list.push_back(&graph.partof_children);
			}		
		}
		else if(rel_list[i] == "regulates"){
			if(direction == "parents")
				edge_hash_list.push_back(&graph.regulates_parents);
			else if(direction == "children")
				edge_hash_list.push_back(&graph.regulates_children);
			else if(direction == "relatives"){
				edge_hash_list.push_back(&graph.regulates_parents);
				edge_hash_list.push_back(&graph.regulates_children);
			}		
		}
		else if(rel_list[i] == "altid"){
			if(direction == "parents")
				edge_hash_list.push_back(&graph.altid_parents);
			else if(direction == "children")
				edge_hash_list.push_back(&graph.altid_children);
			else if(direction == "relatives"){
				edge_hash_list.push_back(&graph.altid_parents);
				edge_hash_list.push_back(&graph.altid_children);
			}
		}
		else if(rel_list[i] == "replacedby"){
			if(direction == "parents")
				edge_hash_list.push_back(&graph.replacedby_parents);
		}
		else if(rel_list[i] == "consider"){
			if(direction == "parents")
				edge_hash_list.push_back(&graph.consider_parents);
		}
		else{
			cerr << "ERROR: <gograph.h>: get_relatives(): invalid par: \'"<<rel_list[i]<<"\'\n";
			exit(1);
		}
	}
	
	if( !(direction == "parents" || direction == "children" || direction == "relatives") ){
		cerr << "ERROR: <gograph.h>: get_relatives(): invalid par: \'"<<direction<<"\'\n";
		exit(1);
	}
	
	
	bool use_filter = false;
	if( !filter.empty() )
		use_filter = true;
	
	queue<unsigned int> q;
	unsigned int id_temp;
	vector<unsigned int> list_of_relatives;
	unordered_map<unsigned int,unsigned int> nodes_visited;
	
	q.push(query_id);
	nodes_visited[query_id]=1;
	
	
	while( !q.empty() && list_of_relatives.size()<max_num){
		id_temp= q.front();
		q.pop();
		
		for(unsigned int i=0; i<edge_hash_list.size(); i++){
			
			if(edge_hash_list[i]->count(id_temp)>0){
				vector<unsigned int> relatives_temp= (*edge_hash_list[i])[id_temp];
				
				for(unsigned int j=0; j<relatives_temp.size(); j++){
					if(nodes_visited.count(relatives_temp[j])==0){
						q.push(relatives_temp[j]);
						nodes_visited[relatives_temp[j]]=1;
					
						if(!use_filter  ||  (use_filter && filter.count(relatives_temp[j])>0) )
							list_of_relatives.push_back(relatives_temp[j]);
					}
				}
			}
		}
	}				
			
	
	// truncate relatives to max_num
	if(list_of_relatives.size() > max_num)
		list_of_relatives.resize(max_num);


	// VERBAL MODE CHECK
	if(verbal && list_of_relatives.size()==0){
		
		if(direction == "children"){ // no children this is normal: do nothing
		}
		else if(direction == "parents"){
			if(graph.nodes.count(query_id) == 0){
				fprintf(stderr,"WARNING: no parents for GO:%07u: node is off the graph\n",query_id);
			}
			else if(graph.nodes[query_id].is_obsolete){
				//fprintf(stderr,"WARNING: no parents for GO:%07u: obsolete node\n",query_id);
				// OK: common
			}
			else if(query_id==gograph::CELLULAR_COMPONENT 
				|| query_id==gograph::BIOLOGICAL_PROCESS 
				|| query_id==gograph::MOLECULAR_FUNCTION){ // OK
			}
			else{
				fprintf(stderr,"WARNING: no parents for GO:%07u: unknown reason\n",query_id);
			}
		}
		else if(direction == "relatives"){
			if(graph.nodes.count(query_id) == 0){
				fprintf(stderr,"WARNING: no relatives for GO:%07u: node is off the graph\n",query_id);
			}
			else if(graph.nodes[query_id].is_obsolete){
				//fprintf(stderr,"WARNING: no relatives for GO:%07u: obsolete node\n",query_id);
				// OK
			}
			else{
				fprintf(stderr,"WARNING: no relatives for GO:%07u: unknown reason\n",query_id);
			}
		}

	}
		
	return list_of_relatives;
}


/*
 * Returns hash map containing lists of parents for each node id in a given graph.
 *
 * graph			go_graph
 * relations			comma-sep list of relations, e.x. "isa,partof,consider,replacedby,altid"
 * include_parents_for_altids	for each alt id will include parets for node, that alt-id is pointing to
 * include_node_itself		ensures that all nodes have at least one parent
 * verbal
 * 
 */ 
unordered_map<unsigned int,vector<unsigned int> > get_parent_map(
					go_graph& graph, 
					string relations, 
					bool include_parents_for_altids	= false, 
					bool include_node_itself	= false, 
					bool verbal			= false){
	
	unordered_map<unsigned int,vector<unsigned int> > pmap;
	
	for (auto it = graph.nodes.begin(); it != graph.nodes.end(); ++it ){
		
		
		vector<unsigned int> parents= get_relatives(graph,
								(it->second).id,
								relations,
								"parents",
								10000,
								verbal);
		if(include_node_itself)
			parents.push_back((it->second).id);
			
		sort(parents.begin(),parents.end());
		pmap[(it->second).id]= parents;
		
		if(include_parents_for_altids){
			for(unsigned int j=0; j< (it->second).altid_list.size(); j++){
				pmap[(it->second).altid_list[j]] = parents;
			}
		}
	}
	
	return pmap;
}
