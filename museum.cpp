#include <iostream>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <limits>
#include "museum.h"
#include <queue>
#include <cstring>
#include <iomanip>
#include <list>
#include <iterator>
#include <algorithm>

using namespace std;

void printHelp();

//the premise of the project is:
/*
 Overview
 Welcome to the Mammoth Museum! This unique museum features mammoths frozen in ice,
 preserved perfectly to this very day! As a curator, you want to plan out the most efficient way to
 peruse every exhibit for your interested visitors.
 Your first task in part A is to plan the pathways connecting the different exhibits efficiently. (That
 is, minimize the total length of the pathways needed.)
 Then in parts B & C, you will plan an efficient daily route to check on each of the exhibits. (That
 is, minimize the distance you must travel to visit all the exhibits.)
 To be clear, these scenarios are separate, your program will create a plan for one or the
 other, but not both in the same run (although you may find that the algorithms from one mode
 help with another mode).
*/

//Generates different paths
void genPerms(deque<Exhibit2> &q, vector <Exhibit2> &s, double &upper_bound, Weight &w, Exhibit2 zero, vector<Exhibit2> &path, vector < vector <double> > &dist_matrix);
//finds endes to connect the path to the MST tree
double lowest_edges(deque<Exhibit2> &mydeque, Exhibit2 &front, Exhibit2 &back, vector < vector <double> > &dist_matrix);
//calculates MST tree
void MST_calc(deque<Exhibit2> mst_d, double &lower_bound, vector < vector <double> > &dist_matrix);

int main(int argc, char * argv[]) {
    
	cout << setprecision(2);
	cout << fixed;
    //takes opt args from the command line
	struct option longOpts[] = {
		{"help", 		no_argument, 		NULL, 'h'},
		{"mode", 		required_argument, 	NULL, 'm'},
		{NULL,			0,					NULL,  0 }
	};
    
    
	opterr = 1;
    bool MST = false;
    bool fast = false;
	int opt = 0, index = 0;
	bool optimal = false;
	while ((opt = getopt_long(argc,argv,"hm:",longOpts,&index)) != -1) {
		switch (opt) {
			case 'h':
				printHelp();
				return 0;
			case 'm':
                if (strcmp(optarg, "MST") == 0) {
                    MST = true;
                }
                if(strcmp(optarg, "FASTTSP") == 0) {
                	fast = true;
                }
                if(strcmp(optarg, "OPTTSP") == 0) {
                	optimal = true;
                }
				break;
			case '?':
			default:
				printHelp();
				exit(1);
		}
	}
    
	vector <Coordinate> inside; //inside the museum
	vector <Coordinate> outside; //outside the museum
    
	priority_queue <Pair, vector<Pair>, PairComp> pairs; //used to print out the output
    
    //the input comes in like this
    //first the boundry is given
    //then the number of coordinates
    //then the 2d coordinates
    /*­
     
     3 0 3 4
     7
     ­4 3
     0 5
     ­1 4
     7 ­2
     3 8
     2 1
     0 2
     
     */
    
    //this MST tree is seperated into outside of the entrance and inside.
    
	Weight end_weight;
	int lx, ly, rx, ry; //boundry variables
	cin >> lx >> ly >> rx >> ry;
	int vertices; // number of vertices
	cin >> vertices;
	int x, y; //
    
    //IF MST OPTION IS GIVEN
    
	if(MST){
        //MST uses a version of Prim's algorithm because for this particular algorithm
        // and tsp we are given a connected weighted undirected graph
        /* prim's works as follows (taken from the wikipedia page on prim's algorithm
         
         1. Initialize a tree with a single vertex, chosen arbitrarily from the graph.
         2. Grow the tree by one edge: of the edges that connect the tree to vertices not yet in the tree, find the minimum-weight edge, and transfer it to the tree.
         3. Repeat step 2 (until all vertices are in the tree).
         */
		int j = 0;
		bool entrance = false;
		while(j < vertices){
			cin >> x >> y;
			if( ((x > lx) and (x < rx)) and ((y > ly) and (y < ry))) {
				inside.push_back(Coordinate(j, x, y));
			}
            //finding the entrance by checking if a coordinate is on the boundry
			else if( (((x == lx) or (x == rx)) and ((y >= ly) and (y <= ry)))
                    or (((y == ly) or (y == ry)) and ((x >= lx) and (x <= rx))) ){
				Coordinate startcoord = Coordinate(j, x, y);
				inside.push_back(startcoord);
				outside.push_back(startcoord);
	            entrance = true;
			}
			else {
				outside.push_back(Coordinate(j, x, y));
			}
	        j++;
		}
		if((entrance != true) and (!inside.empty() and !outside.empty())){
            //if there is no entrance and there are coordiniates inside and outside of the museum
            //then you can't construct a tree.
			cerr << "Cannot construct MST!" << endl;
			exit(1);
		}
        //the algorithm will construct an actual tree only if you have more then one node
        //the entrance appears in both the inside and outside so if there is only one,
        // this is the entrance.
	    if(inside.size() > 1){
            //sets the least distance to the first coordinates distance. it starts out as 0.
	        double least_distance = inside[0].get_distance();
            //then it gets visited.
	        inside[0].visited_now();
	        long least_index = 0;
	        unsigned long i = 0;
	        long newleast = least_index;
            //loops through the vector
	        while(i < inside.size()){
	            inside[least_index].visited_now();
	            for(unsigned long j = 0; j < inside.size(); j++){
	                if(!inside[j].is_visited()){ //checks if the vertex was not visited
	                    inside[j].connection(inside[least_index]);
                        //what the connection function does is:
                        //calculates the distance between the coordinate with index j
                        //and the coordinate with least_index and then saves this
                        //edge with least index being prevcoord (refer to coordinate class).
	                }
	                if((inside[j].get_distance() <= least_distance) and !inside[j].is_visited()) {
                        //if this edge is less then or equal to the least_distance and the vertex is not
                        //visited make the least disance the edge between the index j and the least_index
	                    least_distance = inside[j].get_distance();
	                    newleast = j;
                        //newleast is j
                        //newleast is to save the index but not change it while going through the for loop
                        //I want to check least_index with every j first.
	                }
	            }
	            if(least_index != 0){
                    //skip the first one
	                pairs.push(Pair(inside[least_index].prevcoord, inside[least_index].get_index()));
	                end_weight.add_weight(inside[least_index].get_distance());
                    //save the pair between the least_index's prevcoord and the least_index coordinate.
	            }
	            least_distance = numeric_limits<double>::infinity(); //must make least distance infinity
                //because we are restarting the process.
	            least_index = newleast;
                //now change least to newleast
	            i++;
	        }
	    }
        //deos the same thing as inside. I know, I know repeated code = bad.
		if(outside.size() > 1){
	        double least_distance = outside[0].get_distance();
	        outside[0].visited_now();
	        long least_index = 0;
	        unsigned long i = 0;
	        long newleast = least_index;
	        while(i < outside.size()) {
	            outside[least_index].visited_now();
	            for(unsigned long j = 0; j < outside.size(); j++){
	                if(!outside[j].is_visited()){
	                    outside[j].connection(outside[least_index]);
	                }
	                if((outside[j].get_distance() <= least_distance) and !outside[j].is_visited()) {
	                    least_distance = outside[j].get_distance();
	                    newleast = j;
	                }
	            }
	            if(least_index != 0){
	                pairs.push(Pair(outside[least_index].prevcoord, outside[least_index].get_index()));
	                end_weight.add_weight(outside[least_index].get_distance());
	            }
	            least_distance = numeric_limits<double>::infinity();
	            least_index = newleast;
	            i++;
	        }
		}
		end_weight.print_weight();
        //print the pairs
		while(!pairs.empty()){
			Pair thispair = pairs.top();
			thispair.print_pair();
	        pairs.pop();
		}
	}
	double upper_bound = numeric_limits<double>::infinity();
    //FAST TSP OPTION
	if(fast) {
        //what fast tsp does is a form of the greedy algorithm but of only being able
        //to insert an edge at the end of the path you can insert an edge anywhere.
		list<Exhibit> path;
		cin >> x >> y;
		Exhibit exhibit;
		int i = 0;
		exhibit.set_params(x, y, i);
		path.push_front(exhibit);
		path.push_back(exhibit);
		i++;
 		double min;
		double weight = 0;
        
        //what is happening is that the list starts like this:
        //0<->0 then as you insert edges
        //0<->1<->0 the fiunction calculates the weight of the path
        /*    As you add edges it finds the best place to put the the vertex so that the 
                edge weight is minimal
         
                0-------1
                 \     /
                  \   /   3
                   \ /
                    2

            where can you place that 2 so that the total weight of the path is minimal
            turns out it is pretty close to the optimal solution but much faster yayy!!
         
            the best place is 0<->1<->3<->2<->0
         */
		while(i < vertices){
 			min = numeric_limits<double>::infinity();
 			auto least = next(path.begin()); //using the iterator class here
 			cin >> x >> y;
 			exhibit.set_params(x, y, i);
	 		for(auto it1 = path.begin(); it1 != prev(path.end()); it1++){
	 			auto it2 = next(it1);
	 			double temp = calculate_distance(exhibit.x, exhibit.y, it1->x, it1->y) +
	 			calculate_distance(exhibit.x, exhibit.y, it2->x, it2->y);
	 			double inbetween = calculate_distance(it1->x, it1->y, it2->x, it2->y);
	 			if((temp - inbetween) < min) {
	 				min = temp - inbetween;
	 				least = it2;
	 			}
	 		}
	 		weight += min;
	 		path.insert(least, exhibit);
	 		i++;
		}

        if(fast){
			path.pop_back();
			cout << weight << endl;
			for (auto it = path.begin(); it != path.end(); it++){
				cout << it->index << " ";
			}
		}
	}
	if(optimal){
		Exhibit2 e;
		Exhibit2 zero;
		deque<Exhibit2> d;
		vector<Exhibit2> empt;
		int i = 0;
		cin >> x >> y;
		zero.set_params(x, y, i);
		list<Exhibit2> path;
		Exhibit2 exhibit;
        exhibit.set_params(x, y, i);
        i++;
		path.push_front(exhibit);
		path.push_back(exhibit);
        double min;
		double weight = 0;
		vector < vector <double> > dist_matrix(vertices, vector<double> (vertices));
        
        //to get an excellent lower bound I use the same fast tsp method I used above
		while(i < vertices){
 			min = numeric_limits<double>::infinity();
 			auto least = next(path.begin());
 			cin >> x >> y;
 			exhibit.set_params(x, y, i);
 			e.set_params(x, y, i);
			d.push_back(e);
	 		for(auto it1 = path.begin(); it1 != prev(path.end()); it1++){
	 			auto it2 = next(it1);
	 			double temp = calculate_distance(exhibit.x, exhibit.y, it1->x, it1->y) +
	 			calculate_distance(exhibit.x, exhibit.y, it2->x, it2->y);
	 			double inbetween = calculate_distance(it1->x, it1->y, it2->x, it2->y);
	 			if((temp - inbetween) < min) {
	 				min = temp - inbetween;
	 				least = it2;
	 			}
	 		}
	 		weight += min;
	 		path.insert(least, exhibit);
	 		i++;
        }
		path.pop_back();
		vector<Exhibit2> v(path.begin(), path.end());
		upper_bound = weight;
		d.push_front(zero);
        //ready up the distance matrix for precalculated distances!!
		for(int a = 0; a < vertices; a++){
			for(int b = 0 ; b < vertices; b++){
				dist_matrix[a][b] = calculate_distance(d[a].x, d[a].y, d[b].x, d[b].y);
			}
		}

		Weight total;
		genPerms(d, empt, upper_bound, total, zero, v, dist_matrix);
		total.weight = upper_bound;
		total.print_weight();
        
		for(unsigned int j = 0; j < v.size(); j++){
			cout << v[j].index << " ";
		}
	}
	return 0;
}//end of MAIN


// Generate Permutations: Recursively generates permutations using a deque in place of a queue and a vector in place of a stack to tailor the calculation of permutations to what I need to use them for. generate a path in a vector and the deque is necessary to access both ends to get rid of cyclical permutations.
void genPerms(deque<Exhibit2> &q, vector <Exhibit2> &s, double &upper_bound, Weight &w, Exhibit2 zero, vector<Exhibit2> &path, vector < vector <double> > &dist_matrix){
	unsigned long k, size = q.size();

	if(q.empty()){
        //if the queue is empty that means the permutation is complete so check if the current path's weight is less than the upper bound
		if(w.weight + dist_matrix[0][s.back().index] < upper_bound){
            //if it is then make the upperbound the new path weight
			upper_bound = w.weight + dist_matrix[0][s.back().index];
			path = s;
		}
        
		return;
	}
	else{
        
		for(k = 0; k!= size; k++){
			if(s.empty()){
				w.weight = 0;
                
				s.push_back(zero);
                
				q.pop_front();
			}
			w.add_weight(dist_matrix[q.front().index][s.back().index]);
			s.push_back(q.front());
			q.pop_front();
			double lower_bound = 0;
			MST_calc(q, lower_bound, dist_matrix);
            //mst is as it turns out a good prediction for the whether or not the current state
            // of the permutation is good otherwise known as the lower bound.
			if((s[1].index != q.back().index)/* this condition is to prevent cyclical permutations */ and (w.weight + lower_bound + lowest_edges(q, zero, s.back(), dist_matrix) < upper_bound)){
                //if the current weight and the lower_bound (MST Tree) plus the weight of the
                //connected portion from the path to the tree is less then the upperbound
                //continue with this path as in call genperms again.
				genPerms(q, s, upper_bound, w, zero, path, dist_matrix);
			}
			else {
				//cerr << "pruning!" << endl;
                // ^ debugging left over that I kept for understanding what is happening.
                //if the current place in the permutation doesn't pass the above condition
                //it will be cut short because genPerms wasn't called.
			}
            //how this works is that the algorithm will keep making different permutations of
            //the list of inidces. what is currently in the queue are the unused indices and
            // the stack is the current permutation. if genperms isn't called to finish the
            //permutation the queue will become empty and the next permutation is generated.
			q.push_back(s.back());s.pop_back();
			w.weight -= dist_matrix[q.back().index][s.back().index];
		}
	}
}//end of genPerms

//this function just connects the last vertex currently in the deque with the closest vertex in the mst tree

double lowest_edges(deque<Exhibit2> &mydeque, Exhibit2 &front, Exhibit2 &back, vector < vector <double> > &dist_matrix) {
	double leastdistance1 = numeric_limits<double>::infinity();
	double leastdistance2 = numeric_limits<double>::infinity();
	if(mydeque.empty())
		return 0;
	for(unsigned i = 0; i < mydeque.size(); i++){
		double distance1 = dist_matrix[mydeque[i].index][front.index];
		double distance2 = dist_matrix[mydeque[i].index][back.index];
		if(distance1 < leastdistance1){
			leastdistance1 = distance1;
		}
		if(distance2 < leastdistance2) {
			leastdistance2 = distance2;
		}
	}
	return leastdistance1 + leastdistance2;
}

//another version of the MST tree I used before but it uses the distance matrix and 2 for loops in a for loop.
void MST_calc(deque<Exhibit2> mst_d, double &lower_bound, vector < vector <double> > &dist_matrix){

	if(mst_d.empty()){
		return;
	}
	mst_d[0].distance = 0;
	mst_d[0].prevcoord = 0;
    
	for(unsigned i = 0; i < mst_d.size(); i++){
		unsigned least = 0;
		double least_distance = numeric_limits<double>::infinity();
		for(unsigned j = 0; j < mst_d.size(); j++){
			if(!mst_d[j].visited){
				if(mst_d[j].distance < least_distance){
					least_distance = mst_d[j].distance;
					least = j;
				}
			}
		}
		mst_d[least].visited = true;
		lower_bound += mst_d[least].distance;
		for (unsigned j = 0; j < mst_d.size(); ++j)
		{
            if(!mst_d[j].visited){
                least_distance = dist_matrix[mst_d[j].index][mst_d[least].index];
                if(least_distance < mst_d[j].distance){
                    mst_d[j].distance = least_distance;
                    mst_d[j].prevcoord = mst_d[least].index;

                }
            }
		}
	}
}

//when help is passed in as an arguement.
static const char helpText[] = "You need help.\nBut I won't help you!\nMuaaahahaha!\n";
void printHelp() {
	cout << helpText << flush;
}
//calculate the distance between two vertices.
double calculate_distance(long ax, long ay, long bx, long by){
    double xsq = (ax - bx) * (ax - bx);
    double ysq = (ay - by) * (ay - by);
    double param = xsq + ysq;
    double dist = sqrt(param);
    return dist;
}