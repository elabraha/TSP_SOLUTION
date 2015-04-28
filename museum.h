#ifndef MUSEUM_H_
#define MUSEUM_H_

#include <cmath>
#include <limits>
using namespace std;

//this class is for the MST tree used for the MST option
//I didnt really need a class especially since I have a distance function
//the only really benefit is that everything pertaining to creating the
//mst tree is contained in one class.
class Coordinate {
	bool visited;
	double distance;
	
public:
    long x;
    long y;
    long index;
    long prevcoord;
    
    Coordinate(long index_in, long x_in, long y_in){
        x = x_in;
        y = y_in;
        index = index_in;
        visited = false;
        distance = numeric_limits<double>::infinity();
        prevcoord = -1;
    };
    ~Coordinate(){};
    
    //this function just caculates the distance between two coordinates.
    //then makes sure the vertex is connected to the previous one by setting
    //prevcoord to it's index
    void connection (Coordinate coord){
        double xsq = (x - coord.x) * (x - coord.x);
        double ysq = (y - coord.y) * (y - coord.y);
        double param = xsq + ysq;
        double newdist = sqrt(param);// distance = sqrt((x-dx)^2+(y-dy)^2))
        if(newdist < distance){
            distance = newdist;
            prevcoord = coord.index;
        }
    };
    
    void change_distance(double newdist) {
        distance = newdist;
    };
    
    double get_distance() {
        return distance;
    };
    
    bool is_visited(){
        return visited;
    };
    
    void visited_now(){
        visited = true;
    };
    
    long get_index() const {
        return index;
    };
    
};
//as part of this project we needed to consider
//space optimizations as well so i got rid of
//distance, visited, and prevcoord.
struct Exhibit {
	long index;
	long x;
	long y;
	void set_params(long x_in, long y_in, long index_in){
		index = index_in;
		x = x_in;
		y = y_in;
	};
};

//This is a modification used for mst tree using a struct because
//Coordinate was fairly clunky for what I needed.
struct Exhibit2 {
	long index;
	long x;
	long y;
	bool visited = false;
	long prevcoord = 0;
	double distance = numeric_limits<double>::infinity();
	void set_params(long x_in, long y_in, long index_in){
		index = index_in;
		x = x_in;
		y = y_in;
	};
};
//it wasn't necessary to create a data structure I did it because
//I actually thought I needed more information in my fuctions besides
//weight but I didn't and I ended up not getting rid of it.
struct Weight {
	double weight = 0;
	void add_weight(double more) {
		weight += more;
	};
    
	void print_weight(){
		cout << weight << endl;
	};
};

//this is just for printing pairs of coordinates
//state reason
struct Pair {
	long first;
	long second;
	void print_pair(){
		cout << first << " " << second << "\n";
	};
	Pair(long a, long b){
		if(a < b) {
			first = a;
			second = b;
		}
		if(a > b) {
			first = b;
			second = a;
		}
	};
};
//for ordering pairs
struct PairComp{
    bool operator() (Pair a, Pair b){
    	if(a.first == b.first){
    		return a.second > b.second;
    	}
        return a.first > b.first;
    };
    
};
//declare calculate distance
double calculate_distance(long ax, long ay, long bx, long by);

#endif