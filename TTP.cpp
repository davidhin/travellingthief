#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <numeric>
#include <limits.h>
#include <unordered_map>
#include <math.h> //sqrt & pow

#include "prettyprint.h" // custom .h file for printing
using namespace std;

// Encapsulates a single City
struct City {
    int index;
    int xCoord;
    int yCoord;
    bool added; // If City has been added to
};

// Encapsulates a single Item
struct Item {
    int index;
    double profit;
    double weight;
    int city;
    bool picked = false;
};

// Print an item by printing its index
ostream &operator<<(ostream &os, Item* const &i) {
    return os << i->index;
    //return os << "<id:" << i->index <<",city:" << i->city << ",weight:" << i->weight << "," <<"profit:" << i->profit << ">";
}

// Print an item by printing its index
ostream &operator<<(ostream &os, Item const &i) {
    //return os << i.index;
    return os << "<id:" << i.index <<",city:" << i.city << ",weight:" << i.weight << "," <<"profit:" << i.profit << ">";
}

// Forward Declarations
double Z(vector<int>*, vector<double>*, vector<Item*>*); // Objective function
void insertion(vector<int>*); // Objective function
void bitflip(vector<int>*, double*);
double distance(City*, City*);
vector< pair<double, int> > get_scores(double); // Returns items in "score" order
void parsingFile(char* argv[]);
vector<Item*> packing_plan();

// Globals
int dimension;
double minS;
double maxS;
int ksc;
double rent;
int itemNum;

// Constants
int tau = 2000 ;   // We may not need this because our input sizes are small
                // enough to just try every item in the packing plan
                // However, we could use this when we start optimising tours (?)

// Store order of the cities in the tour and the distance to the next City
vector<int> tour;
vector<double> dist;
vector<City> cities; // Parsed cities
vector<Item> items; // Parsed Items

int main(int argc, char* argv[]) {
    // Parsing the files
    parsingFile(argv);

    // **** Tour Optimization **** //
    // Storing old vectors
    vector<Item*> P_ultimate;
    vector<int> old_tour = tour;
    vector<double> old_dist = dist;
    double Z_ultimate = -0xfffffff;

    // Optimization parameters
    int start = 1; // Optimization parameter
    int numIter = 0; // Number of iterations

    // Main loop
    //for(int tourUpdate=start; tourUpdate<dimension; tourUpdate+=10) {
    for(int tourUpdate=start; tourUpdate<dimension-1; tourUpdate+=1) {
        for(int tour2=tourUpdate+1; tour2<3; tour2+=1) {
        // This is the index at which the 2 opt will split the tour
        int firstEdge = tourUpdate;
        int secondEdge = tour2;

        if(numIter!=0) { // Calculate objective value for original tour first
            // Tour optimization goes starts here
            // Reverse middle segment
            reverse(tour.begin()+firstEdge, tour.begin()+secondEdge); // The whole of 2-opt essentially
            // Tour optimization goes ends here

            // Recalculating distances
            for(int i=firstEdge-1; i<secondEdge; i++) {
                dist[i+1] = distance(&cities[tour[i]], &cities[tour[i-1]]);
            }
        }
        // Calculate distance from last City to first
        dist[dist.size()-1] = distance(&cities[tour[0]-1], &cities[tour[tour.size()-1]-1]);

        // Packing Plan Algorithm Starts Here
        vector< pair<double, int> > score_items = get_scores(1); // Get a list of items ordered by score

        // vector<Item*> P_best ; // Initialise current and best packing plan
        double W_curr = 0, W_cap = ksc ; // Initialise current weight and weight capacity
        int mu = floor(itemNum / tau) ; // Set the frequency of item adding
        double Z_best = -0xfffffff ; // initialise best objective value to -inf
        int k = 0; //, k_last=0; // used to iterate through score_items
        while (W_curr < W_cap && mu > 1 && k < itemNum) {
            Item* kth_item = &items[score_items[k].second-1];
            if (W_curr + kth_item->weight <= W_cap) {
                kth_item->picked = true;
                W_curr += kth_item->weight;
                if (k % 1 == 0) { // Set to 1 if want to check all items
                    vector<Item*> plan = packing_plan();
                    double Z_curr = Z(&tour, &dist, &plan) ;
                    if (Z_curr < Z_best) {
                        kth_item->picked = false;
                        W_curr -= kth_item->weight ;
                        //k = k_last ;
                        //mu /= 2 ; // Truncate / floor
                    }
                    else {
                        // k_last = k ;
                        Z_best = Z_curr ;
                    }
                    //cout << Z_best << "\n";
                    //cout << W_curr << "\n";
                }
            }
            k++;
         cout<<Z_best<<endl;
         cout << k << "\n";
        }

        bitflip(&old_tour, &W_curr);
        //insertion(&old_tour) ;
        // cout<<Z(&old_tour, &dist)<<" after flip"<<endl;

        // Reverting to old tour and distances
        if(Z_best<Z_ultimate) { // worse
            tour = old_tour;
            dist = old_dist;
        }else{// better
            Z_ultimate = Z_best;
            old_tour = tour;
            old_dist = dist;
            // P_ultimate = P_best;
            //cout << "Improved" << "\n";
        }
        numIter ++;
        }
    }

    // Output to file
    ofstream outputFile;
    outputFile.open("fnl_soln.ttp");
    outputFile << old_tour << "\n" ;
    outputFile << packing_plan() << "\n" ;
    outputFile.close();
}

vector<Item*> packing_plan() {
    vector<Item*> plan;
    for(int i=0; i<(int)items.size(); i++)
        if (items[i].picked)
            plan.push_back(&items[i]);
    return plan;
}

// Distance between two cities
double distance(City* a, City* b) {
    double temp1 = 0, temp2 = 0;
    temp1 = a->xCoord - b->xCoord;
    temp2 = a->yCoord - b->yCoord;
    temp1 = pow(temp1, 2);
    temp2 = pow(temp2, 2);
    temp1 = temp1 + temp2;
    return sqrt(temp1);
}

// Objective function
double Z(vector<int>* tour_in, vector<double>* dist, vector<Item*>* plan) {
    long collected = 0;
    double ret = 0;
    double v = ( maxS - minS ) / ksc;

    for (int i = 0; i < (int)(*tour_in).size(); i++) {
        for (int j = 0; j<(int)plan->size(); j++) {
            if ((*plan)[j]->city == (*tour_in)[i]) {
                collected += (*plan)[j]->weight;
                ret += (*plan)[j]->profit;
            }
        }
        ret -= ceil((*dist)[i+1]) * rent / (maxS - v * collected);
    }
    return ret;
}

// Item score calculation
vector< pair<double, int> > get_scores(double alpha) {
    // Create a vector of "distances to end of tour"
    vector<double> end_dists;
    for (int i=0; i<(int)dist.size(); i++) {
        end_dists.push_back(accumulate(dist.begin()+i+1, dist.end(), 0.0));
    }

    // Assign items to cities
    unordered_map<int, vector<Item*>> city_items;
    for (int i=0; i<(int)items.size(); i++) {
        city_items[items[i].city].push_back(&items[i]);
    }

    // Print city/item map
    for (auto i = city_items.begin(); i != city_items.end(); i++) {
        vector<int> temp;
        for (int j=0; j<(int)i->second.size(); j++) {
            temp.push_back(i->second[j]->index);
        }
    }

    // Create the vector of items to read in order of score
    vector< pair<double, int> > score_items;

    for (int i=0; i<(int)tour.size(); i++) {
        // Get list of items at city tour[i]
        unordered_map<int, vector<Item*>>::iterator it = city_items.find(tour[i]+1);
        if (it == city_items.end()) continue;

        // Pointer to list of items at the city tour[i]
        vector<Item*> *current_items = &it->second;

        // Add item to "full" item list
        for (int j=0; j<(int)current_items->size(); j++) {
            Item* current_item = current_items->at(j);
            double score = pow(current_item->profit,alpha) / (pow(current_item->weight,alpha) * end_dists[j]);
            score_items.push_back(make_pair(score,current_item->index));
        }
    }

    // ** b. Sort the items of M in non-decreasing order of their scores ** //
    sort(score_items.begin(), score_items.end());

    return score_items;
}

// Bitflip
void bitflip(vector<int>* tour_in, double* W_curr) {
    vector<Item*> plan_init = packing_plan();
    double Z_score = Z(tour_in, &dist, &plan_init);
    double Z_flipped;

    for(int i=0; i<(int)items.size(); i++) {
        if (items[i].picked) { // If already picked
            items[i].picked = false;
        }else{ // If not in tour
            items[i].picked = true;
            if(*W_curr+items[i].weight>ksc) {
                items[i].picked = false;
                continue;
            }
        }

        vector<Item*> plan = packing_plan();
        Z_flipped = Z(tour_in, &dist, &plan);
        cout<< i << " " << Z_flipped<<" flipped"<<endl;
        // cout<<packing_plan()<<endl;
        if(Z_flipped<Z_score) { // Worse
            if (items[i].picked) {
                items[i].picked = false;
            }else {
                items[i].picked = true;
            }
        }else { // Better when flipped
            Z_score = Z_flipped;
            if (items[i].picked) {
                *W_curr += items[i].weight;
            }else {
                *W_curr -= items[i].weight;
            }

        }

    }
}

// Insertion
void insertion(vector<int>* tour_in) {
    vector<Item*> plan = packing_plan();

    for (int i = tour_in->size()-1; i>=0; i--) {
        double Z_best = Z(tour_in, &dist, &plan);
        cout << "Z_best: " << Z_best << "\n";

        vector<int> tour_starstar = *tour_in;
        vector<int> tour_new = *tour_in;

        for (int j=i; j>=1; j-=1000) {
            // Insert city between xb and xb-1 in tour_new
            swap(tour_new[j], tour_new[j-1]);
            // Compute Z*
            for (int k=max(j-1,1); k<min(j+2, (int)dist.size()-1); k++) {
                dist[k] = distance(&cities[tour_new[k-1]-1], &cities[tour_new[k]-1]);
            }
            dist[dist.size()-1] = distance(&cities[tour_new[tour_new.size()-1]-1], &cities[tour_new[0]-1]);
            double Z_curr = Z(&tour_new, &dist, &plan);
            //cout << Z_curr << "\n";
            // Update if better
            if (Z_curr > Z_best) {
                Z_best = Z_curr;
                tour_starstar = tour_new;
                cout << "updated" << "\n";
            }    
        }

        // Update actual tour to tour_starstar
        *tour_in = tour_starstar;
        for (int k=1; k<(int)dist.size()-1; k++) {
            dist[k] = distance(&cities[(*tour_in)[k-1]-1], &cities[(*tour_in)[k]-1]);
        }
        dist[dist.size()-1] = distance(&cities[(*tour_in)[(*tour_in).size()-1]-1], &cities[(*tour_in)[0]-1]);
    }
    cout << "\n";
}

// Parsing input file
void parsingFile(char* argv[]) {
    ifstream file(argv[1]);
    string line; // Holds each line of the file

    for(int i=0; i<3; i++) {
        getline(file, line);
    }

    // Dimension (number of cities)
    string temp;
    istringstream ss(line);
    ss >> temp >> dimension;
    ss.clear();

    // Number of items
    getline(file, line);
    ss.str(line);
    ss >> temp >> temp >> temp >> itemNum;
    ss.clear();

    // Getting knapsack capaCity
    getline(file, line);
    ss.str(line);
    ss >> temp >> temp >> temp >> ksc;
    ss.clear();

    // Min speed
    getline(file, line);
    ss.str(line);
    ss >> temp >> temp >> minS;
    ss.clear();

    // Max speed
    getline(file, line);
    ss.str(line);
    ss >> temp >> temp >> maxS;
    ss.clear();

    // Rent
    getline(file, line);
    ss.str(line);
    ss >> temp >> temp >> rent;
    ss.clear();

    // Skipping extra lines
    getline(file, line);
    getline(file, line);

    // Parsing City coordinates into vector of City structs
    for(int i=0; i<dimension; i++) {
        City newC;
        getline(file, line);
        ss.str(line);
        ss >> newC.index;
        ss >> newC.xCoord;
        ss >> newC.yCoord;
        newC.added = false;
        cities.push_back(newC);
        ss.clear();
    }

    // **** 2. Nearest Neighbours Algorithm **** //
    // Adding first City to the tour
    tour.push_back(0);
    dist.push_back(0);
    cities[0].added = true;

    // For storing current closest City details
    double lowDist;
    int tempInd = 1;
    double curr_dist;

    // Getting tour
    for(int l=1; l<dimension; l++) {
        lowDist = INT_MAX; // reset distance
        for(int k=1; k<dimension; k++) {
            if(cities[k].added) { // If City is already in tour
                continue;
            }else{
                // Calculate distance
                curr_dist = distance(&cities[tour[l-1]], &cities[k]);
                if(curr_dist<lowDist) {
                    lowDist = curr_dist;
                    tempInd = k;
                }
            }
        }
        tour.push_back(tempInd);
        dist.push_back(lowDist);
        cities[tempInd].added = true;
    }

    // Calculate distance from last City to first
    dist.push_back(distance(&cities[tour[0]], &cities[tour[tour.size()-1]]));

    getline(file, line);

    // Parsing items
    for(int i=0; i<itemNum; i++) {
        Item newI;
        getline(file, line);
        ss.str(line);
        ss >> newI.index;
        ss >> newI.profit;
        ss >> newI.weight;
        ss >> newI.city;
        items.push_back(newI);
        ss.clear();
    }

    // We are given the cities with 1 indexing so lets stick with that
    for (int i=0; i<(int)tour.size(); i++) {
        tour[i]++;
    }
}
