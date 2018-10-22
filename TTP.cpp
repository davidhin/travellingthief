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
};

// Print ItemPointer as the index
//ostream &operator<<(ostream &os, Item* const &i) {
//    return os << "<id:" << i->index <<",city:" << i->city << ",weight:" << i->weight << "," <<"profit:" << i->profit << ">";
//}

// Print with only the index (used for submission)
ostream &operator<<(ostream &os, Item* const &i) {
    return os << i->index;
}

// Output Enabled/disabled
bool output = false;

// Constants
int tau = 1 ; // TODO: Find a better value to set this to

// Globals
int dimension;
double minS;
double maxS;
int ksc;
double rent;
int itemNum;

// Get total weight of all items in city, using packingplan vector
// TODO: Find a more efficient way of doing this
double city_weight(int city, vector<Item*>* packing_plan)
{
    int ret = 0;
    for (int i = 0; i < (int)packing_plan->size(); i++)
        if ((*packing_plan)[i]->city == city)
            ret += (*packing_plan)[i]->weight ;
    return ret ;
}

// Objective function
double Z(vector<int>* tour, vector<double>* dist, vector<Item*>* packing_plan) {
    // Sum of all packed items
    double minuend = 0;
    for (int i = 0; i<(int)packing_plan->size(); i++)
        minuend += (*packing_plan)[i]->profit ;

    // Constant V
    double v = ( maxS - minS ) / ksc ;

    // Amount that the thief pays for the knapsack rent
    // i.e. The total travelling time along tour * Rent
    // TODO: Check if the city weights actually works
    double subtrahend = 0;
    for(int i = 1; i < (int)(*dist).size(); i++)
        subtrahend += (double)(*dist)[i] / (maxS - v*city_weight((*tour)[i], packing_plan)) ;

    return minuend - rent*(subtrahend);
}

int main(int argc, char* argv[]) {
    // **** 1. Parsing **** //

    ifstream file(argv[1]);
    string line; // Holds each line of the file

    char delim = '\n';

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
    vector<City> cities;
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

    // Store order of the cities in the tour and the distance to the next City
    vector<int> tour;
    vector<double> dist;

    // Adding first City to the tour
    tour.push_back(0);
    dist.push_back(0);
    cities[0].added = true;

    // For storing current closest City details
    double lowDist;
    int tempInd = 1;
    // Calculating distance variables
    double temp1;
    double temp2;

    // Distance formula
    // sqrt( (x1-x2)^2 + (y1-y2)^2 )

    // Getting tour
    for(int l=1; l<dimension; l++) {
        lowDist = INT_MAX; // reset distance
        for(int k=1; k<dimension; k++) {
            if(cities[k].added) { // If City is already in tour
                continue;
            }else{
                // Calculate distance
                temp1 = cities[tour[l-1]].xCoord - cities[k].xCoord;
                temp2 = cities[tour[l-1]].yCoord - cities[k].yCoord;
                temp1 = pow(temp1, 2);
                temp2 = pow(temp2, 2);
                temp1 = temp1 + temp2;
                temp1 = sqrt(temp1);
                if(temp1<lowDist) {
                    lowDist = temp1;
                    tempInd = k;
                }
            }
        }
        tour.push_back(tempInd);
        dist.push_back(lowDist);
        cities[tempInd].added = true;
    }

    // Calculate distance from last City to first
    temp1 = cities[tour[0]].xCoord - cities[tour[tour.size()-1]].xCoord;
    temp2 = cities[tour[0]].yCoord - cities[tour[tour.size()-1]].yCoord;
    temp1 = pow(temp1, 2);
    temp2 = pow(temp2, 2);
    temp1 = temp1 + temp2;
    temp1 = sqrt(temp1);
    dist.push_back(temp1); // TODO: Is this correct? We want to push back, not replace the last distance
    // OLD dist[dist.size()-1] = temp1;

    getline(file, line);

    // Parsing items
    vector<Item> items;
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

    // **** 3. Algorithm 1 **** //

    // ** 0. Setup data structures ** //
    // TODO: I tried to do this whole thing without modifying
    // the original code above, meaning some parts are a bit dodgy
    // Calculate distance from City i to end

    //TODO: We are given the cities with 1 indexing so lets
    //stick with that
    for (int i=0; i<(int)tour.size(); i++) {
        tour[i]++;
    }

    // Create a vector of "distances to end of tour"
    vector<double> end_dists;
    for (int i=0; i<(int)dist.size(); i++) {
        end_dists.push_back(accumulate(dist.begin()+i+1, dist.end(), 0.0));
    }

    // Printing
    if (output) cout << "Tour:\n" << tour << "\n\n";
    if (output) cout << "Distances:\n" << dist << "\n\n";
    if (output) cout << "Distance from city i to end of tour: \n" << end_dists << "\n\n";

    // Assign items to cities
    unordered_map<int, vector<Item*>> city_items;
    for (int i=0; i<(int)items.size(); i++) {
        city_items[items[i].city].push_back(&items[i]);
    }

    // Print city/item map
    if (output) cout << "Mapped items to cities: " << "\n";
    for (auto i = city_items.begin(); i != city_items.end(); i++) {
        if (output) cout << "City " << i->first << " items: ";
        vector<int> temp;
        for (int j=0; j<(int)i->second.size(); j++)
            temp.push_back(i->second[j]->index);
        if (output) cout << temp << "\n";
    }
    if (output) cout << "\n";

    // ** a. Compute score for each of the items ** //
    vector< pair<double, Item*> > score_items;

    for (int i=0; i<(int)tour.size(); i++) {
        // Get list of items at city tour[i]
        unordered_map<int, vector<Item*>>::iterator it = city_items.find(tour[i]+1);
        if (it == city_items.end()) continue;

        // Pointer to list of items at the city tour[i]
        vector<Item*> *current_items = &it->second;

        // Add item to "full" item list
        for (int j=0; j<(int)current_items->size(); j++) {
            Item* current_item = current_items->at(j);
            double score = current_item->profit / (current_item->weight * end_dists[j]);
            score_items.push_back(make_pair(score,current_item));
        }
    }

    // ** b. Sort the items of M in non-decreasing order of their scores ** //
    if (output) cout << "Item Queue: " << "\n";
    // My computer can't run this
    // sort(rbegin(score_items), rend(score_items));
    sort(score_items.begin(), score_items.end());
    if (output) cout << score_items << "\n\n";

    // ** c. set the frequency mu = floor(m / tau)
    int mu = floor(itemNum / tau) ;
    if (output) cout << "We set mu = floor(itemNum / tau)=" << mu << " where tau=" << tau << " is a given constant.\n" ;

    // ** d. Initialise the current packing plan P and current weight of bag = 0
    vector<Item*> P_curr ;
    double W_curr = 0, W_cap = ksc ;
    if (output) cout << "P_curr = " << P_curr << " is the initialised packing plan (current)" << "\n" ;
    if (output) cout << "W_curr = " << W_curr << " is the current weight of the packing plan, and" << "\n"
             << "W_cap = " << W_cap << " is the capacity of the bag. \n\n";

    // ** e. Set the best packing plan P_best = empty
    vector<Item*> P_best ;
    if (output) cout << "We initialise an array to hold our best packing plan\nP_best = " << P_best << "\n\n" ;

    // ** f. Set the best objective value Z_best = -inf
    double Z_best = -0xfffffff ;
    if (output) cout << "Z_best is the best objective value (so far)." << "\n";
    if (output) cout << "Z_best = " << Z_best << "\n\n" ;

    // ** g. Set counters to iterate through the score_items
    int k = 0, k_last = 0;
    if (output) cout << "Counters k = k-th item in score_items, it starts at 0.\n" <<
    "k_last also starts at 0, and is the last item that we added to the plan.\n\n" ;

    if (output) cout << "Objective function score: \n" ;
    if (output) cout << Z_best << "\n";
    // ** h. Iterate through items k in the item queue
    while (W_curr < W_cap && mu > 1 && k < itemNum) {
        Item* kth_item = score_items[k].second;
        if (W_curr + kth_item->weight <= W_cap) {
            P_curr.push_back(kth_item);
            W_curr += kth_item->weight;
            if (true) { // TODO We only put the mu limiter in once we use big sets
            // if (k % mu == 0) {
                double Z_curr = Z(&tour, &dist, &P_curr) ;
                if (Z_curr < Z_best) {
                    P_curr = P_best ;
                    k = k_last ;
                    mu /= 2 ; // Truncate / floor
                }
                else {
                    P_best = P_curr ;
                    k_last = k ;
                    Z_best = Z_curr ;
                }
            }
        }
        if (output) cout << Z_best << "\n";
        k++;
    }

    if (output)  cout << "\nFinal Packing Plan: \n";
    ofstream outputFile;
    outputFile.open("fnl_soln.ttp");
    outputFile << tour << "\n" ;
    outputFile << P_best << "\n" ;
    outputFile.close();
    
    // Old couts - still use for testing
    // cout << tour << "\n" ;
    // cout << P_best << "\n" ;
}
