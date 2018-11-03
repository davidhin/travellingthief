// Algorithms where specified were based on those from:
// "Approximate Approaches to the Traveling Thief Problem" 2015, by Faulkner et. al.

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
#include <random>

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
    bool picked = false; // If item is in packing plan
};

// Print an item by printing its index
ostream &operator<<(ostream &os, Item* const &i) {
    return os << i->index;
}

// Print an item by printing its index
ostream &operator<<(ostream &os, Item const &i) {
    return os << "<id:" << i.index <<",city:" << i.city << ",weight:" << i.weight << "," <<"profit:" << i.profit << ">";
}

// Constants
int tau = 1 ; // Used to make packing plan quicker but more inaccurate

// Globals
int dimension; // Number of cities
double minS; // Minimum Speed
double maxS; // Maximum Speed
int ksc; // Knapsack Capacity
double W_curr; // Current weight of knapsack
double rent; // Cost of knapsack rent
int itemNum; // Number of items

vector<int> tour; // Store order of the cities in the tour
vector<double> dist; // Distance of city in tour to next city, starts with 0
vector<City> cities; // Parsed cities
vector<Item> items; // Parsed Items
vector<double> end_dists; // City's distance to the end of the tour
unordered_map<int, vector<Item*>> city_items; // Assigning items to cities

// Forward Declarations
double Z(vector<int>*, vector<double>*, vector<Item*>*); // Objective function / Z-score
vector< pair<double, int> > get_scores(double); // Returns items in "score" order
vector<Item*> packing_plan(); // Returns a packing plan

// Algorithms
void insertion(vector<int>*, bool);
void bitflip(vector<int>*, bool);
void pack_iterative(double c, double delta, int q, bool);
void initialise_plan(bool, double);
void two_opt(bool);
void two_opt_Z(bool output);
void lin_kern(bool, int);
void lin_kern_Z(bool output, int start);

// Helper functions
double distance(City*, City*); // Distance between 2 cities
double tour_dist(); // Return length of tour
int distCount(int, int); // Return distance from one city to another in the tour
void dist_end_tour(); // Setting city distance till end of tour & assigning items to cities
vector<Item*> get_items_in_city(int city); // Return vector of items from given city
void parsingFile(char* argv[]); // Parsing input file

int main(int argc, char* argv[]) {
    // Parsing the files
    parsingFile(argv);
    dist_end_tour();

    // Function plans for separate test cases
    if(dimension==51 && itemNum==50) { // Case 2
        insertion(&tour, false) ;
        for(int i=1; i<dimension-2; i++) {
            lin_kern_Z(true, i);
        }
        two_opt(false);
        pack_iterative(1, 10, 20, false);
        insertion(&tour, false) ;
        pack_iterative(1, 10, 20, false);
    }else if(dimension < 270) { // Case 1 and 3
        two_opt(false);
        pack_iterative(1, 10, 20, false);
        cout << "Done Pack Iterative" << "\n";
        insertion(&tour, false) ;
        two_opt_Z(true);
        pack_iterative(1, 10, 20, false);
        vector<Item*> plan = packing_plan();
        cout << "PackingPlan 1: " << Z(&tour, &dist, &plan) << "\n";
    } else { // Case 4 and 5
        two_opt(false);
        initialise_plan(false, -1.03125);
        two_opt(false);
        initialise_plan(false, -9);
    }

    // Further improvement functions
    for(int i=0; i<10; i++) {
        bitflip(&tour, false);

        if (dimension < 500) {
            if (dimension > 270 && i % 2 == 0) { continue; }
            insertion(&tour, false) ;
            two_opt_Z(false);
        }

        vector<Item*> plan = packing_plan();
        cout << i << ": " << Z(&tour, &dist, &plan) << "\n";

        // Output to file to ensure time out doesn't result in no output
        ofstream outputFile;
        outputFile.open("fnl_soln.ttp");
        outputFile << tour << "\n" ;
        outputFile << packing_plan() << "\n" ;
        outputFile.close();
    }

    // Output to file
    ofstream outputFile;
    outputFile.open("fnl_soln.ttp");
    outputFile << tour << "\n" ;
    outputFile << packing_plan() << "\n" ;
    outputFile.close();
}

// ** as specified from Faulkner et. al. 2015
// Objective function / Z-score
// Parameters:
// tour_in: a TSP tour
// dist: distances from each city to next in tour_in
// plan: packing plan for tour_in
double Z(vector<int>* tour_in, vector<double>* dist, vector<Item*>* plan) {
    long collected = 0;
    double ret = 0;
    double v = ( maxS - minS ) / ksc;

    // Calculates objective score based on formula:
    // (sum of all packed items) - rent * (total travelling time)
    for (int i = 0; i < (int)(*tour_in).size(); i++) {
        vector<Item*> current_items = get_items_in_city((*tour_in)[i]);
        for (int j=0; j<(int)current_items.size(); j++) {
            Item* current_item = current_items.at(j);
            if(current_item->picked==true) {
                collected += current_item->weight;
                ret += current_item->profit;
            }
        }
        ret -= ceil((*dist)[i+1]) * rent / (maxS - v * collected);
    }

    return ret;
}

// Item score calculation
// Parameters:
// alpha: exponent for score calculation
vector< pair<double, int> > get_scores(double alpha) {
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

    // Sort the items of M in non-decreasing order of their scores
    sort(score_items.begin(), score_items.end());

    return score_items;
}

// Return the current packing plan as a vector
vector<Item*> packing_plan() {
    vector<Item*> plan;
    for(int i=0; i<(int)items.size(); i++) {
        if (items[i].picked) {
            plan.push_back(&items[i]);
        }
    }
    return plan;
}

// ** as specified from Faulkner et. al. 2015
// Insertion - Aim: to put "cities where a lot of items are stolen" near
// the end of the tour. This is done through brute force iteration.
// Parameters:
// tour_in: a TSP tour
// output: true for couts to print
void insertion(vector<int>* tour_in, bool output) {
    vector<Item*> plan = packing_plan();

    for (int i = tour_in->size()-1; i>=1; i--) {
        double Z_best = Z(tour_in, &dist, &plan);
        if (output) cout << i << " Z_best: " << Z_best << "\n";

        vector<int> tour_starstar = *tour_in;
        vector<int> tour_new = *tour_in;

        for (int j=i; j>=2; j-=1) {
            // Insert city between xb and xb-1 in tour_new
            swap(tour_new[j], tour_new[j-1]);
            // Compute Z*
            for (int k=max(j-1,1); k<min(j+2, (int)dist.size()-1); k++) {
                dist[k] = distance(&cities[tour_new[k-1]-1], &cities[tour_new[k]-1]);
            }
            dist[dist.size()-1] = distance(&cities[tour_new[tour_new.size()-1]-1], &cities[tour_new[0]-1]);
            double Z_curr = Z(&tour_new, &dist, &plan);
            // Update if better
            if (Z_curr > Z_best) {
                Z_best = Z_curr;
                tour_starstar = tour_new;
                if (output) cout << "updated" << "\n";
            }
        }

        // Update actual tour to tour_starstar
        *tour_in = tour_starstar;
        for (int k=1; k<(int)dist.size()-1; k++) {
            dist[k] = distance(&cities[(*tour_in)[k-1]-1], &cities[(*tour_in)[k]-1]);
        }
        dist[dist.size()-1] = distance(&cities[(*tour_in)[(*tour_in).size()-1]-1], &cities[(*tour_in)[0]-1]);
    }
    if (output) cout << "\n";
}

// ** as specified from Faulkner et. al. 2015
// Bitflip Algorithm - inverts item picked state, keeps change if improved Z-score
// Parameters:
// tour_in: a TSP tour
// output: true for couts to print
void bitflip(vector<int>* tour_in, bool output) {
    vector<Item*> plan_init = packing_plan();
    double Z_score = Z(tour_in, &dist, &plan_init);
    double Z_flipped;

    for (int i=0; i<(int)items.size(); i++) {
        // Invert picked state
        if (items[i].picked) { // If already picked
            items[i].picked = false;
        } else { // If not in tour
            items[i].picked = true;
            if(W_curr+items[i].weight>ksc) {
                items[i].picked = false;
                continue;
            }
        }

        vector<Item*> plan = packing_plan();
        Z_flipped = Z(tour_in, &dist, &plan);
        if (output) cout<< i << " " << Z_flipped<<" flipped"<<endl;
        if (Z_flipped<Z_score) { // Worse Z-score
            if (items[i].picked) {
                items[i].picked = false;
            } else {
                items[i].picked = true;
            }
        } else { // Improved Z-score
            Z_score = Z_flipped;
            if (items[i].picked) {
                W_curr += items[i].weight;
            } else {
                W_curr -= items[i].weight;
            }
        }
    }
}

// ** as specified from Faulkner et. al. 2015
// Search to find best parameters for initialise_plan
// Parameters:
// c: starting alpha value
// delta: size of "jumps"
// q: maximum number of iterations
// output: true for couts to print
void pack_iterative(double c, double delta, int q, bool output) {
    // Checking Z score for different parameters
    initialise_plan(false, c - delta);
    vector<Item*> P_l = packing_plan();
    double Z_l = Z(&tour, &dist, &P_l);
    if (output) cout << "Z_l start: " << Z_l << "\n";

    initialise_plan(false, c);
    vector<Item*> P_m = packing_plan();
    double Z_m = Z(&tour, &dist, &P_m);
    if (output) cout << "Z_m start: " << Z_m << "\n";

    initialise_plan(false, c + delta);
    vector<Item*> P_r = packing_plan();
    double Z_r = Z(&tour, &dist, &P_r);
    if (output) cout << "Z_r start: " << Z_r << "\n";

    vector<Item*> P_best;

    int i = 0;
    while (i <= q) {
        if (output) cout << Z_l << " " << Z_m << " " << Z_r << "\n";

        // Selecting best parameters
        if (Z_l > Z_m && Z_r > Z_m) {
            if (Z_l > Z_r) {
                Z_m = Z_l, c -= delta, P_best = P_l;
            } else {
                Z_m = Z_r, c += delta, P_best = P_r;
            }
        }
        else if (Z_l > Z_m) {
            Z_m = Z_l, c = c - delta, P_best = P_l;
        }
        else if (Z_r > Z_m) {
            Z_m = Z_r, c = c + delta, P_best = P_r;
        }
        delta /= 2;

        // Setting new search space
        initialise_plan(false, c - delta);
        P_l = packing_plan();
        Z_l = Z(&tour, &dist, &P_l);

        initialise_plan(false, c + delta);
        P_r = packing_plan();
        Z_r = Z(&tour, &dist, &P_r);

        // Exit early if difference is negligible
        if (abs(Z_l - Z_m) < 0.1 && abs(Z_r - Z_m) < 0.1) {
            return;
        }

        i++;
        if (output) cout << c << " " << Z_m << "\n";
    }
}

// ** as specified from Faulkner et. al. 2015
// Generate a packing plan for a given tour
// Parameters:
// output: true for couts to print
// alpha: exponent for score calculation
void initialise_plan(bool output, double alpha) {
    // Packing Plan Algorithm Starts Here
    vector< pair<double, int> > score_items = get_scores(alpha); // Get items ordered by score

    // Reset item picked state
    for (int i=0; i<itemNum; i++)
        items[i].picked = false;

    W_curr = 0; // Initialise current weight and weight capacity
    int mu = floor(itemNum / tau) ; // Set the frequency of item adding
    double Z_best = -0xfffffff ; // initialise best objective value to -inf
    int k = 0; //, k_last=0; // used to iterate through score_items
    vector<Item*> plan;
    // Loop to add items if they improve the Z-score
    while (W_curr < ksc && mu > 1 && k < itemNum) {
        Item* kth_item = &items[score_items[k].second-1];

        if (W_curr + kth_item->weight <= ksc) {
            kth_item->picked = true;
            plan.push_back(kth_item);
            W_curr += kth_item->weight;
            if (k % 1 == 0) {
                double Z_curr = Z(&tour, &dist, &plan) ;
                if (Z_curr < Z_best) { // Item decreases Z-score (revert)
                    kth_item->picked = false;
                    W_curr -= kth_item->weight ;
                    plan.erase(plan.end()-1);
                }
                else { // Item improves Z-score
                    Z_best = Z_curr ;
                }
            }
        }
        k++;
        if (output) cout << k << " " << Z_best<<endl;
    }
}

// Two-opt Algorithm - swaps every pair of nodes to check if new connection shortens tour distance
// Parameters:
// output: true for couts to print
void two_opt(bool output) {
    vector<int> old_tour = tour;
    vector<double> old_dist = dist;
    int two_opt_dist;
    int temp_dist;

    for(int two_opt_start=1; two_opt_start<dimension-1; two_opt_start+=1) { // First node of 2-opt
        for(int two_opt_end=two_opt_start+1; two_opt_end<dimension; two_opt_end+=1) { // Second node of 2-opt
            two_opt_dist = distCount(two_opt_start-1, two_opt_end+1);
            reverse(tour.begin()+two_opt_start, tour.begin()+two_opt_end); // Creating new edge

            // Recalculating distances
            for (int i=two_opt_start-1; i<two_opt_end; i++) {
                dist[i+1] = distance(&cities[tour[i+1]-1], &cities[tour[i]-1]);
            }

            // Calculate distance from last City to first
            dist[dist.size()-1] = distance(&cities[tour[0]-1], &cities[tour[tour.size()-1]-1]);
            temp_dist = distCount(two_opt_start-1, two_opt_end+1);

            // Keeping change if new tour is shorter
            if(temp_dist < two_opt_dist) { // Shorter distance
                two_opt_dist = temp_dist;
                old_tour = tour;
                old_dist = dist;
                if (output) cout<<"dist "<<temp_dist<<endl;
            }else { // Longer distance
                tour = old_tour;
                dist = old_dist;
            }
        }
    }
}

// Two-opt-Z - swaps every pair of nodes to check if new connection increases the Z-score
// Parameters:
// output: true for couts to print
void two_opt_Z(bool output) {
    vector<int> old_tour = tour;
    vector<double> old_dist = dist;
    int to_Z;
    int to_Z_temp;

    vector<Item*> plan = packing_plan();
    to_Z = Z(&tour, &dist, &plan);

    for(int two_opt_start=1; two_opt_start<dimension-1; two_opt_start+=1) { // First node of 2-opt
        for(int two_opt_end=two_opt_start+1; two_opt_end<dimension; two_opt_end+=1) { // Second node of 2-opt
            reverse(tour.begin()+two_opt_start, tour.begin()+two_opt_end); // Creating new edge

            // Recalculating distances
            for (int i=two_opt_start-1; i<two_opt_end; i++) {
                dist[i+1] = distance(&cities[tour[i+1]-1], &cities[tour[i]-1]);
            }

            // Calculate distance from last City to first
            dist[dist.size()-1] = distance(&cities[tour[0]-1], &cities[tour[tour.size()-1]-1]);
            to_Z_temp = Z(&tour, &dist, &plan);

            if(to_Z_temp > to_Z) { // Higher Z-score
                to_Z = to_Z_temp;
                old_tour = tour;
                old_dist = dist;
                if (output) cout<<"Z-score: "<<to_Z_temp<<endl;
            }else { // Lower Z-score
                tour = old_tour;
                dist = old_dist;
            }
        }
    }
}

// Lin Kernighan algorithm, shortens tour based on distance
// Parameters:
// output: true for couts to print
// start: starting city for Lin Kernighan
void lin_kern(bool output, int start) {
    vector<int> old_tour = tour;
    vector<double> old_dist = dist;
    int v1 = start;

    double dist1_orig;
    double dist2_orig;
    double dist1;
    double dist2;

    // Loop through possible nodes to swap
    for(int two_opt_end=v1+1; two_opt_end<dimension-1; two_opt_end+=1) {
        // Calculating new and old distances
        dist1_orig = distance(&cities[tour[v1-1]-1], &cities[tour[v1]-1]);
        dist2_orig = distance(&cities[tour[two_opt_end-1]-1], &cities[tour[two_opt_end]-1]);
        dist1 = distance(&cities[tour[v1]-1], &cities[tour[two_opt_end]-1]);
        dist2 = distance(&cities[tour[v1-1]-1], &cities[tour[two_opt_end-1]-1]);

        dist1_orig = dist1_orig + dist2_orig;
        dist1 = dist1 + dist2;

        if(dist1<dist1_orig) { // If swapping would result in a shorter tour
            reverse(tour.begin()+v1, tour.begin()+two_opt_end);

            // Recalculating distances
            for (int i=v1-1; i<two_opt_end; i++) {
                dist[i+1] = distance(&cities[tour[i+1]-1], &cities[tour[i]-1]);
            }

            // Calculate distance from last City to first
            dist[dist.size()-1] = distance(&cities[tour[0]-1], &cities[tour[tour.size()-1]-1]);

            // Set new tour
            old_tour = tour;
            old_dist = dist;

            // Move onto next opt
            v1 = two_opt_end;
            two_opt_end++;
        }else { // If longer tour
            tour = old_tour;
            dist = old_dist;
        }

    }
    // Set new tour
    tour = old_tour;
    dist = old_dist;
}

// Lin Kernighan algorithm, shortens tour based on Z-score
// Parameters:
// output: true for couts to print
// start: starting city for Lin Kernighan
void lin_kern_Z(bool output, int start) {
    vector<int> old_tour = tour;
    vector<double> old_dist = dist;
    int to_Z;
    int to_Z_temp;
    int v1 = start;

    vector<Item*> plan = packing_plan();
    to_Z = Z(&tour, &dist, &plan);

    // Loop through possible nodes to swap
    for(int two_opt_end=v1+1; two_opt_end<dimension; two_opt_end+=1) {
        // Reverse middle segment
        reverse(tour.begin()+v1, tour.begin()+two_opt_end);

        // Recalculating distances
        for (int i=v1-1; i<two_opt_end; i++) {
            dist[i+1] = distance(&cities[tour[i+1]-1], &cities[tour[i]-1]);
        }

        // Calculate distance from last City to first
        dist[dist.size()-1] = distance(&cities[tour[0]-1], &cities[tour[tour.size()-1]-1]);
        to_Z_temp = Z(&tour, &dist, &plan);

        // Keeping change if new tour is shorter
        if(to_Z_temp > to_Z) { // Higher Z score
            to_Z = to_Z_temp;
            // Set new tour
            old_tour = tour;
            old_dist = dist;
            if (output) cout<<"Z-score: "<<to_Z_temp<<endl;
            if (output) cout<<"start: "<<v1<<endl;
            if (output) cout<<"end: "<<two_opt_end<<endl;
            // Move onto next opt
            v1 = two_opt_end;
            two_opt_end++;
        }else {  // If worse Z score
            tour = old_tour;
            dist = old_dist;
        }
    }
    // Set new tour
    tour = old_tour;
    dist = old_dist;
}

// Distance between two cities
double distance(City* a, City* b) {
    return sqrt(pow(a->xCoord - b->xCoord, 2) + pow(a->yCoord - b->yCoord, 2));
}

// Calculating distance of tour
double tour_dist() {
    double temp_dist = 0;
    for(int i=0; i<dist.size(); i++) {
        temp_dist += dist[i];
    }
    return temp_dist;
}

// Calculating distance from one city to another in tour
// Parameters:
// s: first city
// e: end city
int distCount(int s, int e) {
    double total = 0;
    // Adding distances along tour
    for(int i=s; i<e; i++) {
        total += dist[i];
    }
    return total;
}

// Setting city distance till end of tour & assigning items to cities
void dist_end_tour() {
    end_dists = {};
    city_items = {};

    // Create a vector of "distances to end of tour"
    for (int i=0; i<(int)dist.size(); i++) {
        end_dists.push_back(accumulate(dist.begin()+i+1, dist.end(), 0.0));
    }

    // Assign items to cities
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
}

// Return vector of items in a city
// Parameters:
// city: index of city
vector<Item*> get_items_in_city(int city) {
    unordered_map<int, vector<Item*>>::iterator it = city_items.find(city);
    if (it == city_items.end()) return {}; // If city not in map

    // Pointer to list of items at the city at index city
    return it->second;
}

// Parsing input file
// Parameters:
// argv: main input parameters
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
        string doubleCoor;
        ss >> doubleCoor;
        newC.xCoord = stod(doubleCoor);
        ss >> doubleCoor;
        newC.yCoord = stod(doubleCoor);
        newC.added = false;
        cities.push_back(newC);
        ss.clear();
    }

    // Nearest Neighbours Algorithm
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
