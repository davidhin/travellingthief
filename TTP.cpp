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

int main() {
  // **** 1. Parsing **** //

  // Getting file
  ifstream file("simple4_n6.ttp");

  string line; // Holds each line of the file

  for(int i=0; i<3; i++) {
    getline(file, line);
  }

  // Dimension (number of cities)
  string temp;
  int dimension;
  istringstream ss(line);
  ss >> temp >> dimension;
  ss.clear();

  // Number of items
  int itemNum;
  getline(file, line);
  ss.str(line);
  ss >> temp >> temp >> temp >> itemNum;
  ss.clear();

  // Getting knapsack capaCity
  int ksc;
  getline(file, line);
  ss.str(line);
  ss >> temp >> temp >> temp >> ksc;
  ss.clear();

  // Min speed
  float minS;
  getline(file, line);
  ss.str(line);
  ss >> temp >> temp >> minS;
  ss.clear();

  // Max speed
  float maxS;
  getline(file, line);
  ss.str(line);
  ss >> temp >> temp >> maxS;
  ss.clear();

  // Rent
  float rent;
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
  vector<int> dist;

  // Adding first City to the tour
  tour.push_back(0);
  dist.push_back(0);
  cities[0].added = true;

  // For storing current closest City details
  int lowDist;
  int tempInd = 1;
  // Calculating distance variables
  int temp1;
  int temp2;

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
  dist[dist.size()-1] = temp1;

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
  for (int i=0; i<tour.size(); i++) {
    tour[i]++;
  }

  // Create a vector of "distances to end of tour"
  vector<double> end_dists;
  for (int i=0; i<dist.size(); i++) {
    end_dists.push_back(accumulate(dist.begin()+i+1, dist.end(), 0));
  }

  // Printing
  cout << "Tour:\n" << tour << "\n\n";
  cout << "Distances:\n" << dist << "\n\n";
  cout << "Distance from city i to end of tour: \n" << end_dists << "\n\n";

  // Assign items to cities
  unordered_map<int, vector<Item*>> city_items;
  for (int i=0; i<items.size(); i++) {
    city_items[items[i].city].push_back(&items[i]);
  }

  // Print city/item map
  cout << "Mapped items to cities: " << "\n";
  for (auto i = city_items.begin(); i != city_items.end(); i++) {
    cout << "City " << i->first << " items: ";
    vector<int> temp;
    for (int j=0; j<i->second.size(); j++)
      temp.push_back(i->second[j]->index);
    cout << temp << "\n";
  }
  cout << "\n";

  // ** a. Compute score for each of the items ** //
  vector< pair<double, int> > score_items;

  for (int i=0; i<tour.size(); i++) {
    // Get list of items at city tour[i]
    unordered_map<int, vector<Item*>>::iterator it = city_items.find(tour[i]+1);
    if (it == city_items.end()) continue;

    // Pointer to list of items at the city tour[i]
    vector<Item*> *current_items = &it->second;

    // Add item to "full" item list
    for (int j=0; j<current_items->size(); j++) {
      Item* current_item = current_items->at(j);
      double score = current_item->profit / (current_item->weight * end_dists[j]);
      score_items.push_back(make_pair(score,current_item->index));
    }
  }

  // ** b. Sort the items of M in non-decreasing order of their scores ** //
  cout << "Item Queue: " << "\n";
  // My computer can't run this
  // sort(rbegin(score_items), rend(score_items));
  sort(score_items.begin(), score_items.end());
  cout << score_items << "\n";
}
