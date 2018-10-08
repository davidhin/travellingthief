#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <limits.h>
#include <math.h> //sqrt & pow

using namespace std;

// Encapsulates a single city
struct city {
  int index;
  int xCoord;
  int yCoord;
  bool added; // If city has been added to
};

// Encapsulates a single item
struct item {
  int index;
  int profit;
  int weight;
  int city;
};

int main() {
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

  // Getting knapsack capacity
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

  // Parsing City coordinates into vector of city structs
  vector<city> cities;
  for(int i=0; i<dimension; i++) {
    city newC;
    getline(file, line);
    ss.str(line);
    ss >> newC.index;
    ss >> newC.xCoord;
    ss >> newC.yCoord;
    newC.added = false;
    cities.push_back(newC);
    ss.clear();
  }

  // Store order of the cities in the tour and the distance to the next city
  vector<int> tour;
  vector<int> dist;

  // Adding first city to the tour
  tour.push_back(0);
  dist.push_back(0);
  cities[0].added = true;

  // For storing current closest city details
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
      if(cities[k].added) { // If city is already in tour
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
    dist.push_back(0);
    cities[tempInd].added = true;
    dist[l-1] = lowDist;
  }

  // Calculate distance from last city to first
  temp1 = cities[tour[0]].xCoord - cities[tour[tour.size()-1]].xCoord;
  temp2 = cities[tour[0]].yCoord - cities[tour[tour.size()-1]].yCoord;
  temp1 = pow(temp1, 2);
  temp2 = pow(temp2, 2);
  temp1 = temp1 + temp2;
  temp1 = sqrt(temp1);
  dist[dist.size()-1] = temp1;

  getline(file, line);

  // Parsing Items
  vector<item> items;
  for(int i=0; i<itemNum; i++) {
    item newI;
    getline(file, line);
    ss.str(line);
    ss >> newI.index;
    ss >> newI.profit;
    ss >> newI.weight;
    ss >> newI.city;
    items.push_back(newI);
    ss.clear();
  }

}
