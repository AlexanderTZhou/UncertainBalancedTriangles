#define _CRT_SECURE_NO_WARNINGS
#define ld long double

#include<iostream>
#include<fstream>
#include<string>
#include<cstring>
#include<stdio.h>
#include<stdlib.h>
#include<vector>
#include<unordered_map>
#include<random>
#include<utility>
#include<functional>
#include<cstdio>
#include<ctime>
#include<chrono>
#include<algorithm>
#include<cmath>
#include<numeric>

using namespace std;

//Writes information to a file
void writeToFile(string s, char* filename) {
	ofstream writeFile;
	writeFile.open(filename, ios::app);
	if (writeFile.is_open()) {
		writeFile << s << "\n";
		writeFile.flush();
	} 
	writeFile.close();
}



//Classes
class Node {
    public:
        int id;
        vector<pair<int, double> > neighbours; //neighbours with <ID, prob/aev>
        unordered_map<int, double> neighbourMap; //Hashing optimization for performance (O(|E|) complexity)
};

class BasicEdge {
    public:
        int u1;
        int u2;
        double prob;
};

class AdvEdge {
    public:
        int u1;
        int u2;
        double prob;
        double aev;
};


//Globals
double t;
vector<int> vIDs;
vector<BasicEdge> basicEdgeList;
vector<AdvEdge> advEdgeList;
unordered_map<int, Node> nodeMap;

/***
* Calculates Absolute Edge Value
*/
double calcAEV(double val) {
    return fabs(val - 0.5);
}

double calcBalProb(double e1, double e2, double e3) {
    double o1 = 1.0;
    return (e1 * e2 * e3) + ((o1 - e1) * (o1 - e2) * e3) + ((o1 - e1) * e2 * (o1 - e3)) + (e1 * (o1 - e2) * (o1 - e3));   
}

double calcUnbalProb(double e1, double e2, double e3) {
    return (1.0 - calcBalProb(e1, e2, e3));
}

double calcAEVBound(double e) {
    double upper = (2.0 * e * e) - (2.0 * e) + t;
    double lower = (4.0 * e * e) - (4.0 * e) + 1.0;
    return upper/lower;
}

/***
* nCr function
*/
int nChoose2(int n) {
	return (n - 1) * n / 2;
}

/***
* Randomly samples an integer from a list
*/
int rSample(int size, mt19937 &mtGen) {
	return mtGen() % size;
}

/***
* Randomly samples an integer excluding a previously sampled number (i.e. without replacement)
*/
int rSample(int size, unordered_map<int, bool> &map, mt19937 &mtGen) {
	int randVal = mtGen() % size;
	unordered_map<int, bool>::iterator mIt = map.find(randVal);
	if (mIt == map.end()) {
		map.insert(make_pair(randVal, true));
		return randVal;
	} else {
		return rSample(size, map, mtGen);	
	}
}

/***
* Structs used for comparitors
*/
struct NodeOrderComp {
    bool operator()(const Node& n1, const Node& n2) {
        if (n1.neighbours.size() == n2.neighbours.size()) {
            return (n1.id < n2.id);
        }
        return (n1.neighbours.size() < n2.neighbours.size());
    }
    
    bool operator()(const int& n1, const int& n2) {
        unordered_map<int, Node>::iterator nIt1 = nodeMap.find(n1);
        unordered_map<int, Node>::iterator nIt2 = nodeMap.find(n2);
        if (nIt1->second.neighbours.size() == nIt2->second.neighbours.size()) {
            return nIt1->second.id < nIt2->second.id;
        }
        return (nIt1->second.neighbours.size() < nIt2->second.neighbours.size());
    }
    
    bool operator()(const pair<int, double>& n1, const pair<int, double>& n2) {
        unordered_map<int, Node>::iterator nIt1 = nodeMap.find(n1.first);
        unordered_map<int, Node>::iterator nIt2 = nodeMap.find(n2.first);
        if (nIt1->second.neighbours.size() == nIt2->second.neighbours.size()) {
            return nIt1->second.id < nIt2->second.id;
        }
        return (nIt1->second.neighbours.size() < nIt2->second.neighbours.size());
    }
};

struct AbsoluteEdgeOrderComp {
    bool operator()(const AdvEdge& e1, const AdvEdge& e2) {
        if (e1.aev == e2.aev) {
            if (e1.u1 == e2.u2) {
                return (e1.u2 < e2.u2);
            }
            return (e1.u1 < e2.u1);
        }
        return (e1.aev > e2.aev);
    }
};


//Sorts by AEV, assuming the insertion is a neighbour list
struct ProbToAEVComp {
    bool operator()(const pair<int, double>& n1, const pair<int, double>& n2) {
        double n1AEV = calcAEV(n1.second);
        double n2AEV = calcAEV(n2.second);
        if (n1AEV == n2AEV) {
            return n1.first < n2.first;
        }
        return (n1AEV > n2AEV);
    }
};

/***
* Checks Absolute Edge Order
* returns True if e1 before e2
*/
bool checkAEOrder(AdvEdge& e1, pair<int, double> e2) {
    double e2AEV = calcAEV(e2.second);
    if (e1.aev == e2AEV) {
        if (e1.u1 == e2.first) {
            return e1.u2 < e2.first;
        }
        return e1.u1 < e2.first;
    }
    return (e1.aev > e2AEV);
}

/***
* Checks Absolute Edge Order
* returns True if e1 before e2
*/
bool checkAEOrder(AdvEdge& e1, AdvEdge& e2) {
    if (e1.aev == e2.aev) {
        if (e1.u1 == e2.u1) {
            return e1.u2 < e2.u2;
        }
        return e1.u1 < e2.u1;
    }
    return (e1.aev > e2.aev);
}

/***
* Greather than binary search for absolute edge order
*/
int advEdgeOrderGTBS(vector<pair<int, double> >& neighbourList, int start, int end, AdvEdge& target) {
    
    if (start > end) {
        return -1;
    }
    
    if (start == end) {
        if (checkAEOrder(target, neighbourList[start])) {
            if (start == 0) {
                return start;
            }
            if (!checkAEOrder(target, neighbourList[start - 1])) {
                return start;
            }
        }
        return -1;
    }
    
    int mid = (end - start)/2 + start;
    if (checkAEOrder(target, neighbourList[mid])) {
        //Target is before Mid
        if (!checkAEOrder(target, neighbourList[mid - 1])) {
            return mid;
        } else {
            return advEdgeOrderGTBS(neighbourList, start, mid - 1, target);
        }
    } else {
        //Target is after Mid
        return advEdgeOrderGTBS(neighbourList, mid + 1, end, target);
    }
}

/***
* Returns True if n1 is before n2
* Returns False if n1 is after n2
*/
bool checkNodeOrder(Node& n1, int n2) {
    unordered_map<int, Node>::iterator nIt2 = nodeMap.find(n2);
    if (n1.neighbours.size() == nIt2->second.neighbours.size()) {
        return (n1.id < nIt2->second.id);
    }
    return (n1.neighbours.size() < nIt2->second.neighbours.size());
}

/***
* Finds the position of the node in the neighbourlist which starts the sequence of neighbours
* with Node Order after the target node.
* Returns -1 if no such node exists.
*/
int nodeOrderGTBS(vector<pair<int, double> >& neighbourList, int start, int end, Node& target) {
//    cout << mid << "\n";
    if (start > end) {
        return -1;
    }
    
    if (start == end) {
        if (checkNodeOrder(target, neighbourList[start].first)) {
            if (start == 0) {
                return start;
            }
            if (!checkNodeOrder(target, neighbourList[start - 1].first)) {
                return start;
            }
        }
        return -1;
    }
    
    int mid = (end - start)/2 + start;
    
    if (checkNodeOrder(target, neighbourList[mid].first)) {
        //Target is before Mid, go left
        if (mid == 0) {
            return mid;
        }
        if (!checkNodeOrder(target, neighbourList[mid - 1].first)) {
            return mid;
        } else {
            return nodeOrderGTBS(neighbourList, start, mid - 1, target);
        }
    } else {
        //Target is after Mid, go right
        return nodeOrderGTBS(neighbourList, mid + 1, end, target);
    }
}

/***
* Improved local counting for a vertex
*/
pair<double, double> vImprLocal(int u) {
    unordered_map<int, Node>::iterator uIt = nodeMap.find(u);
    unordered_map<int, Node>::iterator vIt;
    unordered_map<int, double>::iterator eIt;
    int iU = 0;
    int iV, iU2;
    double balCount = 0.0;
    double unbalCount = 0.0;
    double triBalProb, aevBound;
    double tAEV = calcAEV(t);
    AdvEdge aEdge;
    for (; iU < uIt->second.neighbours.size(); iU++) {
        aEdge.u1 = min(u, uIt->second.neighbours[iU].first);
        aEdge.u2 = max(u, uIt->second.neighbours[iU].first);
        aEdge.prob = uIt->second.neighbours[iU].second;
        aEdge.aev = calcAEVBound(uIt->second.neighbours[iU].second);
        if (aEdge.aev < tAEV) {
            break;
        }
        vIt = nodeMap.find(uIt->second.neighbours[iU].first);
        aevBound = calcAEV(aEdge.aev);
        if (uIt->second.neighbours.size() >= vIt->second.neighbours.size()) {
            //vIt is smaller
            iV = advEdgeOrderGTBS(vIt->second.neighbours, 0, vIt->second.neighbours.size() - 1, aEdge);
            for (; iV < vIt->second.neighbours.size(); iV++) {
                if (calcAEV(vIt->second.neighbours[iV].second) < tAEV) {
                    break;
                }
                eIt = uIt->second.neighbourMap.find(vIt->second.neighbours[iV].first);
                if (eIt != uIt->second.neighbourMap.end()) {
                    triBalProb = calcBalProb(uIt->second.neighbours[iU].second, vIt->second.neighbours[iV].second, eIt->second);
                    if (triBalProb >= t) {
                        balCount = balCount + 1.0;
                    } else if ((1.0 - triBalProb) > t) {
                        unbalCount = unbalCount + 1.0;
                    }
                }
            }
        } else {
            //uIt is smaller
            iU2 = advEdgeOrderGTBS(uIt->second.neighbours, 0, uIt->second.neighbours.size() - 1, aEdge);
            for (; iU2 < uIt->second.neighbours.size(); iU2++) {
                if (calcAEV(uIt->second.neighbours[iU2].second) < aevBound) {
                    break;
                }
                eIt = vIt->second.neighbourMap.find(uIt->second.neighbours[iU2].first);
                if (eIt != vIt->second.neighbourMap.end()) {
                    triBalProb = calcBalProb(uIt->second.neighbours[iU].second, uIt->second.neighbours[iU2].second, eIt->second);
                    if (triBalProb >= t) {
                        balCount = balCount + 1.0;
                    } else if ((1.0 - triBalProb) > t) {
                        unbalCount = unbalCount + 1.0;
                    }
                }
            }
        }
    }
    return make_pair(balCount, unbalCount);
}

/***
* Sampling framework using improved method for vertexes
*/
void vImpr(int numSamps) {
    if (numSamps > vIDs.size()) {
        numSamps = vIDs.size();
    }
    cout << "Number of Samples: " << numSamps << "\n";
    double estBalCount = 0.0;
    double estUnbalCount = 0.0;
    unordered_map<int, bool> sampledMap;
    mt19937 mtGen(time(nullptr));
    for (int i  = 0; i < numSamps; i++) {
//        cout << i << "\n";
        int sampNode = rSample(vIDs.size(), sampledMap, mtGen);
        pair<double, double> estVals = vImprLocal(vIDs[sampNode]);
        estBalCount = (((double)i * estBalCount) + (estVals.first * (double)vIDs.size() / 3.0)) / ((double)i + 1.0);
        estUnbalCount = (((double)i * estUnbalCount) + (estVals.second * (double)vIDs.size() / 3.0)) / ((double)i + 1.0);
    }
    cout << "Estimated Balanced Count: " << estBalCount << "\n";
    cout << "Estimated Unbalanced Count: " << estUnbalCount << "\n";
}


pair<double, double> eImprLocal(AdvEdge aEdge) {
    unordered_map<int, Node>::iterator uIt = nodeMap.find(aEdge.u1);
    unordered_map<int, Node>::iterator vIt = nodeMap.find(aEdge.u2);
    unordered_map<int, double>::iterator eIt;
    double balCount = 0.0;
    double unbalCount = 0.0;
    double triBalProb;
    double tAEV = calcAEV(t);
    if (aEdge.aev < tAEV) {
        return make_pair(0.0, 0.0);
    }
    if (uIt->second.neighbours.size() >= vIt->second.neighbours.size()) {
        //vIt is smaller
        int iV = 0;
        for (; iV < vIt->second.neighbours.size(); iV++) {
            if (calcAEV(vIt->second.neighbours[iV].second) < tAEV) {
                break;
            }
            eIt = uIt->second.neighbourMap.find(vIt->second.neighbours[iV].first);
            if (eIt != uIt->second.neighbourMap.end()) {
                triBalProb = calcBalProb(aEdge.prob, vIt->second.neighbours[iV].second, eIt->second);
                if (triBalProb >= t) {
                    balCount = balCount + 1.0;
                } else if ((1.0 - triBalProb) > t) {
                    unbalCount = unbalCount + 1.0;
                }
            }
        }
    } else {
        //uIt is smaller
        int iU = 0;
        for (; iU < uIt->second.neighbours.size(); iU++) {
            if (calcAEV(uIt->second.neighbours[iU].second) < tAEV) {
                break;
            }
            eIt = vIt->second.neighbourMap.find(uIt->second.neighbours[iU].first);
            if (eIt != vIt->second.neighbourMap.end()) {
                triBalProb = calcBalProb(aEdge.prob, uIt->second.neighbours[iU].second, eIt->second);
                if (triBalProb >= t) {
                    balCount = balCount + 1.0;
                } else if ((1.0 - triBalProb) > t) {
                    unbalCount = unbalCount + 1.0;
                }
            }
        }
    }
    return make_pair(balCount, unbalCount);
}

void eImpr(int numSamps) {
    if (numSamps > advEdgeList.size()) {
        numSamps = advEdgeList.size();
    }
    cout << "Number of Samples: " << numSamps << "\n";
    unordered_map<int, bool> sampledMap;
    double estBalCount = 0.0;
    double estUnbalCount = 0.0;
    mt19937 mtGen(time(nullptr));
    for (int i = 0; i < numSamps; i++) {
        int sampEdge = rSample(advEdgeList.size(), mtGen);
        pair<double, double> estVals = eImprLocal(advEdgeList[sampEdge]);
        estBalCount = (((double)i * estBalCount) + (estVals.first * (double)advEdgeList.size() / 3.0)) / ((double)i + 1.0);
        estUnbalCount = (((double)i * estUnbalCount) + (estVals.second * (double)advEdgeList.size() / 3.0)) / ((double)i + 1.0);
    }
    cout << "Estimated Balanced Count: " << estBalCount << "\n";
    cout << "Estimated Unbalanced Count: " << estUnbalCount << "\n"; 
}


pair<double, double> eBaseLocal(BasicEdge bEdge) {
    unordered_map<int, Node>::iterator uIt = nodeMap.find(bEdge.u1);
    unordered_map<int, Node>::iterator vIt = nodeMap.find(bEdge.u2);
    unordered_map<int, double>::iterator eIt;
    double balCount = 0.0;
    double unbalCount = 0.0;
    double triBalProb;
    if (uIt->second.neighbours.size() >= vIt->second.neighbours.size()) {
        //vIt is smaller
        int iV = 0;
        for (; iV < vIt->second.neighbours.size(); iV++) {
            eIt = uIt->second.neighbourMap.find(vIt->second.neighbours[iV].first);
            if (eIt != uIt->second.neighbourMap.end()) {
                triBalProb = calcBalProb(bEdge.prob, vIt->second.neighbours[iV].second, eIt->second);
                if (triBalProb >= t) {
                    balCount = balCount + 1.0;
                } else if ((1.0 - triBalProb) > t) {
                    unbalCount = unbalCount + 1.0;
                }
            }
        }
    } else {
        //uIt is smaller
        int iU = 0;
        for (; iU < uIt->second.neighbours.size(); iU++) {
            eIt = vIt->second.neighbourMap.find(uIt->second.neighbours[iU].first);
            if (eIt != vIt->second.neighbourMap.end()) {
                triBalProb = calcBalProb(bEdge.prob, uIt->second.neighbours[iU].second, eIt->second);
                if (triBalProb >= t) {
                    balCount = balCount + 1.0;
                } else if ((1.0 - triBalProb) > t) {
                    unbalCount = unbalCount + 1.0;
                }
            }
        }
    }
    return make_pair(balCount, unbalCount);
}

void eBase(int numSamps) {
    if (numSamps > basicEdgeList.size()) {
        numSamps = basicEdgeList.size();
    }
    cout << "Number of Samples: " << numSamps << "\n";
    unordered_map<int, bool> sampledMap;
    double estBalCount = 0.0;
    double estUnbalCount = 0.0;
    mt19937 mtGen(time(nullptr));
    for (int i = 0; i < numSamps; i++) {
        int sampEdge = rSample(basicEdgeList.size(), mtGen);
        pair<double, double> estVals = eBaseLocal(basicEdgeList[sampEdge]);
        estBalCount = (((double)i * estBalCount) + (estVals.first * (double)basicEdgeList.size() / 3.0)) / ((double)i + 1.0);
        estUnbalCount = (((double)i * estUnbalCount) + (estVals.second * (double)basicEdgeList.size() / 3.0)) / ((double)i + 1.0);
    }
    cout << "Estimated Balanced Count: " << estBalCount << "\n";
    cout << "Estimated Unbalanced Count: " << estUnbalCount << "\n"; 
}


/***
* For a local node, determine the local count using the baseline method
* Returns balanced and unbalanced count in a pair
*/
pair<double, double> vBaseLocal(int u) {
    unordered_map<int, Node>::iterator uIt = nodeMap.find(u);
    unordered_map<int, Node>::iterator vIt;
    unordered_map<int, double>::iterator eIt;
    int iU = 0;
    int iV, iU2;
    double balCount = 0.0;
    double unbalCount = 0.0;
    double triBalProb;
    for (; iU < uIt->second.neighbours.size(); iU++) {
        vIt = nodeMap.find(uIt->second.neighbours[iU].first);
        if (uIt->second.neighbours.size() >= vIt->second.neighbours.size()) {
            //vIt is smaller
            iV = 0;
            for (; iV < vIt->second.neighbours.size(); iV++) {
                eIt = uIt->second.neighbourMap.find(vIt->second.neighbours[iV].first);
                if (eIt != uIt->second.neighbourMap.end()) {
                    triBalProb = calcBalProb(uIt->second.neighbours[iU].second, vIt->second.neighbours[iV].second, eIt->second);
                    if (triBalProb >= t) {
                        balCount = balCount + 1.0;
                    } else if ((1.0 - triBalProb) > t) {
                        unbalCount = unbalCount + 1.0;
                    }
                }
            }
        } else {
            //uIt is smaller
            iU2 = 0;
            for (; iU2 < uIt->second.neighbours.size(); iU2++) {
                eIt = vIt->second.neighbourMap.find(uIt->second.neighbours[iU2].first);
                if (eIt != vIt->second.neighbourMap.end()) {
                    triBalProb = calcBalProb(uIt->second.neighbours[iU].second, uIt->second.neighbours[iU2].second, eIt->second);
                    if (triBalProb >= t) {
                        balCount = balCount + 1.0;
                    } else if ((1.0 - triBalProb) > t) {
                        unbalCount = unbalCount + 1.0;
                    }
                }
            }
        }
    }
    return make_pair(balCount / 2.0, unbalCount / 2.0);
}

/***
* Baseline vertex sampling method to estimate count
*/
void vBase(int numSamps) {
    if (numSamps > vIDs.size()) {
        numSamps = vIDs.size();
    }
    cout << "Number of Samples: " << numSamps << "\n";
    double estBalCount = 0.0;
    double estUnbalCount = 0.0;
    unordered_map<int, bool> sampledMap;
    mt19937 mtGen(time(nullptr));
    for (int i  = 0; i < numSamps; i++) {
//        cout << i << "\n";
        int sampNode = rSample(vIDs.size(), sampledMap, mtGen);
        pair<double, double> estVals = vBaseLocal(vIDs[sampNode]);
        estBalCount = (((double)i * estBalCount) + (estVals.first * (double)vIDs.size() / 3.0)) / ((double)i + 1.0);
        estUnbalCount = (((double)i * estUnbalCount) + (estVals.second * (double)vIDs.size() / 3.0)) / ((double)i + 1.0);
    }
    cout << "Estimated Balanced Count: " << estBalCount << "\n";
    cout << "Estimated Unbalanced Count: " << estUnbalCount << "\n";
}


/***
* Improved algorithm for counting the number of uncertain balanced and unbalanced triangles
*/
void imprCount() {
    int balCount = 0;
    int unbalCount = 0;
    double triBalProb, aevBound;
    double tAEV = calcAEV(t);
    unordered_map<int, Node>::iterator uIt, vIt;
    unordered_map<int, double>::iterator eIt;
    for (AdvEdge &aEdge : advEdgeList) {
        if (aEdge.aev < tAEV) {
            break;
        }
        uIt = nodeMap.find(aEdge.u1);
        vIt = nodeMap.find(aEdge.u2);
        
        aevBound = calcAEV(calcAEVBound(aEdge.prob));
        if (uIt->second.neighbours.size() < vIt->second.neighbours.size()) {
            int iU = advEdgeOrderGTBS(uIt->second.neighbours, 0, uIt->second.neighbours.size() - 1, aEdge);
            for (; iU < uIt->second.neighbours.size(); iU++) {
                if (calcAEV(uIt->second.neighbours[iU].second) < aevBound) {
                    break;
                }
                eIt = vIt->second.neighbourMap.find(uIt->second.neighbours[iU].first);
                if (eIt != vIt->second.neighbourMap.end()) {
                    if (checkAEOrder(aEdge, (*eIt))) {
                        triBalProb = calcBalProb(aEdge.prob, eIt->second, uIt->second.neighbours[iU].second);
                        if (triBalProb >= t) {
                            balCount++;
                        } else if ((1.0 - triBalProb) > t) {
                            unbalCount++;
                        }
                    }

                }


            }
        } else {
            int iV = advEdgeOrderGTBS(vIt->second.neighbours, 0, vIt->second.neighbours.size() - 1, aEdge);
            for (; iV < vIt->second.neighbours.size(); iV++) {  
                if (calcAEV(vIt->second.neighbours[iV].second) < aevBound) {
                    break;
                }
                eIt = uIt->second.neighbourMap.find(vIt->second.neighbours[iV].first);
                if (eIt != uIt->second.neighbourMap.end()) {
                    if (checkAEOrder(aEdge, (*eIt))) {
                        triBalProb = calcBalProb(aEdge.prob, eIt->second, vIt->second.neighbours[iV].second);
                        if (triBalProb >= t) {
                            balCount++;
                        } else if ((1.0 - triBalProb) > t) {
                            unbalCount++;
                        }
                    }

                }

                

            }
        }
    }
    cout << "Balanced Count: " << balCount << "\n";
    cout << "Unbalanced Count: " << unbalCount << "\n";
}

/***
* Baseline algorithm for counting the number of uncertain balanced and unbalanced triangles
*/
void baselineCount() {
    int balCount = 0;
    int unbalCount = 0;
    unordered_map<int, Node>::iterator uIt, vIt;
    unordered_map<int, double>::iterator eIt;
    double triBalProb;
    int iU, iV, iU2;
    for (int u : vIDs) {
        uIt = nodeMap.find(u);
        iU = nodeOrderGTBS(uIt->second.neighbours, 0, uIt->second.neighbours.size() - 1, uIt->second); //Node Order
        for (; iU < uIt->second.neighbours.size(); iU++) {
            vIt = nodeMap.find(uIt->second.neighbours[iU].first);
            iV = nodeOrderGTBS(vIt->second.neighbours, 0, vIt->second.neighbours.size() - 1, vIt->second);
            iU2 = nodeOrderGTBS(uIt->second.neighbours, iU, uIt->second.neighbours.size() - 1, vIt->second);
            
            if ((uIt->second.neighbours.size() - iU2) >= (vIt->second.neighbours.size() - iV)) {
                //vIt is the smaller set
                for (; iV < vIt->second.neighbours.size(); iV++) {
                    eIt = uIt->second.neighbourMap.find(vIt->second.neighbours[iV].first);
                    if (eIt != uIt->second.neighbourMap.end()) {
                        triBalProb = calcBalProb(uIt->second.neighbours[iU].second, vIt->second.neighbours[iV].second, eIt->second);
                        if (triBalProb >= t) {
                            balCount++;
                        } else if ((1.0 - triBalProb) > t) {
                            unbalCount++;
                        }

                    }
                }
            } else {
                //uIt is the smaller set
                for (; iU2 < uIt->second.neighbours.size(); iU2++) {
                    eIt = vIt->second.neighbourMap.find(uIt->second.neighbours[iU2].first);
                    if (eIt != vIt->second.neighbourMap.end()) {
                        triBalProb = calcBalProb(uIt->second.neighbours[iU].second, uIt->second.neighbours[iU2].second, eIt->second);
                        if (triBalProb >= t) {
                            balCount++;
                        } else if ((1.0 - triBalProb) > t) {
                            unbalCount++;
                        }
                    }
                }
            }
            
        }
    }
    cout << "Balanced Count: " << balCount << "\n";
    cout << "Unbalanced Count: " << unbalCount << "\n";
}


void initialiseImprEdges(char* filename) {
    ifstream nameStream(filename);
    if (!nameStream.is_open()) {
		printf("Edges File DNE\n");
		exit(3);
	}
    char line[40];
    char* tempV;
    int v1, v2;
    double signProb;
    unordered_map<int, Node>::iterator nIt1, nIt2;
    
    while (nameStream.good()) {
        nameStream.getline(line, 40);
        if (nameStream.eof()) {
            nameStream.close();
            break;
        }
        
        tempV = strtok(line, "\t");
        v1 = stoi(tempV);
        tempV = strtok(NULL, "\t");
        v2 = stoi(tempV);
        tempV = strtok(NULL, "\t");
		signProb = stod(tempV);
        
        //Ensures v1 is smaller than v2 for ordering purposes
        if (v2 < v1) {
            int temp = v2;
            v1 = v2;
            v2 = temp;
        }
        
        nIt1 = nodeMap.find(v1);
        if (nIt1 == nodeMap.end()) {
            Node n;
            n.id = v1;
            nodeMap.insert(make_pair(v1, n));
            nIt1 = nodeMap.find(v1);
            vIDs.push_back(v1);
        }
        
        nIt2 = nodeMap.find(v2);
        if (nIt2 == nodeMap.end()) {
            Node n;
            n.id = v2;
            nodeMap.insert(make_pair(v2, n));
            nIt2 = nodeMap.find(v2);
            vIDs.push_back(v2);
        }
        
        nIt1->second.neighbours.push_back(make_pair(v2, signProb));
        nIt1->second.neighbourMap.insert(make_pair(v2, signProb));
        nIt2->second.neighbours.push_back(make_pair(v1, signProb));
        nIt2->second.neighbourMap.insert(make_pair(v1, signProb));
        
        //Adds new edge
        AdvEdge e;
        e.u1 = v1;
        e.u2 = v2;
        e.prob = signProb;
        e.aev = calcAEV(signProb);
        advEdgeList.push_back(e);
    }
    
    AbsoluteEdgeOrderComp aeoc;
    sort(advEdgeList.begin(), advEdgeList.end(), aeoc);
    ProbToAEVComp ptac;
    for (auto &a : nodeMap) {
		sort(a.second.neighbours.begin(), a.second.neighbours.end(), ptac);
	}
    
    cout << "Nodes: " << vIDs.size() << "\n";
    cout << "Edges: " << advEdgeList.size() << "\n";
}

void initialiseBasicEdges(char* filename) {
    ifstream nameStream(filename);
    if (!nameStream.is_open()) {
		printf("Edges File DNE\n");
		exit(3);
	}
    char line[40];
    char* tempV;
    int v1, v2;
    double signProb;
    unordered_map<int, Node>::iterator nIt1, nIt2;
        
    while (nameStream.good()) {
        nameStream.getline(line, 40);
        if (nameStream.eof()) {
            nameStream.close();
            break;
        }
        
        tempV = strtok(line, "\t");
        v1 = stoi(tempV);
        tempV = strtok(NULL, "\t");
        v2 = stoi(tempV);
        tempV = strtok(NULL, "\t");
		signProb = stod(tempV);
        
        //Ensures v1 is smaller than v2 for ordering purposes
        if (v2 < v1) {
            int temp = v2;
            v1 = v2;
            v2 = temp;
        }
        
        nIt1 = nodeMap.find(v1);
        if (nIt1 == nodeMap.end()) {
            Node n;
            n.id = v1;
            nodeMap.insert(make_pair(v1, n));
            nIt1 = nodeMap.find(v1);
            vIDs.push_back(v1);
        }
        
        nIt2 = nodeMap.find(v2);
        if (nIt2 == nodeMap.end()) {
            Node n;
            n.id = v2;
            nodeMap.insert(make_pair(v2, n));
            nIt2 = nodeMap.find(v2);
            vIDs.push_back(v2);
        }
        
        nIt1->second.neighbours.push_back(make_pair(v2, signProb));
        nIt1->second.neighbourMap.insert(make_pair(v2, signProb));
        nIt2->second.neighbours.push_back(make_pair(v1, signProb));
        nIt2->second.neighbourMap.insert(make_pair(v1, signProb));
        
        //Adds new edge
        BasicEdge e;
        e.u1 = v1;
        e.u2 = v2;
        e.prob = signProb;
        basicEdgeList.push_back(e);
    }
    
    NodeOrderComp noc;
    sort(vIDs.begin(), vIDs.end(), noc);
    
    for (auto &a : nodeMap) {
		sort(a.second.neighbours.begin(), a.second.neighbours.end(), noc);
	}
    
    cout << "Nodes: " << vIDs.size() << "\n";
    cout << "Edges: " << basicEdgeList.size() << "\n";
    
}

/***
* argv[0]: .\baltri
* argv[1]: Process Number
*   1: Baseline Counting
*   2: Improved Counting
*   3: vBase
*   4: vImpr
*   5: eBase
*   6: eImpr
* argv[2]: Edge File
* argv[3]: t
* argv[4]: Number of Samples (Process 3-6)
*/
int main(int argc, char *argv[]) {
    int processNum = stoi(argv[1]);
    t = stod(argv[3]);
    auto timeStart = chrono::high_resolution_clock().now();
    if (processNum == 1) {
        cout << "Baseline Counting\n";
        initialiseBasicEdges(argv[2]);
        timeStart = chrono::high_resolution_clock().now();
        baselineCount();
    } else if (processNum == 2) {
        cout << "Improved Counting\n";
        initialiseImprEdges(argv[2]);
        timeStart = chrono::high_resolution_clock().now();
        imprCount();
    } else if (processNum == 3) {
        cout << "vBase Sampling\n";
        initialiseBasicEdges(argv[2]);
        timeStart = chrono::high_resolution_clock().now();
        vBase(stoi(argv[4]));
    } else if (processNum == 4) {
        cout << "vImpr Sampling\n";
        initialiseImprEdges(argv[2]);
        timeStart = chrono::high_resolution_clock().now();
        vImpr(stoi(argv[4]));
    } else if (processNum == 5) {
        cout << "eBase Sampling\n";
        initialiseBasicEdges(argv[2]);
        timeStart = chrono::high_resolution_clock().now();
        eBase(stoi(argv[4]));
    } else if (processNum == 6) {
        cout << "eImpr Sampling\n";
        initialiseImprEdges(argv[2]);
        timeStart = chrono::high_resolution_clock().now();
        eImpr(stoi(argv[4]));
    } 
    
    auto timeEnd = chrono::high_resolution_clock().now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(timeEnd - timeStart).count();
	cout << "Total Runtime: " << duration << "\n";
}