//////////////////////////////////////////////////////////////////////
// analyze.h  Copyright (c) 2014 Dario Ghersi                       //
// Version: 20140307                                                //
//          See User's Guide for more details                       //
//                                                                  //
// This file is part of the MOLBLOCKS suite.                        //
// MOLBLOCKS is free software: you can redistribute it and/or       //
// modify it under the terms of the GNU General Public License as   //
// published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.              //
//                                                                  //
// MOLBLOCKS is distributed in the hope that it will be useful,     //
// but WITHOUT ANY WARRANTY; without even the implied warranty of   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    //
// GNU General Public License for more details.                     //
//                                                                  //
// You should have received a copy of the GNU General Public        //
// License along with MOLBLOCKS.                                    //
// If not, see <http://www.gnu.org/licenses/>.                      //
//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// CONSTANTS                                                        //
//////////////////////////////////////////////////////////////////////

const string USAGE = "analyze -i FRAGMENTS [-o OUT_DISTR -e BACKGROUND\n\
                              -c TANIMOTO]\n";
typedef pair<string, unsigned int> stringInt;
typedef pair<unsigned int, unsigned int> intInt;
typedef pair<unsigned int, double> intDouble;

typedef boost::adjacency_list<boost::listS, boost::vecS,\
  boost::undirectedS> MyGraph;
typedef boost::graph_traits<MyGraph>::vertex_descriptor MyVertex;
vector<unsigned int> visitedVertices; // global variable, to
                                     // store BFS results
const unsigned int minFragEnrich = 3; // minimum frequency to carry
                                      // out the hypergeometric test

const double TANIMOTO = 0.8; // default tanimoto clustering threshold

//////////////////////////////////////////////////////////////////////
// CLASSES                                                          //
//////////////////////////////////////////////////////////////////////

class MyVisitor : public boost::default_bfs_visitor
{
  // visitor class for the breadth-first search with Boost Library

public:
  void discover_vertex(MyVertex v, const MyGraph& g) const
  {
    visitedVertices.push_back(v);
    return;
  }
};

//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  char *backgroundFileName;
  char *fragmentFileName;
  char *outFileName;
  bool clusterFragments;
  bool doEnrichment;
  double tanimotoCutoff;

  Parameters(char **, int);
};

//////////////////////////////////////////////////////////////////////

class Fragments {

  // store whether all molecules have names
  bool allMolHaveNames;

  // fragment -> molecules that contain it (indices only,
  // the molecules don't have to have names)
  vector<vector<unsigned int> > frag2IndexMain;
  vector<vector<unsigned int> > frag2IndexBack;

  // fragment -> molecules that contain it (strings, the molecules
  // must have names)
  vector<vector<string> > frag2Mol;

  // strings with the names of the fragments
  vector<string> fragmStr;

  // string and frequency of the main set and background sets
  vector<intInt> fragCounts;
  vector<intInt> backgroundCounts;

  // string of the representative and frequency for the main and
  // background sets
  vector<intInt> fragClustersCounts;
  vector<intInt> backgroundClustersCounts;

  // Babel objects with the main set and the background set
  vector<OBMol> fragments;
  vector<OBMol> background;

  // p-values and FDR adjusted p-values
  vector<double> pvalues;
  vector<double> adjustedP;

  // fingerprints and cluster membership for the clusters
  vector<vector<unsigned int> > fp;
  vector<vector<unsigned int> > clusters;

  // the representative fragments for each cluster
  vector<unsigned int> representatives;

 public:
  Fragments(Parameters);
  void calculateFP(Parameters);
  void cluster(Parameters);
  void countFragments(vector<intInt> &, vector<intInt> &,
		      vector<vector <unsigned int> > &);
  void enrichmentAnalysis(Parameters);
  vector<unsigned int> fdrCorrection();
  void findRepresentatives();
  void getFrequency(bool);
  void printDistribution(Parameters);
  void printEnrichResults(Parameters, vector<unsigned int>,
			  vector<intInt> &, vector<intInt> &,
			  vector<unsigned int> &);
  void sort(bool, bool);
  void storeFragments(char *, vector<intInt> &,
		      vector<vector<unsigned int> > &, bool);
};

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

bool comparator1 (const intInt &, const intInt &);
bool comparator2 (const intDouble &, const intDouble &);
void checkCommandLineArgs(char **, int);
