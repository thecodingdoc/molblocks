//////////////////////////////////////////////////////////////////////
// analyze.C  Copyright (c) 2014 Dario Ghersi                       //
// Version: 20140307                                                //
// Goal:    analyze fragments of small molecules, building          //
//          statistics, clustering fragments and performing         //
//          enrichment analysis                                     //
// Usage:   analyze -i FRAGMENTS [-o OUT_DISTR -e BACKGROUND        //
//                                 -c TANIMOTO]                     //
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

#include <algorithm>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/math/distributions/hypergeometric.hpp>
#include <iostream>
#include <limits>
#include <openbabel/fingerprint.h>
#include <openbabel/obconversion.h>
#include <openbabel/mol.h>
#include <vector>

using namespace boost;
using namespace OpenBabel;
using namespace std;
#include "analyze.h"
#include "utilities.h"

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

bool comparator1(const intInt &pair1, const intInt &pair2)
{
  // comparison function

  return pair1.second > pair2.second;
}

//////////////////////////////////////////////////////////////////////

bool comparator2(const intDouble &pair1, const intDouble &pair2)
{
  // comparison function

  return pair1.second < pair2.second;
}

//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **argv, int argc)
{

  bool err = false;
  if (!cmdOptionExists(argv, argv+argc, "-i")) {
    cerr << "Fragments file missing\n" << USAGE;
    err = true;
  }

  if (err)
    exit(1);
}

//////////////////////////////////////////////////////////////////////

Fragments::Fragments(Parameters p)
{
  // store the fragments (main set and optionally background)
  unsigned int oldPercentage = 0, percentage = 0;

  // store the background (optional)
  bool storeMol = true;
  if (p.doEnrichment) {
    cout << "\n\nStoring the background fragments\n\n";
    storeFragments(p.backgroundFileName, backgroundCounts,
		   frag2IndexBack, false);
    storeMol = true;

    // resize the frag2IndexMain vector
    frag2IndexMain.resize(frag2IndexBack.size());
  }

  // store the main set
  cout << "\nStoring the main fragment set\n\n";
  storeFragments(p.fragmentFileName, fragCounts,
		 frag2IndexMain, storeMol);

  OBConversion conv;
  conv.SetInAndOutFormats("SMI", "CAN");

  // if clustering is chosen, convert the strings to OBMol objects
  if (p.clusterFragments && p.doEnrichment) {
    cout << "\n\nProcessing the fragments\n\n";

    unsigned int i = 0;
    unsigned int backgroundSize = backgroundCounts.size();
    for (i = 0; i < backgroundSize; i++) {
      OBMol mol;
      conv.ReadString(&mol, fragmStr[backgroundCounts[i].first]);
      background.push_back(mol);

      // call the progress bar function every 100 items
      if ((i % 100) == 0 || backgroundSize < 100) {
	percentage = 100.0 * i / backgroundSize;
	if (percentage > oldPercentage) {
	  oldPercentage = percentage;
	  printProgBar(percentage);
	}
      }
    }
    printProgBar(100.0 * i / backgroundSize);
    cout << endl;
  }
  else if (p.clusterFragments) {
    cout << "\n\nProcessing the fragments\n\n";
    unsigned int i = 0;
    unsigned int fragCountsSize = fragCounts.size();

    for (i = 0; i < fragCountsSize; i++) {
      OBMol mol;
      conv.ReadString(&mol, fragmStr[fragCounts[i].first]);
      fragments.push_back(mol);

      // call the progress bar function every 100 items
      if ((i % 100) == 0 || fragCountsSize < 100) {
	percentage = 100.0 * i / fragCountsSize;
	if (percentage > oldPercentage) {
	  oldPercentage = percentage;
	  printProgBar(percentage);
	}
      }
    }
    printProgBar(100.0 * i / fragCountsSize);
    cout << endl;
  }
}

//////////////////////////////////////////////////////////////////////

void Fragments::calculateFP(Parameters p)
{
  // calculate fingerprints of the fragments
  unsigned int oldPercentage = 0, percentage = 0;

  vector<OBMol> *allFrags;
  // if we will do enrichment analysis, then work on the background
  if (p.doEnrichment) {
    allFrags = &background;
  }
  // otherwise just work on the main fragment set
  else {
    allFrags = &fragments;
  }

  // calculate the fingerprints for each fragment
  cout << "\n\nFingerprinting the fragments\n\n";
  unsigned int i = 0;
  unsigned int allFragsSize = (*allFrags).size();
  for (i = 0; i < allFragsSize; i++) {
    vector<unsigned int> currentFp;
    OBFingerprint *pFP = OBFingerprint::FindFingerprint("FP2");
    pFP->GetFingerprint(&(*allFrags)[i], currentFp);
    fp.push_back(currentFp);

    if ((i % 50) == 0 || allFragsSize < 50) {
      percentage = 100.0 * i / allFragsSize;
      if (percentage > oldPercentage) {
	oldPercentage = percentage;
	printProgBar(percentage);
      }
    }
  }
  printProgBar(100.0 * i / allFragsSize);
  cout << endl;
}

//////////////////////////////////////////////////////////////////////

void Fragments::cluster(Parameters p)
{
  // build a graph where there is an edge between two fragments
  // if their Tanimoto coefficient is > tanimotoCutoff.
  // Then, extract the connected components

  unsigned int oldPercentage = 0, percentage = 0;

  // calculate the fingerprint (type FP2) for each fragment
  calculateFP(p);

  // build the fragment graph
  MyGraph g(fp.size());
  OBFingerprint *pFP = OBFingerprint::FindFingerprint("FP2");
  double tanimoto;
  cout << "\n\nCalculating the fragment pairwise similarity\n\n";
  unsigned int i = 0;
  unsigned long fpSize = fp.size();
  unsigned long numPairs = fpSize * (fpSize - 1) / 2, count = 0;
  printProgBar(0);
  for (i = 0; i < fpSize - 1; i++) {
    for (unsigned int j = i + 1; j < fpSize; j++) {

      // calculate the tanimoto coefficient
      tanimoto = pFP->Tanimoto(fp[i], fp[j]);

      if (tanimoto > p.tanimotoCutoff) { // add an edge
	add_edge(i, j, g);
      }

      if ((count % 100) == 0 || numPairs < 100) {
	percentage = 100.0 * count / numPairs;
	if (percentage > oldPercentage) {
	  oldPercentage = percentage;
	  printProgBar(percentage);
	}
      }
      count++;
    }
  }
  printProgBar(100);
  cout << "\n\n";

  // get the vertices in the graph
  typedef adjacency_list<vecS, vecS, undirectedS> Graph;
  typedef property_map<Graph, vertex_index_t>::type IndexMap;
  typedef graph_traits<Graph>::vertex_iterator vertex_iter;
  IndexMap index = get(vertex_index, g);
  std::pair<vertex_iter, vertex_iter> vp;

  // find the connected component
  MyVisitor vis;
  cout << "Calculating the connected components\n\n";
  for (vp = vertices(g); vp.first != vp.second; ++vp.first) {

    // check whether the current vertex is already in a component
    bool found = false; // if we find the element, set to true
    for (unsigned int i = 0; i < clusters.size(); i++) {
      if (found)
	break;
      for (unsigned int j = 0; j < clusters[i].size(); j++) {
	if (index[*vp.first] == clusters[i][j]) {
	  found = true;
	  break;
	}
      }
    }

    // start the search from the element
    if (!found) {
      visitedVertices.clear();
      boost::breadth_first_search(g, vertex(index[*vp.first], g),
				  boost::visitor(vis));

      clusters.push_back(visitedVertices);
    }
  }

  // find the representative fragments for each connected component
  cout << "Finding the representative fragments\n\n";
  findRepresentatives();
}

//////////////////////////////////////////////////////////////////////

void Fragments::countFragments(vector<intInt> &clCounts,
			       vector<intInt> &counts,
			       vector<vector <unsigned int> > &refFrag2Index) {
  // compute the frequency of each cluster

  unsigned int maxIndex = 0; // keep track of the highest index

  // initialize the background counts
  for (unsigned int i = 0; i < clusters.size(); i++) {
    intInt foo;
    foo = make_pair(i, 0);
    clCounts.push_back(foo);
  }

  // fill in the background counts
  for (unsigned int i = 0; i < clusters.size(); i++) {
    vector<unsigned int> alreadyUsed;
    for (unsigned int j = 0; j < clusters[i].size(); j++) {
      for (unsigned int k = 0; k < refFrag2Index[clusters[i][j]].size();
	   k++) {
	unsigned int element = refFrag2Index[clusters[i][j]][k];
	bool found = find(alreadyUsed.begin(), alreadyUsed.end(),
			  element) != alreadyUsed.end();

	// increment the count if the molecule has not contributed
	// to the cluster already
	if (!found) {
	  alreadyUsed.push_back(element);
	  clCounts[i].second++;
	}
      }

      if (clusters[i][j] > maxIndex) {
	maxIndex = clusters[i][j];
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Fragments::enrichmentAnalysis(Parameters p)
{
  // perform enrichment analysis on the fragments using the
  // hypergeometric distribution

  // determine whether we are working with clusters or not
  vector<intInt> &refFragCounts = (p.clusterFragments) ?
    fragClustersCounts : fragCounts;
  vector<intInt> &refBackCounts = (p.clusterFragments) ?
    backgroundClustersCounts : backgroundCounts;

  // calculate the total number of fragments in the main set
  unsigned int t = 0;
  for (unsigned int i = 0; i < refFragCounts.size(); i++) {
    t += refFragCounts[i].second;
  }

  // calculate the total number of fragments in the background
  unsigned int totalFrag = 0;
  for (unsigned int i = 0; i < refBackCounts.size(); i++) {
    totalFrag += refBackCounts[i].second;
  }

  // process each fragment in turn
  unsigned int k, n1;
  double pval;
  vector<unsigned int> indexesBack; // indexes of the background

  for (unsigned int i = 0; i < refFragCounts.size(); i++) {

    // get the number of "white balls drawn from the urn"
    k = refFragCounts[i].second;

    // get the total number of "white balls"
    n1 = 0;

    for (unsigned int j = 0; j < refBackCounts.size(); j++) {
      if (refFragCounts[i].first == refBackCounts[j].first) {
	n1 = refBackCounts[j].second;
	indexesBack.push_back(j);
	break;
      }
    }

    // sanity check
    if (t >= totalFrag) {
      cerr << "\nThe main set should be a proper subset of the";
      cerr << " background!\nAborting...\n";
      exit(1);
    }

    // calculate the p-value
    math::hypergeometric_distribution<double> hyper(n1, t, totalFrag);
    pval = 1.0 - cdf(hyper, k - 1);
    pvalues.push_back(pval);
  }

  // apply the FDR correction
  vector<unsigned int> sortedOrder = fdrCorrection();

  printEnrichResults(p, sortedOrder, refFragCounts, refBackCounts,
		     indexesBack);
}

//////////////////////////////////////////////////////////////////////

vector<unsigned int> Fragments::fdrCorrection()
{
  // adjust the p-values by applying the Benjamini-Hochberg correction

  // combine the p-values wih their index for sorting
  vector<intDouble> pPairs;
  intDouble foo;
  for (unsigned int i = 0; i < pvalues.size(); i++) {
    foo = make_pair(i, pvalues[i]);
    pPairs.push_back(foo);
  }

  // sort the p-values
  std::sort(pPairs.begin(), pPairs.end(), comparator2);

  // apply the Benjamini-Hochberg procedure to correct the p-values
  // for multiple hypothesis testing
  double currMin, value;
  vector<double> minValues;
  unsigned int m = pPairs.size();
  for (unsigned int i = 0; i < m; i++) {
    currMin = pPairs[i].second * m / (i + 1);
    for (unsigned int j = i + 1; j < m; j++) {
      value = pPairs[j].second * m / (j + 1);
      if (value < currMin) {
	currMin = value;
      }
    }
    minValues.push_back(min(currMin, 1.0));
  }

  // put the adjusted values in the original p-value order
  adjustedP.resize(pPairs.size());
  vector <unsigned int> sortedOrder(pPairs.size());
  for (unsigned int i = 0; i < pPairs.size(); i++) {
    adjustedP[pPairs[i].first] = minValues[i];
    sortedOrder[i] = pPairs[i].first;
  }

  return sortedOrder;
}

//////////////////////////////////////////////////////////////////////

void Fragments::findRepresentatives()
{
  // for each cluster, find the fragment with the highest pairwise
  // average similarity to all the other fragments in the same
  // cluster

  OBFingerprint *pFP = OBFingerprint::FindFingerprint("FP2");
  double tanimoto;
  unsigned int oldPercentage = 0, percentage = 0;

  // get the total number of unique fragments
  unsigned int numFrag =
    (fragments.size() > background.size()) ? fragments.size():
    background.size();

  // process each cluster
  unsigned long clustersSize = clusters.size();
  unsigned int i = 0, count = 0;
  printProgBar(0);
  for (i = 0; i < clustersSize; i++) {

    // set the previous total similarity to a negative value
    double previousSim = -1.0;
    unsigned int representative;

    for (unsigned int j = 0; j < clusters[i].size(); j++) {
      double totalSim = 0.0;
      for (unsigned int k = 0; k < clusters[i].size(); k++) {
	if (j != k) {
	  tanimoto = pFP->Tanimoto(fp[clusters[i][j]],
				   fp[clusters[i][k]]);
	  totalSim += tanimoto;
	}
      }
      if (totalSim > previousSim) {
	previousSim = totalSim;
	representative = clusters[i][j];
      }

      // progress bar
      if ((count % 100) == 0 || numFrag < 100) {
	percentage = 100.0 * count / numFrag;
	if (percentage > oldPercentage) {
	  oldPercentage = percentage;
	  printProgBar(percentage);
	}
      }
      count++;
    }

    representatives.push_back(representative);
  }
  printProgBar(100.0 * i / clustersSize);
  cout << endl;
}

//////////////////////////////////////////////////////////////////////

void Fragments::getFrequency(bool doEnrichment)
{
  // compute the frequency of occurrence of each cluster in the
  // main set and (optionally) in the background set

  if (doEnrichment) {
    // calculate the background frequency
    countFragments(backgroundClustersCounts, backgroundCounts,
		   frag2IndexBack);

    // initialize the counts
    vector<unsigned int> sum(clusters.size(), 0);

    // calculate the main set frequency
    for (unsigned int i = 0; i < clusters.size(); i++) {
      vector<unsigned int> alreadyUsed;
      for (unsigned int j = 0; j < clusters[i].size(); j++) {
	for (unsigned int k = 0;
	     k < frag2IndexMain[clusters[i][j]].size(); k++) {
	  unsigned int element = frag2IndexMain[clusters[i][j]][k];
	  bool found = find(alreadyUsed.begin(), alreadyUsed.end(),
			    element) != alreadyUsed.end();
	  if (!found) {
	    sum[i]++;
	    alreadyUsed.push_back(element);
	  }
	}
      }
    }

    // build the main set frequency vector
    for (unsigned int i = 0; i < clusters.size(); i++) {
      if (sum[i] > 0) {
	intInt foo;
	foo = make_pair(backgroundCounts[i].first,
			sum[i]);
	fragClustersCounts.push_back(foo);
      }
    }
  }

  // no enrichment
  else {
    // calculate the main set frequency
    countFragments(fragClustersCounts, fragCounts,
		   frag2IndexMain);
  }
}

//////////////////////////////////////////////////////////////////////

void Fragments::printDistribution(Parameters p)
{
  // print the distribution of fragments, sorted by frequency

  // open the output file
  fstream outFile;
  outFile.open(p.outFileName, fstream::out);

  // print each fragment, with the frequency of occurrence
  if (p.clusterFragments) { // clustering

    for (unsigned int i = 0; i < fragClustersCounts.size(); i++) {

      unsigned int index = fragClustersCounts[i].first;
      string fragString = index < clusters.size() ?
	fragmStr[representatives[index]]: fragmStr[index];
      outFile << fragClustersCounts[i].second << "\t" <<
	fragString;

      string space = "";
      if (allMolHaveNames) {

	outFile << "\t";

	// print for each repres. all the molecules that contain it
	vector<string> alreadyPrinted;
        for (unsigned int j = 0; j < clusters[index].size(); j++) {
	  for (unsigned k = 0;
	       k < frag2Mol[clusters[index][j]].size(); k++) {
	    bool toPrint = true;
	    for (vector<string>::iterator frag = alreadyPrinted.begin();
		 frag != alreadyPrinted.end(); frag++) {
	      if (*frag == frag2Mol[clusters[index][j]][k]) {
		toPrint = false;
	      }
	    }
	    if (toPrint) {
	      outFile << space << frag2Mol[clusters[index][j]][k];
	      space = " ";
	      alreadyPrinted.push_back(frag2Mol[clusters[index][j]][k]);
	    }
	  }
	}
	  outFile << endl;
      }
      else {
	outFile << endl;
      }
    }
  }
  else { // no clustering

    for (unsigned int i = 0; i < fragCounts.size(); i++) {
      outFile << fragCounts[i].second << "\t" <<
	fragmStr[fragCounts[i].first];

      // print for each fragment all the molecules that contain it
      string space = "";
      if (allMolHaveNames) {

	outFile << "\t";

	for (unsigned int j = 0;
	     j < frag2Mol[fragCounts[i].first].size(); j++) {
	  outFile << space << frag2Mol[fragCounts[i].first][j];
	  space = " ";
	}
	outFile << endl;
      }
      else {
	outFile << endl;
      }
    }
  }

  // close the files
  outFile.close();
}

//////////////////////////////////////////////////////////////////////

void Fragments::printEnrichResults(Parameters p,
				   vector<unsigned int> sortedOrder,
				   vector<intInt> &refFragCounts,
				   vector<intInt> &refBackCounts,
				   vector<unsigned int> &indexesBack)
{
  // print the enrichment results
  fstream outFile;
  outFile.open(p.outFileName, fstream::out);
  outFile << std::scientific;

  // print out the header
  outFile <<
    "p-value\tFDR\tFrequency\tBackground\tFragment\tMolecules\n";
  for (unsigned int i = 0; i < refFragCounts.size(); i++) {
    unsigned int index = sortedOrder[i];

    // get the representative fragment, if clustering
    string representative;
    if (p.clusterFragments &&
	refFragCounts[index].first < clusters.size()) {
	representative =
	  fragmStr[representatives[refFragCounts[index].first]];
    }
    else {
      representative = fragmStr[refFragCounts[index].first];
    }

    // round the p-value to epsilon if too small
    double pvalue, adjustedPVal;
    if (pvalues[index] < DBL_EPSILON) {
      pvalue = DBL_EPSILON;
      adjustedPVal = DBL_EPSILON;
    }
    else {
      pvalue = pvalues[index];
      adjustedPVal = adjustedP[index];
    }

    // print the results
    outFile << pvalue << "\t" << adjustedPVal << "\t"\
	    << refFragCounts[index].second << "\t" <<\
	    refBackCounts[indexesBack[index]].second << "\t" <<\
      representative << "\t";

    // print the molecules that contain the fragment(s)
    string space = "";
    vector<string> alreadyPrinted;
      
    if (p.clusterFragments) {
      for (unsigned int j = 0;
	   j < clusters[refFragCounts[index].first].size(); j++) {
	unsigned int currFrag = clusters[refFragCounts[index].first][j];
	if (currFrag < frag2Mol.size()) {
	  for (unsigned int k = 0; k < frag2Mol[currFrag].size(); k++) {
	    bool toPrint = true;
	    for (vector<string>::iterator frag = alreadyPrinted.begin();
		 frag != alreadyPrinted.end(); frag++) {
	      if (*frag == frag2Mol[currFrag][k]) {
		toPrint = false;
	      }
	    }
	    if (toPrint) {
	      outFile << space <<
		frag2Mol[clusters[refFragCounts[index].first][j]][k];
	      space = " ";
	      alreadyPrinted.push_back(frag2Mol[currFrag][k]);
	    }
	  }
	}
      }
      outFile << endl;
    }
    else { // no clustering
      string space = "";
      vector<string> alreadyPrinted;
      for (unsigned int k = 0;
           k < frag2Mol[refFragCounts[index].first].size(); k++) {
        bool toPrint = true;
        for (vector<string>::iterator frag = alreadyPrinted.begin();
             frag != alreadyPrinted.end(); frag++) {
          if (*frag == frag2Mol[refFragCounts[index].first][k]) {
            toPrint = false;
          }
        }
        if (toPrint) {
          outFile << space << frag2Mol[refFragCounts[index].first][k];
          space = " ";
          alreadyPrinted.push_back(frag2Mol[refFragCounts[index].first][k]);
        }
      }
      outFile << endl;
    }
  }
  outFile.close();
}

//////////////////////////////////////////////////////////////////////

void Fragments::sort(bool clusterFragments, bool doEnrichment)
{
  // sort the fragments by frequency of occurrency

  // if the fragments have been clustered, work at the cluster level
  if (clusterFragments) {
    getFrequency(doEnrichment);

    // sort the clusters by frequency
    std::sort(fragClustersCounts.begin(), fragClustersCounts.end(),
	      comparator1);
  }
  // no clustering, just sort the fragments
  else {
    std::sort(fragCounts.begin(), fragCounts.end(), comparator1);
  }
}

//////////////////////////////////////////////////////////////////////

void Fragments::storeFragments(char *fragmentFileName,
			       vector<intInt> &frags,
			       vector<vector<unsigned int> > &refFrag2Index,
			       bool storeMol)
{
  // create pairs with the fragment and its frequency

  // initialize the flag to check whether all molecules have names
  allMolHaveNames = true;

  // open the fragments file and count how many lines are there
  fstream fragFile;
  unsigned long numFrag = 0, count = 0;
  unsigned int oldPercentage = 0, percentage = 0;
  string line;

  fragFile.open(fragmentFileName, fstream::in);

  // complain if the file doesn't exist
  if (! fragFile.good()) {
    cerr << "Can't open " << fragmentFileName << endl;
    exit(1);
  }

  while (getline(fragFile, line)) {
    numFrag++;
  }

  // rewind the file
  fragFile.clear();
  fragFile.seekg(0, ios::beg);

  // declare the variables
  string allFrags, frag, molName;
  vector<unsigned int> fragCount;

  // initialize the fragment counts
  frag2Mol.resize(fragmStr.size());
  for (unsigned int i = 0; i < fragmStr.size(); i++) {
    fragCount.push_back(0);
    vector<string> empty;
    frag2Mol[i] = empty;
  }

  for (unsigned int index = 0; index < numFrag; index++) {
    getline(fragFile, line);

    // get the fragments
    istringstream iss(line);
    iss >> allFrags;

    iss >> molName; // molecule's name
    // check if the name is empty
    if (molName.empty()) {
      allMolHaveNames = false;
    }

    istringstream iss2(allFrags);

    while (getline(iss2, frag, '.')) {
      bool found = false;

      // check if the fragment has already shown up
      for (unsigned int i = 0; i < fragmStr.size(); i++) {
	if (fragmStr[i] == frag) {
	  fragCount[i]++;
	  found = true;
	  refFrag2Index[i].push_back(index);
	  if (storeMol) {
	    frag2Mol[i].push_back(molName);
	  }
	  break;
	}
      }

      if (!found) {
	fragmStr.push_back(frag);
	fragCount.push_back(1);
	vector<string> molNames(1);
	vector<unsigned int> indices(1);
	indices[0] = index;
	molNames[0] = molName;
        refFrag2Index.push_back(indices);
	if (storeMol) {
	  frag2Mol.push_back(molNames);
	}
      }
    }

    count++;

    // call the progress bar function every 50 items
    if ((count % 50) == 0 || numFrag < 50) {
      percentage = 100.0 * count / numFrag;
      if (percentage > oldPercentage) {
	oldPercentage = percentage;
	printProgBar(percentage);
      }
    }
  }
  printProgBar(100.0 * count / numFrag);
  cout << endl;
  fragFile.close();

  // create pairs with the fragment and its frequency
  intInt foo;
  unsigned int fragSize = fragmStr.size();
  for (unsigned int i = 0; i < fragSize; i++) {
    if (fragCount[i] > 0) {
      foo = make_pair(i, fragCount[i]);
      frags.push_back(foo);
    }
  }
}

//////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  backgroundFileName = getCmdOption(argv, argv + argc, "-e");
  fragmentFileName = getCmdOption(argv, argv + argc, "-i");
  outFileName = getCmdOption(argv, argv + argc, "-o");

  // boolean variables
  clusterFragments = cmdOptionExists(argv, argv + argc, "-c");
  doEnrichment = cmdOptionExists(argv, argv + argc, "-e");

  // get the Tanimoto threshold
  if (clusterFragments) {
    stringstream temp(getCmdOption(argv, argv + argc, "-c"));
    temp >> tanimotoCutoff;
    if (tanimotoCutoff < 0.01) {
      cout << "Using default Tanimoto of 0.8\n";
      tanimotoCutoff = TANIMOTO;
    }
  }
}

//////////////////////////////////////////////////////////////////////
// MAIN PROGRAM                                                     //
//////////////////////////////////////////////////////////////////////

int main(int argc, char **argv)
{

  // check the command-line arguments
  checkCommandLineArgs(argv, argc);

  // get the parameters
  Parameters p(argv, argc);

  // store the fragments
  Fragments fragments(p);

  // cluster the fragments (optional)
  if (p.clusterFragments) {
    fragments.cluster(p);
  }

  // sort the fragments by frequency
  fragments.sort(p.clusterFragments, p.doEnrichment);

  // calculate enrichment
  if (p.doEnrichment) {
    fragments.enrichmentAnalysis(p);
  }
  else { // or just print the fragment distribution
    fragments.printDistribution(p);
  }

  cout << endl;

  return 0;
}
