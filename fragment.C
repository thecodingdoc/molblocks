//////////////////////////////////////////////////////////////////////
// fragment.C  Copyright (c) 2014 Dario Ghersi                      //
// Version: 20140307                                                //
// Goal:    fragment small molecules with-user defined cleavage     //
//          rules (e.g., RECAP rules ), encoded as SMARTS patterns  //
// Usage:   fragmenter -r RULES -i SMALL_MOLECULES -o OUTFILE       //
//                     -n MIN_FRAGMENT_SIZE                         //
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
#include <iostream>
#include <iterator>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/parsmart.h>
#include <sstream>
#include <vector>

using namespace OpenBabel;
using namespace std;
#include "fragment.h"
#include "utilities.h"

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **argv, int argc)
{

  bool err = false;
  if (!cmdOptionExists(argv, argv+argc, "-i")) {
    cerr << "Input file missing\n" << USAGE;
    err = true;
  }
  else if (!cmdOptionExists(argv, argv+argc, "-r")) {
    cerr << "Using default RECAP rules\n";
  }
  else if (!cmdOptionExists(argv, argv+argc, "-o")) {
    cerr << "Output file missing\n" << USAGE;
    err = true;
  }
  else if (!cmdOptionExists(argv, argv+argc, "-n")) {
    cerr << "Please specify the minimum fragment size\n" << USAGE;
    err = true;
  }

  if (err)
    exit(1);
}

//////////////////////////////////////////////////////////////////////

void fragAllMols(Parameters p, vector<OBSmartsPattern> rules)
{
  // main function -- fragment one molecule at a time, and print the
  // results to file
  // N.B. the functions expects the small molecules to be represented
  // as SMILES string, one molecule per line

  fstream smiFile, outFile;
  string line, smiles, name = "";
  OBConversion conv;
  unsigned long totNumMol = 0, count = 0;
  unsigned long oldPercentage = 0, percentage = 0;

  // set the input to SMILES, and the output to canonical SMILES
  conv.SetInAndOutFormats("SMI", "CAN");

  // open the small molecules file
  smiFile.open(p.smiFileName, fstream::in);

  // complain if the file doesn't exist
  if (! smiFile.good()) {
    cerr << "Can't open " << p.smiFileName << endl;
    exit(1);
  }

  // count how many molecules (lines) are in the file
  while (getline(smiFile, line)) {
    totNumMol++;
  }

  // rewind the small molecules file
  smiFile.clear();
  smiFile.seekg(0, ios::beg);

  // open the output file
  outFile.open(p.outFileName, fstream::out);
  cout << endl;

  // process each small molecule
  while (getline(smiFile, line)) {
    // reset the name of the molecule to the empty string
    name = "";

    // extract the SMILES string
    istringstream iss(line);
    iss >> smiles;

    // extract the name (if present)
    iss >> name;

    // build the molecule
    Molecule mol(conv, smiles, name);

    // fragment the molecule
    mol.fragment(rules, conv, p);

    // print the results
    mol.printResults(outFile, conv, p.minFrag);

    count++;

    // call the progress bar function every 10 molecules
    if ((count % 10) == 0 || totNumMol < 10) {
      percentage = 100.0 * count / totNumMol;
      if (percentage > oldPercentage) {
	oldPercentage = percentage;
	printProgBar(percentage);
      }
    }
  }
  printProgBar(100.0 * count / totNumMol);
  cout << endl;

  // close the open streams
  smiFile.close();
  outFile.close();
}

//////////////////////////////////////////////////////////////////////

Molecule::Molecule(OBConversion conv, string smiles, string nn)
{
  // constructor for the Molecule class

  conv.ReadString(&mol, smiles);

  // add the name, if not empty
  if (!nn.empty()) {
    name = nn;
  }
  else {
    name = "";
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::bronKerbosch(vector<int> R, vector<int> P,
			    vector<int> X)
{
  // base of the recursion
  if (P.size() == 0 && X.size() == 0) {
    independentSets.push_back(R);
    return;
  }

  // recursive part
  vector<int> elements = P;
  for (vector<int>::iterator it = elements.begin(); it != elements.end(); it++) {

    // add the element to a copy of R
    vector<int> newR = R;
    newR.push_back(*it);

    // find the neighbors of the element
    vector<int> neigh;
    for (unsigned int j = 0; j < matchPairs.size(); j++) {
      if (depMat[*it][j]) {
	neigh.push_back(j);
      }
    }

    // build the newP vector
    vector<int> newP;
    for (unsigned int j = 0; j < P.size(); j++) {
      for (unsigned int k = 0; k < neigh.size(); k++) {
	if (P[j] == neigh[k]) {
	  newP.push_back(neigh[k]);
	  break;
	}
      }
    }

    // build the newX vector
    vector<int> newX;
    for (unsigned int j = 0; j < X.size(); j++) {
      for (unsigned int k = 0; k < neigh.size(); k++) {
	if (X[j] == neigh[k]) {
	  newX.push_back(neigh[k]);
	  break;	
	}
      }
    }

    // recursive call
    bronKerbosch(newR, newP, newX);

    // add the element to R
    X.push_back(*it);

    // erase the element from P
    unsigned int toRemove = 0;
    for (unsigned int i = 0; i < P.size(); i++) {
      if (P[i] == *it) {
	toRemove = i;
      }
    }
    P.erase(P.begin() + toRemove);
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::calculateDepMat(int minFrag)
{
  // calculate the dependence matrix between cleavable bonds
  // A[i, j] = 1 if bond i is independent from bond j
  // A[i, j] = 0 otherwise

  unsigned int n = matchPairs.size();

  // allocate memory for the dependence matrix
  depMat = (bool **) malloc(sizeof(bool *) * n);
  depMat[0] = NULL;
  for (unsigned int i = 0; i < n; i++) {
    depMat[i] = (bool *) malloc(sizeof(bool) * n);
  }

  // process each bonds
  bool isCleaved;
  int atom1, atom2;
  OBMol temp; // intermediate object to store the molecule after the
              // first cut

  // set the diagonal to "false"
  for (unsigned int i = 0; i < n; i++) {
    depMat[i][i] = false;
  }

  for (unsigned int i = 0; i < n; i++) {

    mol = original; // restore the original molecule

    // get the bond atoms
    atom1 = matchPairs[i][0];
    atom2 = matchPairs[i][1];

    // cleave the bond
    isCleaved = cleaveBond(atom1, atom2, minFrag);

    temp = mol;

    for (unsigned int j = i + 1; j < n; j++) {
      mol = temp; // restore the molecule after the first cut

      // get the bond atoms
      atom1 = matchPairs[j][0];
      atom2 = matchPairs[j][1];

      // attempt to cleave
      isCleaved = cleaveBond(atom1, atom2, minFrag);

      // fill in the dependence matrix
      if (isCleaved) {
	depMat[i][j] = depMat[j][i] = true;
      }
      else {
	depMat[i][j] = depMat[j][i] = false;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

bool Molecule::cleaveBond(int atom1, int atom2, int minFrag)
{
  // cleave the bond between two atoms

  int fragSize1 = fragmentSize(atom1, atom2);
  int fragSize2 = fragmentSize(atom2, atom1);

  // make sure that both fragments are not too small
  if (fragSize1 >= minFrag && fragSize2 >= minFrag) {

    // get the corresponding bond
    OBBond *bond = mol.GetBond(atom1, atom2);
    if (bond) {
      mol.DeleteBond(bond);
      return true;
    }
    else {
      return false; // bond not cleaved
    }
  }

  return false;
}

//////////////////////////////////////////////////////////////////////

void Molecule::cut(int minFrag)
{
  // cut the bonds between the atoms that match the patterns

  for (unsigned int i = 0; i < matchPairs.size(); i++) {
    cleaveBond(matchPairs[i][0], matchPairs[i][1], minFrag);
  } 
}

//////////////////////////////////////////////////////////////////////

void Molecule::cut(int minFrag, int singleton)
{
  // cut the singleton bond

  cleaveBond(matchPairs[singleton][0], matchPairs[singleton][1],
	     minFrag);
}

//////////////////////////////////////////////////////////////////////

void Molecule::cut(int minFrag, vector<int> perm)
{
  // cut the bonds between the atoms that match the patterns,
  // cleaving the bonds in the order specified in 'perm',

  for (unsigned int i = 0; i < perm.size(); i++) {
    cleaveBond(matchPairs[perm[i]][0], matchPairs[perm[i]][1],
	       minFrag);
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::dealWithAmbiguousCuts()
{
  // deal with the problem of ambiguous cuts, due to patterns
  // that match more than two atoms (e.g., amine, etc.)
  // N.B. this function is used only if no extensive fragmentation
  // is to be carried out

  // count in how many pairs an atom shows up
  vector<unsigned int> tally(mol.NumAtoms(), 0);
  
  for (unsigned int i = 0; i < matchPairs.size(); i++) {
    tally[matchPairs[i][0]]++;
    tally[matchPairs[i][1]]++;
  }

  // go over each atom pair -- if an atom shows up in multiple
  // fragments, let it end up only in the smallest one
  unsigned int currSize = 0, currPair, size = 0, atom1, atom2;
  bool inAnotherPair = false;

  for (unsigned int i = 0; i < matchPairs.size(); i++) {
    atom1 = matchPairs[i][0];
    atom2 = matchPairs[i][1];

    // one of the atoms is in more than one pair
    if (tally[atom1] > 1 || tally[atom2] > 1) {
      size = fragmentSize(atom1, atom2);
      if (!inAnotherPair) {
	inAnotherPair = true;
	currSize = size;
	currPair = i;
      }
      else {
	// remove the pair that yields the larger fragment
	if (size < currSize) {
	  matchPairs.erase(matchPairs.begin() + i);
	}
	else {
	  matchPairs.erase(matchPairs.begin() + currPair);
	  currPair = i;
	  currSize = size;
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::findCleavableBonds(vector<OBSmartsPattern> rules,
				  int minFrag) {
  // apply the cleavage rules and find the matching atom pairs

  for (vector<OBSmartsPattern>::iterator rule = rules.begin();
       rule != rules.end(); rule++) {

    if ((*rule).Match(mol)) { // the pattern matches

      // find all pairs of matching atoms
      vector<vector<int> > maplist;
      maplist = (*rule).GetMapList();
      int atom1, atom2;

      for (unsigned int i = 0; i < maplist.size(); i++) {
	// get the resulting fragment sizes
	int sizeFrag1 =
	  Molecule::fragmentSize(maplist[i][0], maplist[i][1]);
	int sizeFrag2 =
	  Molecule::fragmentSize(maplist[i][1], maplist[i][0]);

	// fragment sizes are OK
	if (sizeFrag1 >= minFrag && sizeFrag2 >= minFrag) {
          atom1 = maplist[i][0];
	  atom2 = maplist[i][1];

	  // add the atoms to the list of atoms
	  vector<int> pair;
	  pair.push_back(atom1);
	  pair.push_back(atom2);
	  matchPairs.push_back(pair);
	}
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::fragment(vector<OBSmartsPattern> rules,
			OBConversion conv, Parameters p)
{
  // fragment the molecule according to 'rules'

  // find the cleavable bonds
  findCleavableBonds(rules, p.minFrag);

  // proceed if there is at least one cleavable bond
  if (matchPairs.size() > 0) {

    // extensive fragmentation
    if (p.extensive) {

      // first, make a copy of the molecule
      original = mol;
 
      // cluster the bonds to be cleaved into mutually exclusive sets
      getIndependentSets(p.minFrag);

      // get the singleton bonds
      getSingletons();

      // carry out the fragmentation for each independent set
      for (unsigned int i = 0; i < independentSets.size(); i++) {
	mol = original;
	cut(p.minFrag, independentSets[i]);
	storeFragments(conv);
      }

      // Take care of the singleton bonds
      for (unsigned int i = 0; i < singletons.size(); i++) {
	mol = original;
	cut(p.minFrag, singletons[i]);
	storeFragments(conv);
      }

      //  g.free();
    }

    // simple fragmentation
    else {

      // cut bonds, if match
      dealWithAmbiguousCuts();
      cut(p.minFrag);
      storeFragments(conv);
    }
  }

  else { // no cleavable bonds, just print the results
    storeFragments(conv);
  }
}

//////////////////////////////////////////////////////////////////////

int Molecule::fragmentSize(int posAtom1, int posAtom2)
{
  // calculate the size of the fragment comprised between atom1
  // and atom2  

  // get the atoms between atom1 and atom2
  vector<int> children;
  mol.FindChildren(children, posAtom1, posAtom2);

  // make sure we count only heavy atoms (i.e., ignore H)
  OBAtom *atom;
  int fragSize = 0;
  for (unsigned int i = 0; i < children.size(); i++) {
    atom = mol.GetAtom(children[i]);
    if (! atom->IsHydrogen()) {
      fragSize++;
    }
  }
  
  return fragSize + 1;
}

//////////////////////////////////////////////////////////////////////

void Molecule::freeDepMat()
{
  // free the memory allocated for the distance matrix

  unsigned int n = matchPairs.size();

  for (unsigned int i = 0; i < n; i++)
    free(depMat[i]);
  free(depMat);

}

//////////////////////////////////////////////////////////////////////

void Molecule::getIndependentSets(int minFrag)
{
  // calculate the graph theoretic distance (shortest number of
  // bonds) between all atoms, and cluster the atoms pairs into either
  // dependent (if their distance < minFrag) or independent groups
  // using complete linkage

  // build the dependence matrix between the cleavable bonds
  calculateDepMat(minFrag);

  // calculate the independent sets with the Bron-Kerbosch algorithm
  vector<int> R, P, X;
  for (unsigned int i = 0; i < matchPairs.size(); i++) {
    P.push_back(i);
  }
  bronKerbosch(R, P, X);

  // free memory for the distance matrix
  freeDepMat();

}

//////////////////////////////////////////////////////////////////////

void Molecule::getSingletons()
{
  // get the bonds that are not part of any independent set

  vector<bool> isIncluded(matchPairs.size(), false);

  // check which bonds have been included
  for (unsigned int i = 0; i < independentSets.size(); i++) {
    for (unsigned int j = 0; j < independentSets[i].size(); j++) {
      isIncluded[independentSets[i][j]] = true;
    }
  }

  // put the remaining bonds in the "singletons" vector
  for (unsigned int i = 0; i < isIncluded.size(); i++) {
    if (!isIncluded[i]) {
      singletons.push_back(i);
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::printResults(fstream &outFile, OBConversion conv,
			    unsigned int minFrag)
{
  // print the molecule to file

  // further remove redundancy from fragments by converting them
  // to molecules and back to strings

  vector<string> nonRedundant;
  OBMol temp;
  string currFrag;

  for (vector<string>::iterator fragment = fragments.begin();
       fragment != fragments.end(); fragment++) {
    conv.ReadString(&temp, *fragment);
    if (temp.NumAtoms() >= minFrag) {
      currFrag = conv.WriteString(&temp, true);
      if (find(nonRedundant.begin(), nonRedundant.end(), currFrag) ==
	  nonRedundant.end()) {
	nonRedundant.push_back(currFrag);
      }
    }
  }

  bool needDot = false;
  bool notEmpty = false;
  if (nonRedundant.size() > 0) {
    notEmpty = true;
  }
  for (vector<string>::iterator fragment = nonRedundant.begin();
       fragment != nonRedundant.end(); fragment++) {
    if (needDot) {
      outFile << "." << *fragment;
    }
    else {
      outFile << *fragment;
      needDot = true;
    }
  }

  // print the name (if not empty)
  if (notEmpty) {
    if (!name.empty()) {
      outFile << "\t" << name << endl;
    }
    else {
      outFile << endl;
    }
  }
}

//////////////////////////////////////////////////////////////////////

void Molecule::storeFragments(OBConversion conv)
{
  // store the fragments (if not already stored)

  // convert the molecule into a string
  string molString = conv.WriteString(&mol, true);

  // the '.' character separate the fragments...
  istringstream iss(molString);
  string fragment;
  OBMol temp;
  while (getline(iss, fragment, '.')) {
    // include only "new" fragments
    if (find(fragments.begin(), fragments.end(), fragment) ==
	fragments.end()) {

      fragments.push_back(fragment);
    }
  }
}

//////////////////////////////////////////////////////////////////////

Parameters::Parameters(char **argv, int argc)
{
  // parse the command-line arguments

  rulesFileName = getCmdOption(argv, argv + argc, "-r");
  smiFileName = getCmdOption(argv, argv + argc, "-i");
  outFileName = getCmdOption(argv, argv + argc, "-o");
  stringstream temp(getCmdOption(argv, argv + argc, "-n"));
  temp >> minFrag;

  // sanity check
  if (minFrag < 2) {
    cerr << "the minimum fragment size (-n) should be > 1\n";
    exit(1);
  }

  extensive = cmdOptionExists(argv, argv + argc, "-e");
}

//////////////////////////////////////////////////////////////////////

vector<OBSmartsPattern> storeCleavageRules(const char *rulesFileName)
{
  // store the cleavage rules as a vector of SMARTS PATTERN objects

  // declare the pattern vector
  vector<OBSmartsPattern> rules;

  // open the rules file
  fstream rulesFile;
  rulesFile.open(rulesFileName, fstream::in);

  // complain if the file doesn't exist
  if (! rulesFile.good()) {
    cerr << "Can't open " << rulesFileName << endl;
    exit(1);
  }

  // store the rules
  string line, rule;
  unsigned int count = 0;
  while (getline(rulesFile, line)) {

    // add an element to the pattern vectors
    rules.push_back(OBSmartsPattern());

    // process the pattern
    istringstream iss(line);
    iss >> rule;
    rules[count++].Init(rule);
  }

  rulesFile.close();

  return rules;
}

//////////////////////////////////////////////////////////////////////

vector<OBSmartsPattern> storeCleavageRules()
{
  // store the cleavage rules as a vector of SMARTS PATTERN objects

  // declare the pattern vector
  vector<OBSmartsPattern> rules;


  // store the rules
  string line, rule;
  for (unsigned int i = 0; i < NUM_RECAP; i++) {

    // add an element to the pattern vectors
    rules.push_back(OBSmartsPattern());

    // process the pattern
    rules[i].Init(RECAP[i]);
  }

  return rules;
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
  vector<OBSmartsPattern> rules;
  if (p.rulesFileName) {  // load the SMARTS cleavage rulse

    rules = storeCleavageRules(p.rulesFileName);
  }
  else { // use default RECAP rules
    rules = storeCleavageRules();
  }

  // process the molecules
  fragAllMols(p, rules);
  cout << endl;

  return 0;
}
