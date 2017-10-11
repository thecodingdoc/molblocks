//////////////////////////////////////////////////////////////////////
// fragment.h  Copyright (c) 2014 Dario Ghersi                      //
// Version: 20140307                                                //
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

const string USAGE = "fragmenter -r RULES -i SMALL_MOLECULES\n\
                    -o OUTFILE -n MIN_FRAGMENT_SIZE [-e]\n";

// default RECAP rules
const unsigned int NUM_RECAP = 11;
const char *RECAP[NUM_RECAP] = {
  "[$([C!$(C([#7])(=O)[!#1!#6])](=[O]))]!@[#7!$([#7][!#1!#6])]",
  "[$(C=!@O)]!@[$([O;+0])]",
  "[#6]!@[N;!$(N=*);!$(N[#6]=[!#6]);!$(N~[!#1!#6])!X4]",
  "[$(C(=!@O)([#7;+0;D2,D3])!@[#7;+0;D2,D3])]!@[$([#7;+0;D2,D3])]",
  "[O!$(O[#6]~[!#1!#6])]([#6])!@[#6]",
  "[#6]=!@[#6]",
  "[N;+1;D4]!@[#6]",
  "[$([n;+0])]-!@[#6!$([#6]=[!#6])]",
  "[$(N(@-C(=O)))]!@-[#6!$([#6]=[!#6])]",
  "[c]!@[c]",
  "[$([#7;+0;D2,D3])]-!@[$([S](=[O])=[O])]"
};

//////////////////////////////////////////////////////////////////////
// CLASSES                                                          //
//////////////////////////////////////////////////////////////////////

class Parameters {

 public:
  char *rulesFileName;
  char *smiFileName;
  char *outFileName;
  int minFrag;
  bool extensive;

  Parameters(char **, int);
};

//////////////////////////////////////////////////////////////////////

class Molecule {

 private:
  OBMol mol;
  OBMol original; // copy of the original (for extensive fragment.)
  string name;
  vector<vector<int> > matchPairs;
  vector<string> fragments;
  bool **depMat;
  vector< vector<int> > independentSets;
  vector<int> singletons;

  void bronKerbosch(vector<int>, vector<int>, vector<int>);
  void calculateDepMat(int);
  void getIndependentSets(int);
  void getSingletons();
  bool cleaveBond(int, int, int);
  void cut(int);
  void cut(int, int);
  void cut(int, vector<int>);
  void dealWithAmbiguousCuts();
  void findCleavableBonds(vector<OBSmartsPattern>, int);
  int fragmentSize(int, int);
  void freeDepMat();
  void storeFragments(OBConversion);

 public:
  Molecule (OBConversion, string, string);
  void fragment(vector<OBSmartsPattern>, OBConversion,
		Parameters);
  void getMatchAtoms();
  void printResults(fstream &, OBConversion, unsigned int);
};

//////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

void checkCommandLineArgs(char **, int);
void fragAllMols(Parameters, vector<OBSmartsPattern>);
vector<OBSmartsPattern> storeCleavageRules(const char *);
vector<OBSmartsPattern> storeCleavageRules();
