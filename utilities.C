//////////////////////////////////////////////////////////////////////
// utilities.C  Copyright (c) 2014 Dario Ghersi                     //
// Version: 20140307                                                //
// Goal:    utility functions                                       //
//                                                                  //
// This file is part of the BLOCKS suite.                           //
// BLOCKS is free software: you can redistribute it and/or modify   //
// it under the terms of the GNU General Public License as          //
// published by the Free Software Foundation, either version 3 of   //
// the License, or (at your option) any later version.              //
//                                                                  //
// BLOCKS is distributed in the hope that it will be useful,        //
// but WITHOUT ANY WARRANTY; without even the implied warranty of   //
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    //
// GNU General Public License for more details.                     //
//                                                                  //
// You should have received a copy of the GNU General Public        //
// License along with BLOCKS.                                       //
// If not, see <http://www.gnu.org/licenses/>.                      //
//////////////////////////////////////////////////////////////////////

#include <algorithm>
#include <iostream>
#include <string>
using namespace std;

//////////////////////////////////////////////////////////////////////
// DEFINITIONS                                                      //
//////////////////////////////////////////////////////////////////////

bool cmdOptionExists(char **begin, char **end,
		     const std::string& option)
{
  return std::find(begin, end, option) != end;
}

//////////////////////////////////////////////////////////////////////

char *getCmdOption(char **begin, char **end,
		   const string & option)
{
  char **itr = std::find(begin, end, option);
  if (itr != end && ++itr != end) {
    return *itr;
  }

  return 0;
}

//////////////////////////////////////////////////////////////////////

void printProgBar(unsigned int percent) {
  std::string bar;

  for (unsigned int i = 0; i < 50; i++) {
    if (i < (percent / 2)) {
      bar.replace(i, 1, "=");
    }
    else if (i == (percent / 2)) {
      bar.replace(i, 1, ">");
    }
    else {
      bar.replace(i, 1, " ");
    }
  }

  cout << "\r" "[" << bar << "] ";
  cout.width(3);
  cout << percent << "%     " << std::flush;
}

//////////////////////////////////////////////////////////////////////
