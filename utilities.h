//////////////////////////////////////////////////////////////////////
// utilities.h  Copyright (c) 2013 Dario Ghersi                     //
// Version: 20130730                                                //
//          See User's Guide for more details                       //
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

//////////////////////////////////////////////////////////////////////
// PROTOTYPES                                                       //
//////////////////////////////////////////////////////////////////////

bool cmdOptionExists(char **, char **, const string &);
char *getCmdOption(char **, char **, const string &);
void printProgBar(unsigned int);
