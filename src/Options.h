/* options.h - base class for command line options

   Copyright (c) 2006 Frank Dudbridge
   MRC Biostatistics Unit
   Robinson Way
   Cambridge CB2 2SR, UK
   frank.dudbridge@mrc-bsu.cam.ac.uk

   This program is free software; you can redistribute it and/or
   modify it under the terms of the GNU General Public License
   as published by the Free Software Foundation; either version 2
   of the License, or (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
   USA.

*/

#include <iostream>
#include <strstream>
#include "stl.h"
#include <string>

using namespace std;

#ifndef __OPTIONS__

class OptionInfo {
public:
    OptionInfo(int m, bool u = false) {
        nargs = m;
        used = u;
    }
    int nargs;
    bool used;
};

class Options {
public:
    Options(int c, char **v, bool b = true) {
        argc = c;
        for (int i = 0; i < argc; i++) {
            argv.push_back(v[i]);
        }
        programName = argv[0];
        needfiles = b;
    }
    vector<string> inputfile;
protected:

    //  option

    // set up different kinds of option and read their values
    // can have an argv that is the intial substring of an option name
    // should not have two options with one a substring of the other

    // several generic arguments
    template<class Type>
    void option(string s, vector<Type> &arg) {
        OptionInfo o(-1);
        arg.clear();
        for (int i = 1; i < argc; i++) {
            if (s.find(argv[i]) == 0) {
                for (i++; i < argc && argv[i][0] != '-'; i++) {
                    istrstream istr(argv[i].c_str());
                    Type val;
                    istr >> val;
                    arg.push_back(val);
                }
                o.nargs = arg.size();
            }
        }
        optionList.insert(make_pair(s, o));
    }

    // several string arguments
    void option(string s, vector<string> &arg) {
        OptionInfo o(-1);
        arg.clear();
        for (int i = 1; i < argc; i++) {
            if (s.find(argv[i]) == 0) {
                for (i++; i < argc && argv[i][0] != '-'; i++) {
                    arg.push_back(argv[i]);
                }
                o.nargs = arg.size();
            }
        }
        optionList.insert(make_pair(s, o));
    }

    // single argument
    template<class Type1, class Type2>
    void option(string s, Type1 &arg, Type2 initial) {
        OptionInfo o(1);
        arg = (Type1)initial;
        for (int i = 1; i < argc; i++) {
            if (s.find(argv[i]) == 0) {
                i++;
                if (i < argc && argv[i][0] != '-') {
                    istrstream istr(argv[i].c_str());
                    istr >> arg;
                }
            }
        }
        optionList.insert(make_pair(s, o));
    }

    // single string argument
    template<class Type2>
    void option(string s, string &arg, Type2 initial) {
        OptionInfo o(1);
        arg = (string)initial;
        for (int i = 1; i < argc; i++) {
            if (s.find(argv[i]) == 0) {
                i++;
                if (i < argc && argv[i][0] != '-') {
                    arg = argv[i];
                }
            }
        }
        optionList.insert(make_pair(s, o));
    }

    // no argument - boolean flag
    void option(string s, bool &arg, bool initial) {
        OptionInfo o(0);
        arg = initial;
        for (int i = 1; i < argc; i++) {
            if (s.find(argv[i]) == 0) {
                arg = true;
            }
        }
        optionList.insert(make_pair(s, o));
    }



    //  parse

    // parse the command line and read arguments
    void parse() {
        //  inputfile.clear();
        for (int i = 1; i < argc; i++) {
            // print usage message if -help selected
            if ((new string("-help"))->find(argv[i]) == 0) {
                usage();
            }
            if (argv[i][0] == '-') {
                string s;
                s = argv[i];
                // find an exact match, or failing that an intial match
                vector<string> matches;
                OptionInfo *o;
                for (map<string, OptionInfo>::iterator j = optionList.begin(); j != optionList.end(); j++) {
                    if (j->first.find(argv[i]) == 0) {
                        matches.push_back(j->first);
                        o = &j->second;
                    }
                }
                if (matches.size() == 0) {
                    cout << endl << "ERROR: unrecognised option: " << argv[i] << endl;
                    usage();
                }
                if (matches.size() > 1) {
                    cout << endl << "ERROR: ambiguous option: " << argv[i] << " (matches:";
                    for (int j = 0; j < matches.size(); j++) {
                        cout << " " << matches[j];
                    }
                    cout << ")" << endl;
                    usage();
                }
                if (o->used) {
                    cout << endl << "ERROR: multiple use of option: " << matches[0] << endl;
                    usage();
                }

                // skip the arguments for this option
                o->used = true;
                for (int j = i + 1; j < i + o->nargs + 1; j++)
                    if (j >= argc || argv[j][0] == '-') {
                        cout << endl << "ERROR: argument missing for " << matches[0] << endl;
                        usage();
                    }
                i += o->nargs;
            } else {
                inputfile.push_back(argv[i]);
            }
        }
        if (inputfile.size() == 0 && needfiles) {
            cout << endl << "ERROR: no input files" << endl;
            usage();
        }
    }



    //  usage

    // print usage message
    void usage() {
        cout << endl << "Usage: " << programName;
        if (needfiles) {
            cout << " <inputfiles>";
        }
        cout << " [options]" << endl;
        cout << "Options:" << endl;
        int ix = 1;
        for (map<string, OptionInfo>::iterator i = optionList.begin(); i != optionList.end(); i++) {
            cout << "\t" << i->first;
            switch (i->second.nargs) {
            case 1: {
                cout << " <n>";
                break;
            }
            case 0: {
                break;
            }
            default: {
                cout << " <...>";
                break;
            }
            }
            if (!(ix % 4)) {
                cout << endl;
            }
            ix++;
        }
        cout << endl;
        exit(-1);

    }



    string programName;
    bool needfiles; // true if input file required
    map<string, OptionInfo> optionList;
    int argc;
    vector<string> argv;
};

#endif
#define __OPTIONS__

