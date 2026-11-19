/***************************************************************************
*                                                                          *
*   Copyright (c) 2017                                                     *
*   FastFieldSolvers S.R.L.  http://www.fastfieldsolvers.com               *
*                                                                          *
*   This program is free software; you can redistribute it and/or modify   *
*   it under the terms of the GNU Lesser General Public License (LGPL)     *
*   as published by the Free Software Foundation; either version 2 of      *
*   the License, or (at your option) any later version.                    *
*   for detail see the LICENCE text file.                                  *
*                                                                          *
*   This program is distributed in the hope that it will be useful,        *
*   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
*   GNU Library General Public License for more details.                   *
*                                                                          *
*   You should have received a copy of the GNU Library General Public      *
*   License along with this program; if not, write to the Free Software    *
*   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307   *
*   USA                                                                    *
*                                                                          *
***************************************************************************/


#include "FasterCapConsole.h"

#ifdef FCG_HEADLESS

#include <iostream>
#include <wx/debug.h>

// failed asserts are printed, not shown in a dialog: wxWidgets opens one on Windows,
// also in console programs, and nobody closes it in a batch run
static void ConsoleAssertHandler(const wxString& file, int line, const wxString& func,
                                 const wxString& cond, const wxString& msg)
{
    std::cerr << file << "(" << line << "): assert \"" << cond << "\" failed in "
              << func << "(): " << msg << std::endl;
}

int main(int argc, char **argv)
{
    int i, ret;
    char bOption;

    wxSetAssertHandler(ConsoleAssertHandler);

    // the character after '-b' ('-b?', '-bv'), as FasterCapApp::OnInit() finds it for the GUI
    bOption = '\0';
    for(i=1; i<argc; i++) {
        if(argv[i][0] == '-' && argv[i][1] == 'b') {
            bOption = argv[i][2];
        }
    }

    FasterCapConsole console;

    ret = console.main(argc, argv, bOption);

    return ret;
}

#endif // FCG_HEADLESS
