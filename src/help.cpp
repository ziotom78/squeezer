/*
 * Squeezer - compress LFI detector pointings and differenced data
 * Copyright (C) 2013 Maurizio Tomasi (Planck collaboration)
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA.
 */

#include <vector>
#include <string>
#include <iostream>

#include "common_defs.hpp"
#include "help.hpp"

const char * help_text_general =
    "Usage: squeezer COMMAND [parameters...]\n"
    "\n"
    "An invocation to \"squeezer\" can have one of the following forms:\n"
    "\n"
    "    squeezer compress [options] RADIOMETER OD INPUT_FILE OUTPUT_FILE\n"
    "    squeezer compress [options] PARAMETER_FILE\n"
    "    squeezer decompress [options] INPUT_FILE OUTPUT_FILE\n"
    "    squeezer statistics [options] INPUT_FILE OUTPUT_FILE\n"
    "    squeezer help [COMMAND]\n"
    "    squeezer help version\n"
    "\n"
    "To get more help about COMMAND, run \"squeezer help COMMAND\".\n"
    "\n"
    "This program has been written by Maurizio Tomasi <tomasi@lambrate.inaf.it>.\n";

const char * help_text_compress =
    "Usage: squeezer compress [options] RADIOMETER OD INPUT_FILE OUTPUT_FILE\n"
    "   or: squeezer compress [options] PARAMETER_FILE\n"
    "\n"
    "Compress data and save them into a binary file.\n"
    "\n"
    "The name of the RADIOMETER must be in the form LFInna,\n"
    "with \"nn\" in the 18..28 range and \"a\" either \"M\" or \"S\".\n"
    "\n"
    "The number OD is the number of the operational day.\n"
    "\n"
    "In the form using PARAMETER_FILE, the values for RADIOMETER, OD, \n"
    "INPUT_FILE, and OUTPUT_FILE are specified in a text file. Every\n"
    "line that does not start with \"#\" and is not empty is interpreted\n"
    "as a sequence of the form \"RADIOMETER OD INPUT_FILE OUTPUT_FILE\".\n"
    "\n"
    "The following is an example of a text file:\n"
    "\n"
    "   # This is a comment and is ignored\n"
    "   LFI18M 91 LFI18M_0091_pointings.fits 18M_0091.pntz\n"
    "   LFI18M 92 LFI18M_0092_pointings.fits 18M_0092.pntz\n"
    "\n"
    "Data are read from INPUT_FILE file, which can be either a FITS\n"
    "file or a DMC object. The code determines the data source\n"
    "depending on the following rules:\n"
    "\n"
    "   * If it begins with TOODI%, it is a DMC object\n"
    "   * In any other case, it is a FITS file. CFITSIO extended\n"
    "     syntax will work (e.g. appending [N] to the file name\n"
    "     will open the HDU number N).\n"
    "\n"
    "OUTPUT_FILE can be a minus sign (\"-\"), in which case the file\n"
    "is written to standard output. This allows redirection and piping.\n"
    "\n"
    "Possible options are:\n"
    "\n"
    "   --pointings Assume that the input data are detector pointings (default).\n"
    "   --datadiff  Assume that the input data are differenced voltages.\n"
    "   -n NUM      When compressing angles, this specifies the number of\n"
    "               elements in a \"frame\". This value must always be\n"
    "               greater than the one specified using -p.\n"
    "   -p NUM      When compressing angles, this specifies the order of\n"
    "               the interpolating polynomial. This value must always\n"
    "               be smaller than the one specified using -n.\n"
    "   -s NUM      When compressing angles, this specifies the maximum\n"
    "               error in arcseconds between the angles and the interpolating\n"
    "               polynomial. Every time this value is overcame in a frame,\n"
    "               compression will be turned off for that frame. This\n"
    "               prevents compression errors from getting too big.\n"
    "   -v       Be verbose.\n";

const char * help_text_decompress =
    "Usage: squeezer decompress INPUT_FILE OUTPUT_FITS_FILE\n"
    "\n"
    "Decompress a binary file into a FITS file. This is the\n"
    "opposite of \"squeezer compress\". A few caveats:\n"
    "\n"
    "  1. OUTPUT_FITS_FILE cannot be \"-\" (redirection to standard\n"
    "     output), because of limitations in the FITS file format.\n"
    "  2. Compressing a file and then decompressing it will not\n"
    "     produce the same file. Some of the compression algorithms\n"
    "     used by \"squeezer\" are lossy, and therefore some information\n"
    "     gets lost.\n"
    "  3. To quantify the amount of compression and how much information\n"
    "     has been lost, use \"squeezer statistics\". (Run the command\n"
    "     \"squeezer help statistics\" for more information.)\n"
    "\n"
    "Possible options are:\n"
    "\n"
    "   -v      Be verbose.\n";

const char * help_text_statistics =
    "Usage: squeezer statistics [options] BINARY_FILE\n"
    "\n"
    "Show some statistics about a compressed binary file (created\n"
    "using \"squeezer compress\"). If BINARY_FILE is a minus sign (\"-\")\n"
    "then the file is read from standard input. This allows to use\n"
    "redirection and piping.\n"
    "\n"
    "Possible options are:\n"
    "\n"
    "   -html    Output a report in HTML format.\n";

const char * help_text_help =
    "Print command-line help.\n";

//////////////////////////////////////////////////////////////////////

void
print_help_on_compress_command()
{
    std::cout << help_text_compress;
}

//////////////////////////////////////////////////////////////////////

void
print_help_on_decompress_command()
{
    std::cout << help_text_decompress;
}

//////////////////////////////////////////////////////////////////////

void
print_help_on_statistics_command()
{
    std::cout << help_text_statistics;
}

//////////////////////////////////////////////////////////////////////

void
print_help_on_help_command()
{
    std::cout << help_text_help;
}

//////////////////////////////////////////////////////////////////////

void
print_help(const std::vector<std::string> & list_of_arguments)
{
    if(list_of_arguments.size() < 2) {
	// Generic help
	std::cout << help_text_general;
    } else {
	// More help on a given command
	if(list_of_arguments.at(1) == "compress") {

	    print_help_on_compress_command();

	} else if(list_of_arguments.at(1) == "decompress") {

	    print_help_on_decompress_command();

	} else if(list_of_arguments.at(1) == "statistics") {

	    print_help_on_statistics_command();

	} else if(list_of_arguments.at(1) == "help") {

	    print_help_on_help_command();

	} else if(list_of_arguments.at(1) == "version") {

	    std::cout << PROGRAM_NAME
		      << ' '
		      << reinterpret_cast<int>((PROGRAM_VERSION & 0xFF00) >> 8)
		      << '.'
		      << reinterpret_cast<int>(PROGRAM_VERSION & 0xFF)
		      << '\n';

	} else {

	    std::cerr << PROGRAM_NAME 
		      << ": unknown command \""
		      << list_of_arguments.at(1)
		      << "\"\n";
	    std::exit(1);

	}
    }
}
