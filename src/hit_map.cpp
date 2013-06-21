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

#include "config.hpp"

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <hpixlib/hpix.h>

#include "decompress.hpp"
#include "detpoint.hpp"

#ifdef PROGRAM_NAME
#undef PROGRAM_NAME
#endif

#define PROGRAM_NAME "hit_map"

//////////////////////////////////////////////////////////////////////

struct Configuration_t {
    hpix_nside_t nside;
    std::string output_map_name;
    std::vector<std::string> list_of_pointings;
    
    void parse_from_command_line(int argc, const char * argv[]) {
	hpix_nside_t user_nside;
	std::stringstream ss(argv[1]);
	ss >> user_nside;

	if(! hpix_valid_nside(user_nside)) {
	    std::cerr << PROGRAM_NAME 
		      << ": invalid value for NSIDE (" 
		      << user_nside 
		      << ")\n";
	    std::exit(1);
	} else
	    nside = user_nside;

	output_map_name = argv[2];
	for(int idx = 3; idx < argc; ++idx)
	    list_of_pointings.push_back(argv[idx]);
    }
};

//////////////////////////////////////////////////////////////////////

void
add_hits_from_pointings(hpix_map_t * map,
			const std::string & input_file_name)
{
    std::unique_ptr<Detector_pointings_t> detpoints;
    Decompression_parameters_t params;

    std::cerr << PROGRAM_NAME << ": reading file " << input_file_name << "\n";

#ifdef HAVE_TOODI
    if(input_file_name.substr(0, 6) == "TOODI%") {
        detpoints.reset(new Detector_pointings_t());
	detpoints->read_from_database(input_file_name);
    } else
#endif
    {
	FILE * input_file = std::fopen(input_file_name.c_str(), "rb");
	detpoints.reset(dynamic_cast<Detector_pointings_t *>(decompress_from_file(input_file, params)));
	std::fclose(input_file);
    }

    hpix_nside_t nside = hpix_map_nside(map);
    double * map_pixels = hpix_map_pixels(map);
    for(size_t idx = 0; idx < detpoints->theta.size(); ++idx) {
    
	hpix_pixel_num_t cur_pixel =
	    hpix_angles_to_ring_pixel(nside, 
				      detpoints->theta[idx],
				      detpoints->phi[idx]);
	
	map_pixels[cur_pixel]++;
	
    }
}

//////////////////////////////////////////////////////////////////////

int
main(int argc, const char * argv[])
{
    if(argc < 4) {

	std::cerr << "Usage: hit_map NSIDE MAP_FITS_FILE POINTING_FILE1 [...]\n";
	return 1;

    }

    Configuration_t config;
    config.parse_from_command_line(argc, argv);

    hpix_map_t * hit_map = hpix_create_map(config.nside, HPIX_ORDER_SCHEME_RING);

#ifdef HAVE_TOODI
    ObjectHandle init_handle;
    toodiInitializeLLIO("TOODI%file", &init_handle);
#endif

    for(auto & detpoints_file : config.list_of_pointings) {

	add_hits_from_pointings(hit_map, detpoints_file);

    }

    int status = 0;
    hpix_save_fits_component_to_file(config.output_map_name.c_str(),
				     hit_map, TDOUBLE, "Hits", &status);
    if(status != 0) {

	char fitsio_err_msg[FLEN_ERRMSG];
	fits_read_errmsg(fitsio_err_msg);
	std::cerr << PROGRAM_NAME 
		  << ": unable to save the map, reason is \""
		  << fitsio_err_msg
		  << "\"\n";
	return 1;

    }

#if HAVE_TOODI
    toodiCloseLLIO(init_handle);
#endif
}
