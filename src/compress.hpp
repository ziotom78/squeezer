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

#ifndef COMPRESS_HPP
#define COMPRESS_HPP

#include <cstdio>
#include <cstdint>
#include "common_defs.hpp"

class Detector_pointings_t;

struct Compression_parameters_t {
    Squeezer_file_type_t file_type;
    Radiometer_t radiometer;
    uint16_t od_number;
    size_t elements_per_frame;
    unsigned int number_of_poly_terms;
    double max_abs_error;
    bool read_calibrated_data;
    bool verbose_flag;

    Compression_parameters_t()
	: file_type(SQZ_NO_DATA),
	  radiometer(),
	  od_number(0),
	  elements_per_frame(25),
	  number_of_poly_terms(3),
	  max_abs_error(10.0 / 3600.0 * M_PI / 180.0),
	  read_calibrated_data(false),
	  verbose_flag(false) {}
};

void
compress_file_to_file(const std::string & input_file_name,
		      FILE * output_file,
		      const Compression_parameters_t & params);

#endif
