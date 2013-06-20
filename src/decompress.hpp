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

#ifndef DECOMPRESS_HPP
#define DECOMPRESS_HPP

#include <cstdio>
#include <string>

struct Decompression_parameters_t {
    bool verbose_flag;

    Decompression_parameters_t() {
	verbose_flag = false;
    }
};

void decompress_detpoints_from_file(FILE * input_file,
				    const std::string & output_file_name,
				    const Decompression_parameters_t & params);

#endif
