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

#ifndef POLY_FIT_ENCODING_HPP
#define POLY_FIT_ENCODING_HPP

#include <vector>
#include <cstdint>
#include <cstddef>

#include "byte_buffer.hpp"

struct Frame_t {
    uint8_t num_of_elements;
    std::vector<double> parameters;

    Frame_t()
	: num_of_elements(0),
	  parameters() {}

    Frame_t(uint8_t a_num_of_elements,
	    std::vector<double> a_parameters)
	: num_of_elements(a_num_of_elements),
	  parameters(a_parameters) {}

    void write_to_buffer(Byte_buffer_t & output_buffer);
    void read_from_buffer(Byte_buffer_t & input_buffer);

    bool is_encoded_as_a_polynomial() const {
	return num_of_elements > parameters.size();
    }
};

typedef std::vector<Frame_t> Vector_of_frames_t;

void poly_fit_encode(const std::vector<double> & values,
		     size_t elements_per_frame,
		     unsigned int num_of_parameters,
		     double max_abs_error,
		     Byte_buffer_t & output_buffer,
		     size_t & num_of_frames,
		     size_t & num_of_frames_encoded_directly);
void poly_fit_decode(size_t num_of_elements_to_decode,
		     Byte_buffer_t & input_buffer,
		     std::vector<double> & values);

#endif
