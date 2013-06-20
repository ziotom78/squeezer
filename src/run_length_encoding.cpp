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

#include <algorithm>
#include <vector>

#include "byte_buffer.hpp"
#include "run_length_encoding.hpp"

void
rle_compression(const uint32_t * input_stream,
		size_t input_size,
		Byte_buffer_t & output_stream)
{
    const uint32_t * ptr_to_cur_dword = input_stream;

    while((size_t) (ptr_to_cur_dword - input_stream + 1) <= input_size) {
	uint32_t first_dword_in_the_run = *ptr_to_cur_dword++;

	uint32_t count = 1;
	while(*ptr_to_cur_dword == first_dword_in_the_run &&
	      count < UINT32_MAX &&
	      (size_t) (ptr_to_cur_dword - input_stream + 1) < input_size) {
	    ptr_to_cur_dword++;
	    count++;
	}

	output_stream.append_uint32(count);
	output_stream.append_uint32(first_dword_in_the_run);
    }    
}

//////////////////////////////////////////////////////////////////////

void
rle_decompression(Byte_buffer_t & input_stream,
		  size_t output_size,
		  std::vector<uint32_t> & output)
{
    output.resize(output_size);
    size_t output_idx = 0;
    while(output_idx < output_size) {
	uint32_t count = input_stream.read_uint32();
	uint32_t value = input_stream.read_uint32();

	std::fill(output.begin() + output_idx,
		  output.begin() + output_idx + count,
		  value);

	output_idx += count;
    }
}
