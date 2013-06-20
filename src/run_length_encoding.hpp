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

#ifndef RUN_LENGTH_ENCODING
#define RUN_LENGTH_ENCODING

#include <cstdint>
#include <vector>
#include <cstddef>

#include "byte_buffer.hpp"

void rle_compression(const uint32_t * input_stream,
		     size_t input_size,
		     Byte_buffer_t & output_stream);

void rle_decompression(Byte_buffer_t & input_stream,
		       size_t output_size,
		       std::vector<uint32_t> & output);

#endif
