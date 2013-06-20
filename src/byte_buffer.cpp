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
#include <stdexcept>
#include <cassert>
#include "byte_buffer.hpp"

//////////////////////////////////////////////////////////////////////

Byte_buffer_t::Byte_buffer_t(size_t length, const uint8_t * raw_buffer)
{
    buffer.resize(length);
    std::copy(raw_buffer, raw_buffer + length,
	      buffer.begin());

    cur_position = 0;
}

//////////////////////////////////////////////////////////////////////

uint16_t
Byte_buffer_t::read_uint16()
{
    uint8_t byte1 = read_uint8();
    uint8_t byte2 = read_uint8();
    return (((uint16_t) byte1) << 8) + byte2;
}

//////////////////////////////////////////////////////////////////////

uint32_t
Byte_buffer_t::read_uint32()
{
    uint16_t word1 = read_uint16();
    uint16_t word2 = read_uint16();
    return (((uint32_t) word1) << 16) + word2;
}

//////////////////////////////////////////////////////////////////////

uint64_t
Byte_buffer_t::read_uint64()
{
    uint32_t dword1 = read_uint32();
    uint32_t dword2 = read_uint32();
    return (((uint64_t) dword1) << 32) + dword2;
}

//////////////////////////////////////////////////////////////////////

float
Byte_buffer_t::read_float()
{
    static_assert(sizeof(float) == 4, 
		  "This code assumes that single-precision floating "
		  "point values are 4 bytes wide.");

    uint32_t value = read_uint32();
    return *(reinterpret_cast<float *>(&value));
}

//////////////////////////////////////////////////////////////////////

double
Byte_buffer_t::read_double()
{
    static_assert(sizeof(double) == 8, 
		  "This code assumes that double-precision floating "
		  "point values are 8 bytes wide.");

    uint64_t value = read_uint64();
    return *(reinterpret_cast<double *>(&value));
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::read_buffer(size_t length, uint8_t * raw_buffer)
{
    if(cur_position + length >= buffer.size()) {
	throw std::out_of_range("Byte_buffer_t::read_buffer asked "
				"for too much data");
    }

    std::copy(buffer.begin() + cur_position,
	      buffer.begin() + cur_position + length,
	      raw_buffer);

    cur_position += length;
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::append_uint16(uint16_t value)
{
    uint8_t byte1 = (value & 0xFF00U) >> 8;
    uint8_t byte2 = value & 0xFFU;

    append_uint8(byte1);
    append_uint8(byte2);
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::append_uint32(uint32_t value)
{
    uint16_t word1 = (value & 0xFFFF0000U) >> 16;
    uint16_t word2 = value & 0xFFFFU;

    append_uint16(word1);
    append_uint16(word2);
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::append_uint64(uint64_t value)
{
    uint32_t dword1 = (value & 0xFFFFFFFF00000000UL) >> 32;
    uint32_t dword2 = value & 0xFFFFFFFFUL;

    append_uint32(dword1);
    append_uint32(dword2);
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::append_float(float value)
{
    static_assert(sizeof(float) == 4, 
		  "This code assumes that single-precision floating "
		  "point values are 4 bytes wide.");

    uint32_t int_value = *(reinterpret_cast<uint32_t *>(&value));
    append_uint32(int_value);
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::append_double(double value)
{
    static_assert(sizeof(double) == 8, 
		  "This code assumes that single-precision floating "
		  "point values are 8 bytes wide.");

    uint64_t int_value = *(reinterpret_cast<uint64_t *>(&value));
    append_uint64(int_value);
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::append_data_from_buffer(size_t length, const uint8_t * buffer)
{
    for(size_t idx = 0; idx < length; ++idx)
	append_uint8(buffer[idx]);
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::append_data_from_file(FILE * input, size_t length)
{
    for(size_t idx = 0; idx < length; ++idx) {
	int character = fgetc(input);
	if(character == EOF)
	    throw std::runtime_error("unexpected end of file");

	append_uint8(character);
    }
}

//////////////////////////////////////////////////////////////////////

void
Byte_buffer_t::write_to_file(FILE * out)
{
    assert(out != NULL);
    if(fwrite(buffer.data(), 1, buffer.size(), out) < buffer.size()) {
	throw std::runtime_error("unable to write the byte buffer to the file");
    }
}
