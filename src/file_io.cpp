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

#include <cstdio>
#include <cstdint>
#include <cstring>
#include <stdexcept>

#include "file_io.hpp"

//////////////////////////////////////////////////////////////////////

uint8_t
read_uint8(FILE * in)
{
    int result = fgetc(in);
    if(result == EOF)
	throw std::runtime_error(std::strerror(errno));

    return (uint8_t) result;
}

//////////////////////////////////////////////////////////////////////

uint16_t
read_uint16(FILE * in)
{
    uint8_t byte1 = read_uint8(in);
    uint8_t byte2 = read_uint8(in);

    return (((uint16_t) byte1) << 8) + byte2;
}

//////////////////////////////////////////////////////////////////////

uint32_t
read_uint32(FILE * in)
{
    uint16_t word1 = read_uint16(in);
    uint16_t word2 = read_uint16(in);

    return (((uint32_t) word1) << 16) + word2;
}

//////////////////////////////////////////////////////////////////////

uint64_t
read_uint64(FILE * in)
{
    uint32_t double_word1 = read_uint32(in);
    uint32_t double_word2 = read_uint32(in);

    return (((uint64_t) double_word1) << 32) + double_word2;
}

//////////////////////////////////////////////////////////////////////

double
read_double(FILE * in)
{
    double value;
    if(fread((void *) &value, sizeof(value), 1, in) < 1) {
	throw std::runtime_error(std::strerror(errno));
    }

    return value;
}

//////////////////////////////////////////////////////////////////////

void
write_uint8(FILE * out, uint8_t value)
{
    if(fputc(value, out) == EOF)
	throw std::runtime_error(std::strerror(errno));
}

//////////////////////////////////////////////////////////////////////

void
write_uint16(FILE * out, uint16_t value)
{
    write_uint8(out, (value & 0xFF00) >> 8); 
    write_uint8(out, (uint8_t) (value & 0xFF)); 
}

//////////////////////////////////////////////////////////////////////

void
write_uint32(FILE * out, uint32_t value)
{
    write_uint16(out, (value & 0xFFFF0000) >> 16); 
    write_uint16(out, (uint16_t) (value & 0xFFFF)); 
}

//////////////////////////////////////////////////////////////////////

void
write_uint64(FILE * out, uint64_t value)
{
    write_uint32(out, (value & 0xFFFFFFFF00000000) >> 32);
    write_uint32(out, (uint32_t) (value & 0xFFFFFFFF)); 
}

//////////////////////////////////////////////////////////////////////

void
write_double(FILE * out, double value)
{
    if(fwrite((const void *) &value, sizeof(value), 1, out) < 1) {
	throw std::runtime_error(std::strerror(errno));
    }
}
