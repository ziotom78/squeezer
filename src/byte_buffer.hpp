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

#ifndef BYTE_BUFFER_HPP
#define BYTE_BUFFER_HPP

#include <cstdio>
#include <cstdint>
#include <cstddef>
#include <vector>

class Byte_buffer_t {
public:
    std::vector<uint8_t> buffer;
    size_t cur_position;

    Byte_buffer_t() 
	: buffer(),
	  cur_position(0) {}

    Byte_buffer_t(size_t length, const uint8_t * raw_buffer);

    size_t size() const {
	return buffer.size();
    }

    // Number of bytes that can still be read
    size_t items_left() const {
	return buffer.size() - cur_position;
    }

    uint8_t read_uint8() {
	return buffer.at(cur_position++);
    }

    uint16_t read_uint16();
    uint32_t read_uint32();
    uint64_t read_uint64();
    float read_float();
    double read_double();
    void read_buffer(size_t length, uint8_t * buffer);

    // Append operations do not update the current position
    void append_uint8(uint8_t value) {
	buffer.push_back(value);
    }

    void append_uint16(uint16_t value);
    void append_uint32(uint32_t value);
    void append_uint64(uint64_t value);
    void append_float(float value);
    void append_double(double value);
    void append_data_from_buffer(size_t length, const uint8_t * buffer);
    void append_data_from_file(FILE * input, size_t length);

    // Write the *whole* buffer to disk, irrespective of the current position
    void write_to_file(FILE * out);
};

#endif
