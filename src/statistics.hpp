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

#ifndef STATISTICS_HPP
#define STATISTICS_HPP

#include <cstdint>
#include <vector>
#include <map>
#include <cstddef>

typedef uint8_t byte_t;
typedef std::vector<byte_t> bytestream_t;
typedef std::map<byte_t, size_t> frequency_table_t;

void build_frequency_table(const bytestream_t & bytestream,
			   frequency_table_t & freq_table);
double entropy_from_frequency_table(const frequency_table_t & freq_table);

template<typename T> void
vector_to_bytestream(const std::vector<T> & vector,
		     bytestream_t & result)
{
    const size_t type_size = sizeof(T);
    result.resize(vector.size() * type_size);

    for(size_t idx = 0; idx < vector.size(); ++idx)
    {
	size_t byte_idx = idx * type_size;
	T cur_element = vector[idx];
	byte_t * chunk_bytes = reinterpret_cast<byte_t *>(&cur_element);
	for(size_t i = 0; i < type_size; ++i) {
	    result[byte_idx + i] = chunk_bytes[i];
	}
    }
}

template<typename T> double
calc_entropy(const std::vector<T> & vector)
{
    bytestream_t bytestream;
    vector_to_bytestream(vector, bytestream);

    frequency_table_t freq_table;
    build_frequency_table(bytestream, freq_table);

    return entropy_from_frequency_table(freq_table);
}

#endif
