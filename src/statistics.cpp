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

#include <vector>
#include <cmath>
#include "statistics.hpp"

////////////////////////////////////////////////////////////////////////

void
build_frequency_table(const bytestream_t & bytestream,
		      frequency_table_t & freq_table)
{
    for(auto cur_byte : bytestream) {

	if(freq_table.find(cur_byte) != freq_table.end()) {
	    freq_table[cur_byte]++;
	} else {
	    freq_table[cur_byte] = 1;
	}
    }
}

////////////////////////////////////////////////////////////////////////

static size_t
num_of_symbols(const frequency_table_t & freq_table)
{
    size_t number_of_symbols = 0;
    for(auto & i : freq_table) {

	number_of_symbols += i.second;
    }
    
    return number_of_symbols;
}

////////////////////////////////////////////////////////////////////////

double entropy_from_frequency_table(const frequency_table_t & freq_table)
{
    size_t num = num_of_symbols(freq_table);
    if(num == 0)
	return 0.0;

    const double inv_num = 1.0 / num;
    double entropy = 0.0;
    for(auto & i : freq_table) {

	double prob = i.second * inv_num;
	entropy -= prob * std::log2(prob);
    }
    
    return entropy;
}
