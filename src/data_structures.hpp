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

#ifndef DATA_STRUCTURES_HPP
#define DATA_STRUCTURES_HPP

#include <cstdint>
#include <cstdio>
#include "common_defs.hpp"

extern const uint16_t program_version;

struct Squeezer_file_header_t {
    uint8_t file_type_mark[4];
    double floating_point_check;

    uint16_t program_version;

    uint16_t date_year;
    uint8_t date_month;
    uint8_t date_day;

    uint8_t time_hour;
    uint8_t time_minute;
    uint8_t time_second;

    Radiometer_t radiometer;

    uint16_t od;

    // Since these values are needed to decompress the SCET chunk, it
    // would be better to gather them into a dedicated sub-structure.
    double first_obt;
    double last_obt;
    double first_scet_in_ms;
    double last_scet_in_ms;

    uint32_t number_of_chunks;

    Squeezer_file_header_t(Squeezer_file_type_t type);

    void read_from_file(FILE * in);
    void write_to_file(FILE * out) const;

    bool is_valid() const;
    bool is_compatible_version() const;
};

//////////////////////////////////////////////////////////////////////

struct Error_t {
    double min_abs_error;
    double max_abs_error;
    double mean_abs_error;
    double mean_error;

    Error_t();

    void read_from_file(FILE * in);
    void write_to_file(FILE * out) const;

    bool is_valid() const;
};

//////////////////////////////////////////////////////////////////////

struct Squeezer_chunk_header_t {
    uint8_t chunk_mark[4];
    uint64_t number_of_bytes;
    uint32_t number_of_samples;

    uint32_t chunk_type;

    Error_t compression_error;

    Squeezer_chunk_header_t();

    void read_from_file(FILE * in);
    void write_to_file(FILE * out) const;

    bool is_valid() const;
};

#endif
