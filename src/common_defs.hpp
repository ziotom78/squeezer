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

#ifndef COMMON_DEFS
#define COMMON_DEFS

#include <string>
#include <cstdint>

#define PROGRAM_NAME "squeezer"
#define PROGRAM_VERSION 0x0100

#define MAJOR_VERSION_FROM_UINT16(x) ((int) ((x) & 0xFF00) >> 8)
#define MINOR_VERSION_FROM_UINT16(x) ((int) (x) & 0xFF)

#define MAJOR_PROGRAM_VERSION MAJOR_VERSION_FROM_UINT16(PROGRAM_VERSION)
#define MINOR_PROGRAM_VERSION MINOR_VERSION_FROM_UINT16(PROGRAM_VERSION)

struct Radiometer_t {
    uint8_t horn;
    uint8_t arm;

    void parse_from_name(const std::string & name);
    bool is_valid() const {
	return (horn >= 18 && horn <= 28) && 
	    (arm == 0 || arm == 1);
    }

    std::string to_str() const;
};

enum Squeezer_file_type_t { 
    SQZ_NO_DATA,
    SQZ_DETECTOR_POINTINGS, 
    SQZ_DIFFERENCED_DATA 
};

enum Chunk_type_t { 
    CHUNK_DELTA_OBT = 10,
    CHUNK_SCET_ERROR = 11,
    CHUNK_THETA = 12,
    CHUNK_PHI = 13,
    CHUNK_PSI = 14,
    CHUNK_DIFFERENCED_DATA = 15,
    CHUNK_QUALITY_FLAGS = 16
};

#endif
