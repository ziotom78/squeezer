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
#include <stdexcept>
#include <cstring>
#include <ctime>

#include "file_io.hpp"
#include "common_defs.hpp"
#include "data_structures.hpp"

//////////////////////////////////////////////////////////////////////

Detpoint_file_header_t::Detpoint_file_header_t()
{
    file_type_mark[0] = 'P';
    file_type_mark[1] = 'D';
    file_type_mark[2] = 'P';
    file_type_mark[3] = 0;

    program_version = PROGRAM_VERSION;
    floating_point_check = 2.3125e+5;

    std::time_t time_now = std::time(NULL);
    std::tm * tm = std::gmtime(&time_now);

    date_year = tm->tm_year + 1900;
    date_month = tm->tm_mon;
    date_day = tm->tm_mday;

    time_hour = tm->tm_hour;
    time_minute = tm->tm_min;
    time_second = tm->tm_sec;

    radiometer.horn = 0;
    radiometer.arm = 0;
    od = 0;
    first_obt = 0.0;
    last_obt = 0.0;
    first_scet_in_ms = 0.0;
    last_scet_in_ms = 0.0;
    number_of_chunks = 0;
}

//////////////////////////////////////////////////////////////////////

void 
Detpoint_file_header_t::read_from_file(FILE * in)
{
    file_type_mark[0] = read_uint8(in);
    file_type_mark[1] = read_uint8(in);
    file_type_mark[2] = read_uint8(in);
    file_type_mark[3] = read_uint8(in);

    floating_point_check = read_double(in);

    date_year = read_uint16(in);
    date_month = read_uint8(in);
    date_day = read_uint8(in);

    time_hour = read_uint8(in);
    time_minute = read_uint8(in);
    time_second = read_uint8(in);

    radiometer.horn = read_uint8(in);
    radiometer.arm = read_uint8(in);
    od = read_uint16(in);
    first_obt = read_double(in);
    last_obt = read_double(in);
    first_scet_in_ms = read_double(in);
    last_scet_in_ms = read_double(in);
    number_of_chunks = read_uint32(in);
}

//////////////////////////////////////////////////////////////////////

void 
Detpoint_file_header_t::write_to_file(FILE * out) const
{
    write_uint8(out, file_type_mark[0]);
    write_uint8(out, file_type_mark[1]);
    write_uint8(out, file_type_mark[2]);
    write_uint8(out, file_type_mark[3]);

    write_double(out, floating_point_check);

    write_uint16(out, date_year);
    write_uint8(out, date_month);
    write_uint8(out, date_day);

    write_uint8(out, time_hour);
    write_uint8(out, time_minute);
    write_uint8(out, time_second);

    write_uint8(out, radiometer.horn);
    write_uint8(out, radiometer.arm);
    write_uint16(out, od);
    write_double(out, first_obt);
    write_double(out, last_obt);
    write_double(out, first_scet_in_ms);
    write_double(out, last_scet_in_ms);
    write_uint32(out, number_of_chunks);
}

//////////////////////////////////////////////////////////////////////

bool
Detpoint_file_header_t::is_valid() const
{
    if(file_type_mark[0] != 'P' ||
       file_type_mark[1] != 'D' ||
       file_type_mark[2] != 'P' ||
       file_type_mark[3] != 0 ||
       date_year < 2013 ||
       date_month < 1 || date_month > 12 ||
       date_day < 1 || date_day > 31 ||
       time_hour > 23 ||
       time_minute > 59 ||
       time_second > 59 ||
       floating_point_check != 231250.0 ||
       (! radiometer.is_valid()) ||
       first_obt >= last_obt ||
       first_scet_in_ms >= last_scet_in_ms ||
       number_of_chunks == 0)
	return false;

    return true;
}

//////////////////////////////////////////////////////////////////////

bool
Detpoint_file_header_t::is_compatible_version() const
{
    if(MAJOR_VERSION_FROM_UINT16(program_version) <= MAJOR_PROGRAM_VERSION)
	return true;

    if(MINOR_VERSION_FROM_UINT16(program_version) <= MINOR_PROGRAM_VERSION)
	return true;
    else
	return false;
}

//////////////////////////////////////////////////////////////////////

Error_t::Error_t()
{
    min_abs_error = 0.0;
    max_abs_error = 0.0;
    mean_error = 0.0;
    mean_abs_error = 0.0;
}

//////////////////////////////////////////////////////////////////////

void
Error_t::read_from_file(FILE * in)
{
    min_abs_error = read_double(in);
    max_abs_error = read_double(in);
    mean_abs_error = read_double(in);
    mean_error = read_double(in);
}

//////////////////////////////////////////////////////////////////////

void
Error_t::write_to_file(FILE * out) const
{
    write_double(out, min_abs_error);
    write_double(out, max_abs_error);
    write_double(out, mean_abs_error);
    write_double(out, mean_error);
}

//////////////////////////////////////////////////////////////////////

bool
Error_t::is_valid() const
{
    return (min_abs_error >= 0.0) && 
	(mean_abs_error >= 0.0) &&
	(min_abs_error <= max_abs_error);
}

//////////////////////////////////////////////////////////////////////

Detpoint_chunk_header_t::Detpoint_chunk_header_t() 
{
    chunk_mark[0] = 'C';
    chunk_mark[1] = 'N';
    chunk_mark[2] = 'K';
    chunk_mark[3] = 0;

    number_of_bytes = 0;
    number_of_samples = 0;

    chunk_type = 0;
}

//////////////////////////////////////////////////////////////////////

void
Detpoint_chunk_header_t::read_from_file(FILE * in)
{
    chunk_mark[0] = read_uint8(in);
    chunk_mark[1] = read_uint8(in);
    chunk_mark[2] = read_uint8(in);
    chunk_mark[3] = read_uint8(in);

    number_of_bytes = read_uint64(in);
    number_of_samples = read_uint32(in);

    chunk_type = read_uint32(in);

    compression_error.read_from_file(in);
}

//////////////////////////////////////////////////////////////////////

void
Detpoint_chunk_header_t::write_to_file(FILE * out) const
{
    write_uint8(out, chunk_mark[0]);
    write_uint8(out, chunk_mark[1]);
    write_uint8(out, chunk_mark[2]);
    write_uint8(out, chunk_mark[3]);

    write_uint64(out, number_of_bytes);
    write_uint32(out, number_of_samples);

    write_uint32(out, chunk_type);

    compression_error.write_to_file(out);
}

//////////////////////////////////////////////////////////////////////

bool
Detpoint_chunk_header_t::is_valid() const
{
    if(chunk_mark[0] != 'C' ||
       chunk_mark[1] != 'N' ||
       chunk_mark[2] != 'K' ||
       chunk_mark[3] != 0 ||
       number_of_bytes == 0 ||
       number_of_samples == 0 ||
       chunk_type < 10 || chunk_type > 14)
	return false;

    return true;
}
