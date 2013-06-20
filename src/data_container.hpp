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

#ifndef DATA_CONTAINER_HPP
#define DATA_CONTAINER_HPP

#include <cstdint>
#include <string>
#include <fitsio.h>

#include "config.hpp"
#include "common_defs.hpp"

#if HAVE_TOODI
#include <LowLevelIO.h>
#endif

struct Data_container_t {
    std::vector<double> obt_times;
    std::vector<double> scet_times;

    Data_container_t() {}
    virtual ~Data_container_t() {}

    virtual double first_obt() const  { return obt_times.front(); }
    virtual double last_obt() const   { return obt_times.back(); }
    virtual double first_scet() const { return scet_times.front(); }
    virtual double last_scet() const  { return scet_times.back(); }

    virtual size_t number_of_columns() const = 0;

#if HAVE_TOODI
    virtual void read_from_database(const std::string & obj_name) = 0;
#endif

    virtual void read_from_fits_file(const std::string & file_name) = 0;
    virtual void write_to_fits_file(fitsfile * fptr,
				    const Radiometer_t & radiometer,
				    uint16_t od,
				    int & status) = 0;
};

#endif
