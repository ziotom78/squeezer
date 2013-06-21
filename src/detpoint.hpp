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

#ifndef DETPOINT_HPP
#define DETPOINT_HPP

#include <vector>
#include <stdexcept>

#include <fitsio.h>

#include "config.hpp"
#include "common_defs.hpp"
#include "data_container.hpp"

#if HAVE_TOODI
#include <LowLevelIO.h>
#endif

struct Detector_pointings_t : public Data_container_t {
    std::vector<double> theta;
    std::vector<double> phi;
    std::vector<double> psi;

    virtual size_t number_of_columns() const { return 5; }

#if HAVE_TOODI
    virtual void read_from_database(const std::string & obj_name);
#endif

    virtual void read_from_fits_file(const std::string & file_name);
    virtual void write_to_fits_file(fitsfile * fptr, int & status);
};

#endif
