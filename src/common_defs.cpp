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

#include <sstream>
#include <stdexcept>
#include "common_defs.hpp"

#include <iostream>
void Radiometer_t::parse_from_name(const std::string & name)
{
    std::string horn_str;
    char arm_char;

    if(name.size() == 6) {

	// Long name, e.g., LFI18M
	horn_str = name.substr(3, 2);
	arm_char = name[5];

    } else if(name.size() == 3) {

	// Short name, e.g. 18M
	horn_str = name.substr(0, 2);
	arm_char = name[2];

    } else goto invalid_name;

    {
	std::stringstream ss(horn_str);
	int horn_int;
	ss >> horn_int;
	if(ss.fail())
	    goto invalid_name;

	horn = horn_int;
    }

    if(arm_char == 'M')
	arm = 0;
    else
	arm = 1;

    return;

invalid_name:
    throw std::runtime_error("invalid radiometer name \"" + name + "\"");
}

std::string Radiometer_t::to_str() const 
{
    std::stringstream ss;
    ss << "LFI" << ((int) horn);

    if(arm == 0)
	ss << "M";
    else
	ss << "S";

    return ss.str();
}
