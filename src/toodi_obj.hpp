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

#ifndef TOODI_OBJ_HPP
#define TOODI_OBJ_HPP

#include <string>
#include <vector>

#if HAVE_TOODI
#include <LowLevelIO.h>
#endif

//////////////////////////////////////////////////////////////////////

#if HAVE_TOODI
class ToodiObject {
public:
    ObjectHandle object_handle;

    ToodiObject(SessionHandle session_handle,
		const std::string & type,
		const std::string & object_name) {

        toodiOpenPersistentObject(session_handle,
				  type.c_str(),
				  object_name.c_str(),
				  &object_handle);

    }

    virtual ~ToodiObject() {
        toodiCloseObject(object_handle);
    }

    void read_column_of_double(const std::string & column_name,
			       std::vector<double> & vector) {
        toodiInt32 column_idx;
	toodiGetIndexOfColumn(object_handle, column_name.c_str(), &column_idx);

        toodiInt64 num_of_rows;
	toodiGetSizeOfColumn(object_handle, column_idx, &num_of_rows);

	vector.resize(num_of_rows);
	toodiGetDoubleValues(object_handle, column_idx, 0, num_of_rows, vector.data());
    }

    void read_column_of_int32(const std::string & column_name,
			       std::vector<uint32_t> & vector) {
        toodiInt32 column_idx;
	toodiGetIndexOfColumn(object_handle, column_name.c_str(), &column_idx);

        toodiInt64 num_of_rows;
	toodiGetSizeOfColumn(object_handle, column_idx, &num_of_rows);

	vector.resize(num_of_rows);
	toodiGetInt32Values(object_handle, column_idx, 0, num_of_rows, 
			    reinterpret_cast<toodiInt32 *>(vector.data()));
    }
};

#endif

#endif
