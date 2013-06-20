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

#include <fitsio.h>
#include <cstring>

#include "common_defs.hpp"
#include "datadiff.hpp"
#include "toodi_obj.hpp"

//////////////////////////////////////////////////////////////////////

#if HAVE_TOODI
void
Differenced_data_t::read_from_database(const std::string & obj_name)
{
    SessionHandle session;
    ObjectHandle init_handle = 0;
    toodiOpenDiffDbStores(init_handle, &session, "");
    toodiBeginTransaction(session);

    {
        ToodiObject detpoint_obj(session, "toi.science.LFI_DataDiff", obj_name);

	detpoint_obj.read_column_of_double("sampleOBT", obt_times);
	detpoint_obj.read_column_of_double("sampleSCET", scet_times);
	detpoint_obj.read_column_of_double("sky_load", sky_load);
	detpoint_obj.read_column_of_uint32("qualityFlag", quality_flags);
    }

    toodiCommitTransaction(session);
    toodiCloseStores(session);
}
#endif

//////////////////////////////////////////////////////////////////////

static int
read_double_vector_from_fits(fitsfile * fptr,
			     LONGLONG num_of_rows,
			     const std::string & column_name,
			     std::vector<double> & vector)
{
    vector.resize(num_of_rows);

    int status = 0;
    int column_number = 0;
    char * column_name_asciiz = strdup(column_name.c_str());
    fits_get_colnum(fptr, CASEINSEN, column_name_asciiz, 
		    &column_number, &status);
    free(column_name_asciiz);

    fits_read_col(fptr, TDOUBLE, column_number,
		  1, 1, num_of_rows, NULL,
		  vector.data(), NULL, &status);
    return status;
}

//////////////////////////////////////////////////////////////////////

static int
read_uint32_vector_from_fits(fitsfile * fptr,
			     LONGLONG num_of_rows,
			     const std::string & column_name,
			     std::vector<uint32_t> & vector)
{
    int status = 0;
    int column_number = 0;
    char * column_name_asciiz = strdup(column_name.c_str());
    fits_get_colnum(fptr, CASEINSEN, column_name_asciiz, 
		    &column_number, &status);
    free(column_name_asciiz);

    std::vector<unsigned long> long_vector;
    long_vector.resize(num_of_rows);
    fits_read_col_ulng(fptr, column_number,
		       1, 1, num_of_rows, 0UL,
		       long_vector.data(), NULL, &status);

    if(status == 0) {
	vector.resize(num_of_rows);
	std::copy(long_vector.begin(), long_vector.end(),
		  vector.begin());
    }

    return status;
}

//////////////////////////////////////////////////////////////////////

void
Differenced_data_t::read_from_fits_file(const std::string & file_name)
{
    fitsfile * fptr;
    int status = 0;
    fits_open_table(&fptr, const_cast<char *>(file_name.c_str()),
		    READONLY, &status);
    if(status != 0)
	goto throw_error;

    LONGLONG num_of_rows;
    fits_get_num_rowsll(fptr, &num_of_rows, &status);
    if(status != 0)
	goto close_and_throw_error;

    if(read_double_vector_from_fits(fptr, num_of_rows, "OBT", obt_times) != 0 ||
       read_double_vector_from_fits(fptr, num_of_rows, "SCET", scet_times) != 0 ||
       read_double_vector_from_fits(fptr, num_of_rows, "skyLoad", sky_load) != 0 ||
       read_uint32_vector_from_fits(fptr, num_of_rows, "flag", quality_flags) != 0)
	goto close_and_throw_error;

    fits_close_file(fptr, &status);
    if(status != 0)
	goto throw_error;

    return;

close_and_throw_error:
    status = 0; // Reset the error status
    fits_close_file(fptr, &status);
    
throw_error:
    char error_msg[80];
    std::string error_string;
    while(fits_read_errmsg(error_msg) != 0) {
	error_string += error_msg;
	error_string += '\n';
    }
    throw std::runtime_error(error_string);
}

//////////////////////////////////////////////////////////////////////

void
Differenced_data_t::write_to_fits_file(fitsfile * fptr,
				       const Radiometer_t & radiometer,
				       uint16_t od,
				       int & status)
{
    // Since we're going to use a few gotos, it is better to declare
    // all the variables before the first goto
    char * ttype[] = { "OBT", "SCET", "skyLoad", "flags" };
    char * tform[] = { "1D", "1D", "1D", "1J" };
    char * tunit[] = { "Clock ticks", "ms", "V", "dimensionless" };

    double firstobt = obt_times.front();
    double lastobt = obt_times.back();
    double firstsct = scet_times.front();
    double lastsct = scet_times.back();

    char extname[30];

    std::strncpy(extname, radiometer.to_str().c_str(), sizeof(extname));
    if(fits_create_tbl(fptr, BINARY_TBL, obt_times.size(), 4, 
		       ttype, tform, tunit, extname, &status) != 0)
	return;

    std::vector<unsigned long> ulong_flags;
    ulong_flags.resize(quality_flags.size());
    std::copy(quality_flags.begin(), quality_flags.end(),
	      ulong_flags.begin());

    if(fits_write_col(fptr, TDOUBLE, 1, 1, 1, 
		      obt_times.size(), obt_times.data(), &status) != 0 ||
       fits_write_col(fptr, TDOUBLE, 2, 1, 1, 
		      scet_times.size(), scet_times.data(), &status) != 0 ||
       fits_write_col(fptr, TDOUBLE, 3, 1, 1, 
		      sky_load.size(), sky_load.data(), &status) != 0 ||
       fits_write_col(fptr, TULONG, 4, 1, 1, 
		      ulong_flags.size(), ulong_flags.data(), &status) != 0)
	return;

    if(fits_write_key(fptr, TDOUBLE, "FIRSTOBT", 
		      (void *) &firstobt, "First OBT time", &status) != 0 ||
       fits_write_key(fptr, TDOUBLE, "LASTOBT", 
		      (void *) &lastobt, "Last OBT time", &status) != 0 ||
       fits_write_key(fptr, TDOUBLE, "FIRSTSCT", 
		      (void *) &firstsct, "First SCET time [ms]", &status) != 0 ||
       fits_write_key(fptr, TDOUBLE, "LASTSCT", 
		      (void *) &lastsct, "Last SCET time [ms]", &status) != 0 ||
       fits_write_key(fptr, TSHORT, "OD",
		      (void *) &od, "Operational day", &status) != 0 ||
       fits_write_key(fptr, TBYTE, "HORN",
		      (void *) &radiometer.horn, "Horn number (18..28)", &status) != 0 ||
       fits_write_key(fptr, TBYTE, "RAD",
		      (void *) &radiometer.arm, "Radiometer number (0, 1)", &status) != 0 ||
       fits_write_date(fptr, &status) != 0)
	return;
}
