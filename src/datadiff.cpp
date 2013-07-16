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
#include <algorithm>
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
	std::string obj_type;
	if(calibrated)
	    obj_type = "toi.science.LFI_DataDiffReduced";
	else
	    obj_type = "toi.science.LFI_DataDiff";

        ToodiObject detpoint_obj(session, obj_type, obj_name);

	detpoint_obj.read_column_of_double("sampleOBT", obt_times);
	detpoint_obj.read_column_of_double("sampleSCET", scet_times);
	detpoint_obj.read_column_of_double("skyLoad", sky_load);
	detpoint_obj.read_column_of_int32("qualityFlag", quality_flags);
    }

    toodiCommitTransaction(session);
    toodiCloseStores(session);
}
#endif

//////////////////////////////////////////////////////////////////////

static int
load_data(long total_num,
	  long offset,
	  long first_num,
	  long num_of_values,
	  int num_of_arrays,
	  iteratorCol * data,
	  void * user_ptr)
{
    auto datadiff = reinterpret_cast<Differenced_data_t *>(user_ptr);

    std::copy_n(((double *) data[0].array) + 1,
		num_of_values,
		datadiff->obt_times.data() + first_num - 1);

    std::copy_n(((double *) data[1].array) + 1,
		num_of_values,
		datadiff->scet_times.data() + first_num - 1);

    std::copy_n(((double *) data[2].array) + 1,
		num_of_values,
		datadiff->sky_load.data() + first_num - 1);

    std::copy_n(((unsigned long *) data[3].array) + 1,
		num_of_values,
		datadiff->quality_flags.data() + first_num - 1);

    return 0;
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

    char radiometer_name[FLEN_KEYWORD];
    fits_read_key_str(fptr, "EXTNAME", &radiometer_name[0], NULL, &status);
    if(status != 0)
	goto close_and_throw_error;

    obt_times.resize(num_of_rows);
    scet_times.resize(num_of_rows);
    sky_load.resize(num_of_rows);
    quality_flags.resize(num_of_rows);

    {
	iteratorCol cols[4];
	int n_cols = 4;
	fits_iter_set_by_name(&cols[0], fptr, (char *) "OBT", TDOUBLE, InputCol);
	fits_iter_set_by_name(&cols[1], fptr, (char *) "SCET", TDOUBLE, InputCol);
	fits_iter_set_by_name(&cols[2], fptr, (char *) radiometer_name, TDOUBLE, InputCol);
	fits_iter_set_by_name(&cols[3], fptr, (char *) "flag", TULONG, InputCol);

	fits_iterate_data(n_cols, cols, 0, 0, &load_data, this, &status);
    }
    if(status != 0)
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

static int
save_data(long total_num,
	  long offset,
	  long first_num,
	  long num_of_values,
	  int num_of_arrays,
	  iteratorCol * data,
	  void * user_ptr)
{
    auto datadiff = reinterpret_cast<Differenced_data_t *>(user_ptr);

    ((double *) data[0].array)[0] = DOUBLENULLVALUE;
    ((double *) data[1].array)[0] = DOUBLENULLVALUE;
    ((double *) data[2].array)[0] = DOUBLENULLVALUE;
    ((unsigned long *) data[3].array)[0] = 0;

    std::copy_n(datadiff->obt_times.data() + first_num - 1,
		num_of_values,
		((double *) data[0].array) + 1);

    std::copy_n(datadiff->scet_times.data() + first_num - 1,
		num_of_values,
		((double *) data[1].array) + 1);

    std::copy_n(datadiff->sky_load.data() + first_num - 1,
		num_of_values,
		((double *) data[2].array) + 1);

    std::copy_n(datadiff->quality_flags.data() + first_num - 1,
		num_of_values,
		((unsigned long *) data[3].array) + 1);

    return 0;
}

//////////////////////////////////////////////////////////////////////

void
Differenced_data_t::write_to_fits_file(fitsfile * fptr, int & status)
{
    // Since we're going to use a few gotos, it is better to declare
    // all the variables before the first goto
    std::vector<char *> ttype { (char *) "OBT", 
                                (char *) "SCET", 
                                (char *) "",
                                (char *) "FLAG" };
    char * tform[] = { (char *) "1D", 
		       (char *) "1D", 
		       (char *) "1D", 
		       (char *) "1J" };
    char * tunit[] = { (char *) "Clock ticks", 
		       (char *) "ms", 
		       (char *) "V", 
		       (char *) "dimensionless" };

    double firstobt = obt_times.front();
    double lastobt = obt_times.back();
    double firstsct = scet_times.front();
    double lastsct = scet_times.back();

    char extname[30];

    ttype[2] = (char *) alloca(12);
    std::strncpy(ttype[2], radiometer.to_str().c_str(), 12);
    std::strncpy(extname, radiometer.to_str().c_str(), sizeof(extname));
    if(fits_create_tbl(fptr, BINARY_TBL, obt_times.size(), 4, 
		       ttype.data(), tform, tunit, extname, &status) != 0)
	return;

    {
	iteratorCol cols[4];
	int n_cols = 4;
	fits_iter_set_by_num(&cols[0], fptr, 1, TDOUBLE, OutputCol);
	fits_iter_set_by_num(&cols[1], fptr, 2, TDOUBLE, OutputCol);
	fits_iter_set_by_num(&cols[2], fptr, 3, TDOUBLE, OutputCol);
	fits_iter_set_by_num(&cols[3], fptr, 4, TULONG, OutputCol);

	fits_iterate_data(n_cols, cols, 0, 0, &save_data, this, &status);
    }
    if(status != 0)
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
