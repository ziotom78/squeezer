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

#include <cassert>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <memory>

#include <gsl/gsl_math.h>

#include "common_defs.hpp"
#include "data_structures.hpp"
#include "byte_buffer.hpp"
#include "run_length_encoding.hpp"
#include "poly_fit_encoding.hpp"
#include "datadiff.hpp"
#include "detpoint.hpp"
#include "decompress.hpp"

//////////////////////////////////////////////////////////////////////

void
decompress_obt_times(Byte_buffer_t & buffer,
		     double first_obt,
		     size_t num_of_samples,
		     std::vector<double> & dest)
{
    std::vector<uint32_t> obt_delta_values;
    rle_decompression(buffer, num_of_samples, obt_delta_values);

    dest.resize(obt_delta_values.size() + 1);
    dest[0] = first_obt;

    for(size_t idx = 1; idx < dest.size(); ++idx) {
	dest[idx] = dest[idx - 1] + obt_delta_values[idx - 1];
    }
}

//////////////////////////////////////////////////////////////////////

void
decompress_scet_times(Byte_buffer_t & buffer,
		      const Squeezer_file_header_t & file_header,
		      const std::vector<double> & obt_times,
		      std::vector<double> & dest)
{
    dest.resize(obt_times.size());

    const double slope = 
	(file_header.last_scet_in_ms - file_header.first_scet_in_ms) / 
	(file_header.last_obt - file_header.first_obt);

    for(size_t idx = 0; idx < obt_times.size(); ++idx) {
	double interpolated_scet =
	    file_header.first_scet_in_ms + slope * (obt_times[idx] - file_header.first_obt);
	double scet_correction = buffer.read_float();

	dest[idx] = interpolated_scet + scet_correction;
    }
}

//////////////////////////////////////////////////////////////////////

void
decompress_angles(Byte_buffer_t & buffer,
		  size_t num_of_samples,
		  std::vector<double> & dest,
		  const Decompression_parameters_t & params)
{
    dest.resize(num_of_samples);

    poly_fit_decode(num_of_samples, buffer, dest);

    // Clip angles within [0, 2pi]
    double offset = 0.0;
    for(size_t idx = 0; idx < dest.size(); ++idx) {
	if(dest[idx] + offset < 0.0)
	    offset += 2 * M_PI;
	else if(dest[idx] + offset >= 2 * M_PI)
	    offset -= 2 * M_PI;

	dest[idx] += offset;
    }
}

//////////////////////////////////////////////////////////////////////

void
decompress_scientific_data(Byte_buffer_t & buffer,
			   size_t num_of_samples,
			   std::vector<double> & dest)
{
    dest.resize(num_of_samples);
    for(size_t idx = 0; idx < num_of_samples; ++idx)
	dest[idx] = buffer.read_float();
}

//////////////////////////////////////////////////////////////////////

void
decompress_quality_flags(Byte_buffer_t & buffer,
			 size_t num_of_samples,
			 std::vector<uint32_t> & dest)
{
    rle_decompression(buffer, num_of_samples, dest);
}

//////////////////////////////////////////////////////////////////////

void
decompress_chunk(size_t chunk_idx,
		 const Squeezer_file_header_t & file_header,
		 const Squeezer_chunk_header_t & chunk_header,
		 FILE * input_file,
		 const Decompression_parameters_t & params,
		 Data_container_t * data_container)
{
    if(! chunk_header.is_valid()) {

	std::cerr << PROGRAM_NAME
		  << ": the input file seems to have been corrupted.\n";
	return;

    } else if(params.verbose_flag) {

	std::cerr << PROGRAM_NAME
		  << ": reading data chunk #"
		  << chunk_idx + 1
		  << " (";

	switch(chunk_header.chunk_type) {
	case CHUNK_DELTA_OBT: std::cerr << "OBT times"; break;
	case CHUNK_SCET_ERROR: std::cerr << "SCET times"; break;
	case CHUNK_THETA: std::cerr << "theta angle"; break;
	case CHUNK_PHI: std::cerr << "phi angle"; break;
	case CHUNK_PSI: std::cerr << "psi angle"; break;
	default: std::cerr << "unknown chunk";
	}

	std::cerr << ")\n";

    }

    Byte_buffer_t chunk_data;
    chunk_data.buffer.resize(chunk_header.number_of_bytes);
    if(fread(chunk_data.buffer.data(), 
	     1, 
	     chunk_header.number_of_bytes, 
	     input_file) < chunk_header.number_of_bytes) {

	std::cerr << PROGRAM_NAME
		  << ": unable to read the contents of chunk #"
		  << chunk_idx + 1
		  << ", perhaps the file is corrupted or disappeared "
		  << "during reading\n";
	return;

    }

    switch(chunk_header.chunk_type) {
    case CHUNK_DELTA_OBT:
	decompress_obt_times(chunk_data, 
			     file_header.first_obt,
			     chunk_header.number_of_samples,
			     data_container->obt_times);
	break;
    case CHUNK_SCET_ERROR:
	if(data_container->obt_times.empty()) {
	    std::cerr << PROGRAM_NAME
		      << ": malformed chunk #"
		      << chunk_idx + 1
		      << ", SCET times have been found here but no "
		      << "OBT times have been read yet\n";

	    return;
	}

	decompress_scet_times(chunk_data, 
			      file_header,
			      data_container->obt_times, 
			      data_container->scet_times);
	break;
    case CHUNK_THETA:
    {
	Detector_pointings_t * detpoints =
	    dynamic_cast<Detector_pointings_t *>(data_container);
	decompress_angles(chunk_data, 
			  chunk_header.number_of_samples,
			  detpoints->theta,
			  params);
	break;
    }
    case CHUNK_PHI:
    {
	Detector_pointings_t * detpoints =
	    dynamic_cast<Detector_pointings_t *>(data_container);
	decompress_angles(chunk_data, 
			  chunk_header.number_of_samples,
			  detpoints->phi,
			  params);
	break;
    }
    case CHUNK_PSI:
    {
	Detector_pointings_t * detpoints =
	    dynamic_cast<Detector_pointings_t *>(data_container);
	decompress_angles(chunk_data, 
			  chunk_header.number_of_samples,
			  detpoints->psi,
			  params);
	break;
    }
    case CHUNK_DIFFERENCED_DATA:
    {
	Differenced_data_t * datadiff =
	    dynamic_cast<Differenced_data_t *>(data_container);
	decompress_scientific_data(chunk_data, 
				   chunk_header.number_of_samples,
				   datadiff->sky_load);
	break;
    }
    case CHUNK_QUALITY_FLAGS:
    {
	Differenced_data_t * datadiff =
	    dynamic_cast<Differenced_data_t *>(data_container);
	decompress_quality_flags(chunk_data, 
				 chunk_header.number_of_samples,
				 datadiff->quality_flags);
	break;
    }
    default:
	abort();
    }
}

//////////////////////////////////////////////////////////////////////

Data_container_t *
decompress_from_file(FILE * input_file,
		     const Decompression_parameters_t & params)
{
    Squeezer_file_header_t file_header(SQZ_NO_DATA);
    file_header.read_from_file(input_file);
    if(! file_header.is_valid()) {
	std::cerr << PROGRAM_NAME
		  << ": the input file does not seem to have been "
		  << "created by \"squeezer\". It might have been damaged.\n";
	return nullptr;
    }

    if(! file_header.is_compatible_version()) {
	std::cerr << PROGRAM_NAME
		  << ": the input file seems to have been created by "
		  << "a version of \"squeezer\" ("
		  << MAJOR_VERSION_FROM_UINT16(file_header.program_version)
		  << '.'
		  << MINOR_VERSION_FROM_UINT16(file_header.program_version)
		  << ") than this executable ("
		  << MAJOR_PROGRAM_VERSION
		  << '.'
		  << MINOR_PROGRAM_VERSION
		  << "). I will cowardly stop here.";

	return nullptr;
    }

    Squeezer_file_type_t file_type = file_header.get_type();
    Data_container_t * file_data;
    switch(file_type) {
    case SQZ_DETECTOR_POINTINGS:
	file_data = new Detector_pointings_t();
	break;

    case SQZ_DIFFERENCED_DATA:
	file_data = new Differenced_data_t(false);
	break;

    default:
	abort();
    }

    file_data->radiometer = file_header.radiometer;
    file_data->od = file_header.od;

    if(params.verbose_flag) {
	std::string data_type;
	switch(file_type) {
	case SQZ_DETECTOR_POINTINGS: data_type = "detector pointings"; break;
	case SQZ_DIFFERENCED_DATA: data_type = "differenced data"; break;
	default:
	    data_type = "unknown data"; break;
	}

	std::cerr << PROGRAM_NAME << ": the file contains "
		  << data_type
		  << " for radiometer "
		  << file_header.radiometer.to_str()
		  << ", OD "
		  << file_header.od
		  << '\n';
    }

    for(size_t idx = 0; idx < file_header.number_of_chunks; ++idx) {

	Squeezer_chunk_header_t chunk_header;
	chunk_header.read_from_file(input_file);

	decompress_chunk(idx, 
			 file_header,
			 chunk_header, 
			 input_file, 
			 params, 
			 file_data);

    }

    return file_data;
}

//////////////////////////////////////////////////////////////////////

void
decompress_file_from_file(FILE * input_file,
			  const std::string & output_file_name,
			  const Decompression_parameters_t & params)
{
    std::unique_ptr<Data_container_t> file_data(decompress_from_file(input_file, params));

    if(params.verbose_flag) {
	std::cerr << PROGRAM_NAME
		  << ": writing detector pointings to file "
		  << output_file_name
		  << '\n';
    }

    fitsfile * fptr;
    int status = 0;

    fits_create_file(&fptr,
		     const_cast<char *>(output_file_name.c_str()),
		     &status);

    if(status == 0) {

	file_data->write_to_fits_file(fptr, status);
	fits_close_file(fptr, &status);
    
    }

    if(status != 0) {

	char error_msg[80];
	std::string error_string;
	while(fits_read_errmsg(error_msg) != 0) {
	    error_string += error_msg;
	    error_string += '\n';
	}
	throw std::runtime_error(error_string);

    }
}
