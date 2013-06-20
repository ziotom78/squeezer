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

#include <iostream>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstring>

#include <gsl/gsl_math.h>

#include "common_defs.hpp"
#include "file_io.hpp"
#include "run_length_encoding.hpp"
#include "poly_fit_encoding.hpp"
#include "compress.hpp"
#include "detpoint.hpp"
#include "data_structures.hpp"

const uint32_t OPTIMAL_CHUNK_SIZE = 128 * 1024 * 1024; // 128 MB

//////////////////////////////////////////////////////////////////////

double
rad2arcmin(double x)
{
    return x * 180.0 / M_PI * 3600.0;
}

//////////////////////////////////////////////////////////////////////

void
initialize_file_header(Detpoint_file_header_t & file_header,
		       const Detector_pointings_t & detpoints,
		       const Compression_parameters_t & params)
{
    file_header.radiometer = params.radiometer;
    file_header.od = params.od_number;

    file_header.first_obt = detpoints.obt_times.front();
    file_header.last_obt = detpoints.obt_times.back();

    file_header.first_scet_in_ms = detpoints.scet_times.front();
    file_header.last_scet_in_ms = detpoints.scet_times.back();
    file_header.number_of_chunks = 5;
}

//////////////////////////////////////////////////////////////////////

void
compress_obt(const std::vector<double> & obt,
	     FILE * output_file,
	     const Compression_parameters_t & params)
{
    std::vector<uint32_t> obt_delta;
    obt_delta.resize(obt.size() - 1);
    for(size_t idx = 0; idx < obt_delta.size(); ++idx) {
	obt_delta[idx] = obt[idx + 1] - obt[idx];
    }

    Byte_buffer_t obt_delta_buffer;
    rle_compression(obt_delta.data(),
		    obt_delta.size(),
		    obt_delta_buffer);

    Detpoint_chunk_header_t chunk_header;
    chunk_header.number_of_bytes = obt_delta_buffer.buffer.size();
    chunk_header.number_of_samples = obt_delta.size();
    chunk_header.chunk_type = CHUNK_DELTA_OBT;

    chunk_header.compression_error.min_abs_error = 0.0;
    chunk_header.compression_error.max_abs_error = 0.0;
    chunk_header.compression_error.mean_abs_error = 0.0;
    chunk_header.compression_error.mean_error = 0.0;

    chunk_header.write_to_file(output_file);
    obt_delta_buffer.write_to_file(output_file);

    if(params.verbose_flag) {
	std::cerr << PROGRAM_NAME
		  << ": the size of the OBT times shrunk from "
		  << obt.size() * sizeof(obt[0])
		  << " to "
		  << obt_delta_buffer.buffer.size()
		  << " bytes (using run-length encoding)\n";
	
	std::cerr << PROGRAM_NAME
		  << ":     the overall compression factor for OBT times is "
		  << (obt.size() * sizeof(obt[0])) * (1.0 / obt_delta_buffer.buffer.size())
		  << '\n';
    }

}

//////////////////////////////////////////////////////////////////////

void
estimate_scet_reconstruction_error(const std::vector<double> & scet, 
				   const std::vector<double> & obt, 
				   const std::vector<float> & scet_interp_error,
				   Error_t & compression_error)
{
    const double slope = (scet.back() - scet.front()) / (obt.back() - obt.front());

    compression_error.mean_abs_error = 0.0;
    compression_error.mean_error = 0.0;

    for(size_t idx = 0; idx < scet_interp_error.size(); ++idx) {
	double reconstructed_scet =
	    scet.front() + slope * (obt[idx] - obt.front())
	    + scet_interp_error[idx];

	double error = reconstructed_scet - scet[idx];
	double abs_error = fabs(error);

	compression_error.mean_abs_error += abs_error;
	compression_error.mean_error += error;

	if(idx == 0) {
	    compression_error.min_abs_error = abs_error;
	    compression_error.max_abs_error = abs_error;
	} else {
	    if(abs_error < compression_error.min_abs_error)
		compression_error.min_abs_error = abs_error;

	    if(abs_error > compression_error.max_abs_error)
		compression_error.max_abs_error = abs_error;
	}
    }

    compression_error.mean_abs_error /= scet_interp_error.size();
    compression_error.mean_error /= scet_interp_error.size();
}

//////////////////////////////////////////////////////////////////////

void
compress_scet(const std::vector<double> & scet,
	      const std::vector<double> & obt,
	      FILE * output_file,
	      const Compression_parameters_t & params)
{
    std::vector<float> scet_interp_error;
    scet_interp_error.resize(scet.size());

    const double slope = (scet.back() - scet.front()) / (obt.back() - obt.front());
    for(size_t idx = 0; idx < scet.size(); ++idx) {

	double interpolated_scet =
	    scet.front() + slope * (obt[idx] - obt.front());
	scet_interp_error[idx] = scet[idx] - interpolated_scet;

    }

    Byte_buffer_t buffer;
    for(auto value : scet_interp_error) {
	buffer.append_float(value);
    }

    Detpoint_chunk_header_t chunk_header;
    chunk_header.number_of_bytes = buffer.buffer.size();
    chunk_header.number_of_samples = scet_interp_error.size();
    chunk_header.chunk_type = CHUNK_SCET_ERROR;

    estimate_scet_reconstruction_error(scet, 
				       obt, 
				       scet_interp_error,
				       chunk_header.compression_error);

    chunk_header.write_to_file(output_file);
    buffer.write_to_file(output_file);

    if(params.verbose_flag) {
	std::cout << PROGRAM_NAME
		  << ": the size of the SCET times shrunk from "
		  << scet.size() * sizeof(scet[0])
		  << " to "
		  << buffer.buffer.size()
		  << " bytes\n";

	std::cout << PROGRAM_NAME
		  << ":     the overall compression factor is "
		  << (scet.size() * sizeof(scet[0])) * (1.0 / buffer.buffer.size())
		  << '\n';

	std::cout << PROGRAM_NAME
		  << ":     the maximum error is "
		  << chunk_header.compression_error.max_abs_error
		  << " ms\n";

	std::cerr << PROGRAM_NAME
		  << ":     the average absolute error is "
		  << chunk_header.compression_error.mean_abs_error
		  << " ms\n";
    }
}

//////////////////////////////////////////////////////////////////////

void
estimate_angle_reconstruction_error(const std::vector<double> angle,
				    const std::vector<double> reconstructed_angle,
				    Error_t & compression_error)
{
    compression_error.mean_abs_error = 0.0;
    compression_error.mean_error = 0.0;

    for(size_t idx = 0; idx < reconstructed_angle.size(); ++idx) {

	double error = angle[idx] - reconstructed_angle[idx];
	if(error > M_PI)
	    error -= M_PI * 2;
	else if(error < -M_PI)
	    error += M_PI * 2;

	double abs_error = std::fabs(error);

	compression_error.mean_abs_error += abs_error;
	compression_error.mean_error += error;

	if(idx == 0) {
	    compression_error.min_abs_error = abs_error;
	    compression_error.max_abs_error = abs_error;
	} else {
	    if(compression_error.min_abs_error > abs_error)
		compression_error.min_abs_error = abs_error;
	    if(compression_error.max_abs_error < abs_error)
		compression_error.max_abs_error = abs_error;
	}
    }

    compression_error.mean_abs_error /= angle.size();
    compression_error.mean_error /= angle.size();
}

//////////////////////////////////////////////////////////////////////

void
compress_angle(const std::vector<double> & angle,
	       Chunk_type_t chunk_type,
	       FILE * output_file,
	       const Compression_parameters_t & params)
{
    Byte_buffer_t output_buffer;
    size_t num_of_frames = 0;
    size_t num_of_frames_encoded_directly = 0;
    poly_fit_encode(angle, 
		    params.elements_per_frame,
		    params.number_of_poly_terms,
		    params.max_abs_error,
		    output_buffer,
		    num_of_frames,
		    num_of_frames_encoded_directly);
    
    Detpoint_chunk_header_t chunk_header;
    chunk_header.number_of_bytes = output_buffer.size();
    chunk_header.number_of_samples = angle.size();
    chunk_header.chunk_type = chunk_type;

    std::vector<double> reconstructed_angle;
    poly_fit_decode(angle.size(),
		    output_buffer,
		    reconstructed_angle);

    estimate_angle_reconstruction_error(angle, 
					reconstructed_angle,
					chunk_header.compression_error);

    chunk_header.write_to_file(output_file);
    output_buffer.write_to_file(output_file);

    if(params.verbose_flag) {
	std::cerr << PROGRAM_NAME
		  << ": the size of the angle vector shrunk from "
		  << angle.size() * sizeof(angle[0])
		  << " to "
		  << output_buffer.size()
		  << " bytes (using polynomial encoding)\n";

	std::cerr << PROGRAM_NAME
		  << ":     "
		  << num_of_frames
		  << " frames written, of which "
		  << num_of_frames_encoded_directly
		  << " were uncompressed ("
		  << (num_of_frames_encoded_directly * 100) / num_of_frames
		  << "%)\n";

	std::cerr << PROGRAM_NAME
		  << ":     the overall compression factor is "
		  << (angle.size() * sizeof(angle[0])) * (1.0 / output_buffer.size())
		  << '\n';

	std::cerr << PROGRAM_NAME
		  << ":     the absolute error ranges from "
		  << rad2arcmin(chunk_header.compression_error.min_abs_error)
		  << " to "
		  << rad2arcmin(chunk_header.compression_error.max_abs_error)
		  << " arcsec\n";

	std::cerr << PROGRAM_NAME
		  << ":     the average absolute error is "
		  << rad2arcmin(chunk_header.compression_error.mean_abs_error)
		  << " arcsec\n";
    }
}

//////////////////////////////////////////////////////////////////////

void
compress_detpoints_to_file(const std::string & input_file_name,
			   FILE * output_file,
			   const Compression_parameters_t & params)
{
    Detector_pointings_t detpoints;
    if(params.verbose_flag) {
	std::cerr << PROGRAM_NAME << ": reading pointings from " 
		  << input_file_name << '\n';
    }

    bool detpoints_read = false;

#if HAVE_TOODI
    if(input_file_name.compare(0, 6, "TOODI%") == 0) {
        detpoints.read_from_database(input_file_name);
        detpoints_read = true;
    }
#endif

    if(! detpoints_read) {
        detpoints.read_from_fits_file(input_file_name);
    }

    Detpoint_file_header_t file_header;
    initialize_file_header(file_header, detpoints, params);

    file_header.write_to_file(output_file);

    compress_obt(detpoints.obt_times, output_file, params);
    compress_scet(detpoints.scet_times, 
		  detpoints.obt_times,
		  output_file,
		  params);
    compress_angle(detpoints.theta, CHUNK_THETA, output_file, params);
    compress_angle(detpoints.phi,   CHUNK_PHI, output_file, params);
    compress_angle(detpoints.psi,   CHUNK_PSI, output_file, params);
}
