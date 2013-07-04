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

#include "config.hpp"

#include <assert.h>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <cstring>
#include <cmath>
#include <cstdlib>
#include <stdexcept>

#include <fitsio.h>

#include "data_structures.hpp"
#include "compress.hpp"
#include "decompress.hpp"
#include "detpoint.hpp"
#include "common_defs.hpp"
#include "help.hpp"

//////////////////////////////////////////////////////////////////////

void
run_compression_task_for_one_file(const std::string radiometer_str,
				  const std::string od_str,
				  const std::string input_file_name,
				  const std::string output_file_name,
				  Compression_parameters_t & params)
{
    FILE * output_file = NULL;
    bool write_to_stdout = false;
    if(output_file_name == "-") {
        write_to_stdout = true;
	output_file = stdout;
    } else {
        output_file = std::fopen(output_file_name.c_str(), "wb");
    }

    Radiometer_t radiometer;
    radiometer.parse_from_name(radiometer_str);
    params.radiometer = radiometer;

    std::stringstream ss(od_str);
    ss >> params.od_number;

    compress_file_to_file(input_file_name,
			  output_file,
			  params);

    if(! write_to_stdout) {
        std::fclose(output_file);
    }
}

//////////////////////////////////////////////////////////////////////

void
compress_using_a_parameter_file(const std::string & file_name,
				Compression_parameters_t & params)
{
    std::ifstream input_stream(file_name);

    if(params.verbose_flag) {
      std::cerr << PROGRAM_NAME
		<< ": reading the list of files from "
		<< file_name
		<< '\n';
    }

    size_t num_of_processed_files = 0;
    while(input_stream.good()) {
        std::string cur_line;
	getline(input_stream, cur_line);

	if(cur_line.empty() || cur_line.at(0) == '#')
	  continue;

	std::string radiometer_str;
	std::string od_str;
	std::string input_file_name;
	std::string output_file_name;

	std::stringstream ss (cur_line);
	ss >> radiometer_str;
	ss >> od_str;
	ss >> input_file_name;
	ss >> output_file_name;

	run_compression_task_for_one_file(radiometer_str,
					  od_str,
					  input_file_name,
					  output_file_name,
					  params);
	num_of_processed_files++;
    }

    if(params.verbose_flag) {
      std::cerr << PROGRAM_NAME
		<< ": "
		<< num_of_processed_files
		<< " objects specified in file "
		<< file_name
		<< " have been processed.\n";
    }
}

//////////////////////////////////////////////////////////////////////

void
run_compression_task(const std::vector<std::string> & list_of_arguments)
{
    Compression_parameters_t params;
    size_t cur_argument = 0;

    params.file_type = SQZ_DETECTOR_POINTINGS;

    while((! list_of_arguments.empty()) &&
	  list_of_arguments.at(cur_argument)[0] == '-') {

	if(list_of_arguments.at(cur_argument) == "-v") {

	    params.verbose_flag = true;
	    cur_argument++;

	} else if(list_of_arguments.at(cur_argument) == "--pointings") {

	    params.file_type = SQZ_DETECTOR_POINTINGS;
	    cur_argument++;

	} else if(list_of_arguments.at(cur_argument) == "--datadiff") {

	    params.file_type = SQZ_DIFFERENCED_DATA;
	    cur_argument++;

	} else if(list_of_arguments.at(cur_argument) == "--calibrated") {

	    params.read_calibrated_data = true;
	    cur_argument++;

	} else if(list_of_arguments.at(cur_argument) == "--uncalibrated") {

	    params.read_calibrated_data = false;
	    cur_argument++;

	} else if(list_of_arguments.at(cur_argument) == "-n") {

	    std::stringstream ss(list_of_arguments.at(++cur_argument));
	    size_t number;
	    ss >> number;
	    if(number > UINT8_MAX) {
		std::cerr << PROGRAM_NAME
			  << ": the maximum value allowed for the -n parameter is "
			  << UINT8_MAX
			  << " (you provided "
			  << number
			  << ")\n";
	    } else {
		params.elements_per_frame = number;
	    }

	    ++cur_argument;

	} else if(list_of_arguments.at(cur_argument) == "-p") {

	    std::stringstream ss(list_of_arguments.at(++cur_argument));
	    unsigned number;
	    ss >> number;
	    if(number > UINT8_MAX) {
		std::cerr << PROGRAM_NAME
			  << ": the maximum value allowed for the -p parameter is "
			  << UINT8_MAX
			  << " (you provided "
			  << number
			  << ")\n";
	    } else {
		params.number_of_poly_terms = number;
	    }

	    ++cur_argument;

	} else if(list_of_arguments.at(cur_argument) == "-s") {

	    std::stringstream ss(list_of_arguments.at(++cur_argument));
	    double number;
	    ss >> number;
	    if(number < 0.0) {
		std::cerr << PROGRAM_NAME
			  << ": the value passed to -s is negative. "
			  << "This will disable compression for angles.\n";
	    }

	    params.max_abs_error = number / 3600.0 * M_PI / 180.0;

	    ++cur_argument;

	} else {

	    std::cerr << PROGRAM_NAME
		      << ": unknown flag \"" 
		      << list_of_arguments.at(cur_argument) 
		      << "\"\n";
	    std::exit(1);

	}
    }

    if(params.number_of_poly_terms >= params.elements_per_frame) {
	std::cerr << PROGRAM_NAME
		  << ": invalid numbers specified using -n ("
		  << params.elements_per_frame
		  << ") and -p ("
		  << params.number_of_poly_terms
		  << ")\n";
    }

    if(list_of_arguments.size() - cur_argument != 4 &&
       list_of_arguments.size() - cur_argument != 1) {

	std::cerr << PROGRAM_NAME
		  << ": wrong number of arguments. Run \"squeezer compress help\".\n";
	std::exit(1);

    }

    if(list_of_arguments.size() - cur_argument == 4) {
        run_compression_task_for_one_file(list_of_arguments.at(cur_argument),
					  list_of_arguments.at(cur_argument + 1),
					  list_of_arguments.at(cur_argument + 2),
					  list_of_arguments.at(cur_argument + 3),
					  params);
    } else {
      compress_using_a_parameter_file(list_of_arguments.at(cur_argument),
				      params);
    }
}

//////////////////////////////////////////////////////////////////////

void
run_decompression_task(const std::vector<std::string> & list_of_arguments)
{
    Decompression_parameters_t params;
    size_t cur_argument = 0;

    while((! list_of_arguments.empty()) &&
	  list_of_arguments.at(cur_argument)[0] == '-') {

	if(list_of_arguments.at(cur_argument) == "-v") {

	    params.verbose_flag = true;
	    cur_argument++;

	} else {

	    if(list_of_arguments.at(cur_argument) == "-")
		break;

	    std::cerr << PROGRAM_NAME
		      << ": unknown flag \"" 
		      << list_of_arguments.at(cur_argument) 
		      << "\"\n";
	    std::exit(1);

	}
    }

    if(list_of_arguments.size() - cur_argument != 2) {

	std::cerr << PROGRAM_NAME
		  << ": wrong number of arguments. Run \"squeezer decompress help\".\n";
	std::exit(1);

    }

    const std::string input_file_name = list_of_arguments.at(cur_argument);
    const std::string output_file_name = list_of_arguments.at(cur_argument + 1);
    
    FILE * input_file = NULL;
    bool read_from_stdin = false;
    if(input_file_name == "-") {
	read_from_stdin = true;
	input_file = stdin;
    } else {
	input_file = std::fopen(input_file_name.c_str(), "rb");
    }

    decompress_file_from_file(input_file,
			      output_file_name,
			      params);

    if(! read_from_stdin) {
	std::fclose(input_file);
    }
}

//////////////////////////////////////////////////////////////////////

void
dump_file_header_to_stdout(const Squeezer_file_header_t & file_header)
{
    std::printf("File format version: %d.%d (0x%04x)\n",
		MAJOR_VERSION_FROM_UINT16(file_header.program_version),
		MINOR_VERSION_FROM_UINT16(file_header.program_version),
		file_header.program_version);

    std::printf("Creation date: %04d-%02d-%02d %02d:%02d:%02d\n",
		file_header.date_year,
		file_header.date_month,
		file_header.date_day,
		file_header.time_hour,
		file_header.time_minute,
		file_header.time_second);

    std::printf("Radiometer: %s\n", file_header.radiometer.to_str().c_str());
    std::printf("Operational day: %d\n", file_header.od);
}

//////////////////////////////////////////////////////////////////////

std::string
sensible_size(uint32_t size)
{
    std::stringstream ss;
    std::vector<std::string> measure_units { "bytes", "kB", "MB", "GB", "TB" };

    uint32_t scaled_size = size;
    auto cur_unit = measure_units.begin();
    while(1) {
	if(cur_unit + 1 == measure_units.end() || scaled_size < 1024) {

	    ss << scaled_size << ' ' << *cur_unit;
	    break;

	}

	scaled_size /= 1024;
	++cur_unit;
    }

    if(scaled_size != size) {

	ss << " (" << size << " bytes)";
    }

    return ss.str();
}

//////////////////////////////////////////////////////////////////////

void
dump_chunk_header_to_stdout(size_t index,
			    const Squeezer_chunk_header_t & chunk_header)
{
    std::printf("Chunk #%lu: ", index + 1);
    switch(chunk_header.chunk_type) {
    case CHUNK_DELTA_OBT:
	std::printf("OBT times (consecutive differences)\n");
	break;
    case CHUNK_SCET_ERROR:
	std::printf("SCET times (deviation from linear interpolation with OBT times)\n");
	break;
    case CHUNK_THETA:
	std::printf("theta angle (polynomial compression)\n");
	break;
    case CHUNK_PHI:
	std::printf("phi angle (polynomial compression)\n");
	break;
    case CHUNK_PSI:
	std::printf("psi angle (polynomial compression)\n");
	break;
    case CHUNK_DIFFERENCED_DATA:
	std::printf("Scientific data (differenced)\n");
	break;
    case CHUNK_QUALITY_FLAGS:
	std::printf("Scientific flags\n");
	break;
    default:
	std::printf("Unknown chunk type, I will skip it.\n");
	return;
    }

    std::printf("    Size of the chunk: %s\n",
		sensible_size(chunk_header.number_of_bytes).c_str());
    std::printf("    Number of samples: %u\n",
		chunk_header.number_of_samples);
}


//////////////////////////////////////////////////////////////////////

void
run_statistics_task(const std::vector<std::string> & list_of_arguments)
{
    if(list_of_arguments.empty()) {

	std::cerr << PROGRAM_NAME
		  << ": you must supply at least one file name on the command line.\n"
		  << PROGRAM_NAME
		  << ": run \"squeezer help statistics\" for more information.\n";
	std::exit(1);

    }
    const std::string input_file_name = list_of_arguments.at(0);

    FILE * input_file = std::fopen(input_file_name.c_str(), "rb");
    if(input_file == NULL) {

	std::cerr << PROGRAM_NAME
		  << ": unable to open file \""
		  << input_file_name
		  << "\", reason:\n";
	std::cerr << PROGRAM_NAME
		  << ": "
		  << std::strerror(errno)
		  << '\n';
	std::exit(1);

    }

    Squeezer_file_header_t file_header(SQZ_NO_DATA);
    file_header.read_from_file(input_file);
    if(! file_header.is_valid()) {
	std::cerr << PROGRAM_NAME
		  << ": file \""
		  << input_file_name
		  << "\" does not seem to have been created by \"squeezer\".\n";
	std::exit(1);
    }

    dump_file_header_to_stdout(file_header);

    for(size_t chunk_idx = 0; 
	chunk_idx < file_header.number_of_chunks; 
	++chunk_idx) {

	Squeezer_chunk_header_t chunk_header;
	chunk_header.read_from_file(input_file);
	if(! chunk_header.is_valid()) {
	    std::cerr << PROGRAM_NAME
		      << ": file \""
		      << input_file_name
		      << "\" seems to have been damaged, chunk headers are inconsistent.\n";
	    std::exit(1);
	}

	dump_chunk_header_to_stdout(chunk_idx, chunk_header);

	if(std::fseek(input_file, chunk_header.number_of_bytes, SEEK_CUR) < 0) {

	    std::cerr << PROGRAM_NAME
		      << ": unable to move within file \""
		      << input_file_name
		      << "\", reason:\n";
	    std::cerr << PROGRAM_NAME
		      << ": "
		      << std::strerror(errno)
		      << '\n';
	    std::exit(1);

	}
    }

    std::fclose(input_file);
}

//////////////////////////////////////////////////////////////////////

void
run_program(const std::vector<std::string> & list_of_arguments)
{
    if(list_of_arguments.size() == 0 ||
	list_of_arguments.at(0) == "help" ||
	list_of_arguments.at(0) == "--help" ||
	list_of_arguments.at(0) == "-h") {
	print_help(list_of_arguments);
	std::exit(0);
    }

    std::vector<std::string> command_args = list_of_arguments;
    command_args.erase(command_args.begin());

    if(list_of_arguments.at(0) == "compress") {

	run_compression_task(command_args);

    } else if(list_of_arguments.at(0) == "decompress") {

	run_decompression_task(command_args);

    } else if(list_of_arguments.at(0) == "statistics") {

	run_statistics_task(command_args);

    } else {
	std::cerr << PROGRAM_NAME 
		  << ": unknown command \"" 
		  << list_of_arguments.at(0)
		  << "\"\n";
	std::exit(1);
    }
}

//////////////////////////////////////////////////////////////////////

int
main(int argc, const char *argv[])
{
    std::vector<std::string> list_of_arguments;
    for(int i = 1; i < argc; ++i) {
	list_of_arguments.push_back(argv[i]);
    }

#if HAVE_TOODI
    ObjectHandle init_handle;
    toodiInitializeLLIO("TOODI%file", &init_handle);
#endif

    run_program(list_of_arguments);

#if HAVE_TOODI
    toodiCloseLLIO(init_handle);
#endif
}
