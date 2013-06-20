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

#include <gsl/gsl_math.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_statistics_double.h>

#include "poly_fit_encoding.hpp"

//////////////////////////////////////////////////////////////////////

void
Frame_t::write_to_buffer(Byte_buffer_t & output_buffer)
{
    output_buffer.append_uint8(num_of_elements);
    output_buffer.append_uint8(parameters.size());

    for(auto cur_parameter : parameters) {
	output_buffer.append_float(cur_parameter);
    }
}

//////////////////////////////////////////////////////////////////////

void
Frame_t::read_from_buffer(Byte_buffer_t & input_buffer)
{
    num_of_elements = input_buffer.read_uint8();

    size_t num_of_parameters = input_buffer.read_uint8();
    parameters.resize(num_of_parameters);
    for(size_t idx = 0; idx < num_of_parameters; ++idx) {
	parameters[idx] = input_buffer.read_float();
    }
}

//////////////////////////////////////////////////////////////////////

struct Multifit_workspace {
    gsl_matrix * X, * cov;
    gsl_vector * y, * c;
    gsl_vector * residuals;
    gsl_multifit_linear_workspace * gsl_workspace;

    size_t num_of_elements;
    size_t num_of_parameters;

    Multifit_workspace() {
	X = cov = NULL;
	y = c = NULL;
	residuals = NULL;
	gsl_workspace = NULL;

	num_of_elements = 0;
	num_of_parameters = 0;
    }

    ~Multifit_workspace() {
	free();
    }

    void free() {
	if(gsl_workspace != NULL)
	    gsl_multifit_linear_free(gsl_workspace);

	if(X != NULL)
	    gsl_matrix_free(X);

	if(cov != NULL)
	    gsl_matrix_free(cov);

	if(y != NULL)
	    gsl_vector_free(y);

	if(c != NULL)
	    gsl_vector_free(c);

	if(residuals != NULL)
	    gsl_vector_free(residuals);
    }

    void set_sizes(size_t a_num_of_elements,
		   size_t a_num_of_parameters) {

	if(a_num_of_elements == num_of_elements &&
	   a_num_of_parameters == num_of_parameters) {
	    return;
	}

	free();

	num_of_elements = a_num_of_elements;
	num_of_parameters = a_num_of_parameters;

	X = gsl_matrix_alloc(num_of_elements, num_of_parameters);
	cov = gsl_matrix_alloc(num_of_parameters, num_of_parameters);

	y = gsl_vector_alloc(num_of_elements);
	c = gsl_vector_alloc(num_of_parameters);

	residuals = gsl_vector_alloc(num_of_elements);

	for(size_t idx = 0; idx < num_of_elements; ++idx) {

	    for(size_t param_idx = 0; param_idx < num_of_parameters; ++param_idx) {
		gsl_matrix_set(X, idx, param_idx, gsl_pow_int(idx, param_idx));
	    }

	}

	gsl_workspace = gsl_multifit_linear_alloc(num_of_elements,
						  num_of_parameters);
    }

};

//////////////////////////////////////////////////////////////////////

void
encode_elements_in_frame(Multifit_workspace & workspace,
			 Frame_t & frame,
			 const std::vector<double>::const_iterator first_value,
			 double & abs_error)
{
    /* The code assumes that frame.num_of_elements has already been
     * initialized, and that frame.parameters has been resized to the
     * number of parameters. The only purpose of this function is to
     * fill frame.parameters with meaningful floating-point values. */

    const size_t num_of_elements = frame.num_of_elements;
    const size_t num_of_parameters = frame.parameters.size();

    /* We're solving the system
     *
     *  y = X c
     *
     * where "X" is a n times p matrix, "y" is a n-vector and "c" a
     * p-vector (containing the unknown coefficients of the fit).
     */

    /* In this loop we remove abrupt changes in the angles when
     * running across the \pm 360° boundary. In this way polynomial
     * interpolation does not have to deal with sharp
     * discontinuities. */

    gsl_vector_set(workspace.y, 0, *first_value);

    double offset = 0.0;
    for(size_t idx = 1; idx < num_of_elements; ++idx) {
	double diff_with_previous = *(first_value + idx) - *(first_value + idx - 1);
	if(diff_with_previous > M_PI) {
	    offset -= M_PI * 2.0;
	} else if(diff_with_previous < -M_PI) {
	    offset += M_PI * 2.0;
	}

	gsl_vector_set(workspace.y, idx, *(first_value + idx) + offset);
    }

    double chi_squared;
    gsl_multifit_linear(workspace.X, 
			workspace.y,
			workspace.c, 
			workspace.cov, 
			&chi_squared, 
			workspace.gsl_workspace);

    gsl_multifit_linear_residuals(workspace.X,
				  workspace.y,
				  workspace.c,
				  workspace.residuals);

    abs_error = 0.0;
    for(size_t idx = 0; idx < num_of_elements; ++idx) {
	double error = std::fabs(gsl_vector_get(workspace.residuals, idx));
	if(idx == 0 || error > abs_error)
	    abs_error = error;
    }


    for(size_t idx = 0; idx < num_of_parameters; ++idx) {
	frame.parameters[idx] = gsl_vector_get(workspace.c, idx);
    }
}

//////////////////////////////////////////////////////////////////////

void
poly_fit_encode(const std::vector<double> & values,
		size_t elements_per_frame,
		unsigned int num_of_parameters,
		double max_abs_error,
		Byte_buffer_t & output_buffer,
		size_t & num_of_frames,
		size_t & num_of_frames_encoded_directly)
{
    Vector_of_frames_t frames;
    Multifit_workspace workspace;

    num_of_frames = 0;
    num_of_frames_encoded_directly = 0;
    for(size_t cur_idx = 0;
	cur_idx < values.size();
	cur_idx += elements_per_frame, ++num_of_frames) {

	Frame_t cur_frame;
	cur_frame.num_of_elements = std::min(elements_per_frame,
					     values.size() - cur_idx);

	double direct_encoding = true;

	if(cur_frame.num_of_elements > num_of_parameters) {

	    // Apply the polynomial fitting compression
	    workspace.set_sizes(cur_frame.num_of_elements,
				num_of_parameters);

	    cur_frame.parameters.resize(num_of_parameters);

	    double abs_error;
	    encode_elements_in_frame(workspace, 
				     cur_frame, 
				     values.begin() + cur_idx,
				     abs_error);

	    direct_encoding = abs_error >= max_abs_error;
	}

	if(direct_encoding) {

	    cur_frame.parameters.resize(cur_frame.num_of_elements);
	    // There are too few elements left, just copy them as they are
	    std::copy(values.begin() + cur_idx,
		      values.begin() + cur_idx + cur_frame.num_of_elements,
		      cur_frame.parameters.begin());

	    ++num_of_frames_encoded_directly;

	}

	cur_frame.write_to_buffer(output_buffer);
    }
}

//////////////////////////////////////////////////////////////////////

void
poly_fit_decode(size_t num_of_elements_to_decode,
		Byte_buffer_t & input_buffer,
		std::vector<double> & values)
{
    values.resize(num_of_elements_to_decode);

    size_t cur_idx = 0;
    while(cur_idx < num_of_elements_to_decode) {
	Frame_t cur_frame;
	cur_frame.read_from_buffer(input_buffer);

	if(cur_frame.num_of_elements > cur_frame.parameters.size()) {

	    for(size_t i = 0; i < cur_frame.num_of_elements; ++i) {
		values[cur_idx + i] =
		    gsl_poly_eval(cur_frame.parameters.data(),
				  cur_frame.parameters.size(),
				  i);
	    }

	} else {

	    std::copy(cur_frame.parameters.begin(),
		      cur_frame.parameters.end(),
		      values.begin() + cur_idx);

	}

	cur_idx += cur_frame.num_of_elements;
    }
}
