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

#ifndef FILE_IO_HPP
#define FILE_IO_HPP

uint8_t read_uint8(FILE * in);
uint16_t read_uint16(FILE * in);
uint32_t read_uint32(FILE * in);
uint64_t read_uint64(FILE * in);
double read_double(FILE * in);

void write_uint8(FILE * out, uint8_t value);
void write_uint16(FILE * out, uint16_t value);
void write_uint32(FILE * out, uint32_t value);
void write_uint64(FILE * out, uint64_t value);
void write_double(FILE * out, double value);

#endif

