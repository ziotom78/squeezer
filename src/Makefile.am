# Squeezer - compress LFI detector pointings and differenced data
# Copyright (C) 2013 Maurizio Tomasi (Planck collaboration)
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
# 02110-1301, USA.

AM_CPPFLAGS = --pedantic
PROGRAMS_TO_BUILD = squeezer

######################################################################

if HPIXLIB_PRESENT

PROGRAMS_TO_BUILD += hit_map

hit_map_SOURCES = \
	byte_buffer.cpp \
	common_defs.cpp \
	data_structures.cpp \
	decompress.cpp \
	datadiff.cpp \
	detpoint.cpp \
	file_io.cpp \
	poly_fit_encoding.cpp \
	run_length_encoding.cpp \
	hit_map.cpp

hit_map_CPPFLAGS = $(GSL_CFLAGS) $(HPIXLIB_CFLAGS)
hit_map_LDADD = $(GSL_LDFLAGS) $(HPIXLIB_LIBS)

if TOODI_PRESENT

hit_map_CPPFLAGS += $(TOODI_CFLAGS)
hit_map_LDADD += $(TOODI_LIBS)

endif

endif

######################################################################

bin_PROGRAMS = $(PROGRAMS_TO_BUILD)

squeezer_SOURCES = \
	byte_buffer.cpp \
	common_defs.cpp \
	compress.cpp \
	data_structures.cpp \
	datadiff.cpp \
	decompress.cpp \
	detpoint.cpp \
	file_io.cpp \
	help.cpp \
	main.cpp \
	poly_fit_encoding.cpp \
	run_length_encoding.cpp \
	statistics.cpp

squeezer_CPPFLAGS = $(GSL_CFLAGS)
squeezer_LIBS = $(GSL_LDFLAGS)
squeezer_LDFLAGS = $(GSL_LDFLAGS)

if TOODI_PRESENT

squeezer_CPPFLAGS += $(TOODI_CFLAGS)
squeezer_LIBS += $(TOODI_LIBS)
squeezer_LDFLAGS += $(TOODI_LDFLAGS)

endif

######################################################################

if CPPUNIT_PRESENT

TESTS = check_program

check_PROGRAMS = $(TESTS)

check_program_SOURCES = \
	byte_buffer.cpp \
	check_program.cpp \
	common_defs.cpp \
	data_structures.cpp \
	file_io.cpp \
	poly_fit_encoding.cpp \
	run_length_encoding.cpp \
	statistics.cpp

check_program_CPPFLAGS = $(GSL_CFLAGS) $(CPPUNIT_CFLAGS)
check_program_LIBS = $(GSL_LIBS) $(CPPUNIT_LIBS)
check_program_LDFLAGS = $(GSL_LDFLAGS) $(CPPUNIT_LDFLAGS)

endif
