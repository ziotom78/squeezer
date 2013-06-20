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

bin_PROGRAMS = squeezer

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

if TOODI_PRESENT

AM_CPPFLAGS += $(TOODI_CFLAGS)
LIBS += $(TOODI_LIBS)

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

AM_CPPFLAGS += $(CPPUNIT_CFLAGS)
LIBS += $(CPPUNIT_LIBS)

endif
