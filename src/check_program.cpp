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

#include <cmath>
#include <vector>
#include <regex>
#include <algorithm>

#include <cppunit/TestAssert.h>
#include <cppunit/TestCaller.h>
#include <cppunit/TestFixture.h>
#include <cppunit/ui/text/TestRunner.h>

#include "common_defs.hpp"
#include "statistics.hpp"
#include "run_length_encoding.hpp"
#include "poly_fit_encoding.hpp"
#include "byte_buffer.hpp"
#include "data_structures.hpp"

//////////////////////////////////////////////////////////////////////

class Radiometer_test : public CppUnit::TestFixture {
public:
    static void check_radiometer(uint8_t reference_horn, 
			  uint8_t reference_arm,
			  Radiometer_t radiometer) {
	CPPUNIT_ASSERT_EQUAL((int) reference_horn, (int) radiometer.horn);
	CPPUNIT_ASSERT_EQUAL((int) reference_arm, (int) radiometer.arm);
    }

    void testFullRadiometerNameParsing() {
	Radiometer_t radiometer;

	radiometer.parse_from_name("LFI24S");
	check_radiometer(24, 1, radiometer);

	radiometer.parse_from_name("LFI18M");
	check_radiometer(18, 0, radiometer);

	radiometer.parse_from_name("LFI28S");
	check_radiometer(28, 1, radiometer);
    }

    void testShortRadiometerNameParsing() {
	Radiometer_t radiometer;

	radiometer.parse_from_name("24S");
	check_radiometer(24, 1, radiometer);

	radiometer.parse_from_name("18M");
	check_radiometer(18, 0, radiometer);

	radiometer.parse_from_name("28S");
	check_radiometer(28, 1, radiometer);
    }

    void testWrongRadiometerNameParsing() {
	Radiometer_t radiometer;

	CPPUNIT_ASSERT_THROW(radiometer.parse_from_name("this_is_not_valid"),
			     std::runtime_error);
    }

    static CppUnit::Test * suite() {
	CppUnit::TestSuite * suite = new CppUnit::TestSuite("Radiometer_test");
	suite->addTest(new CppUnit::TestCaller<Radiometer_test>(
			   "testFullRadiometerNameParsing", 
			   &Radiometer_test::testFullRadiometerNameParsing));
	suite->addTest(new CppUnit::TestCaller<Radiometer_test>(
			   "testShortRadiometerNameParsing", 
			   &Radiometer_test::testShortRadiometerNameParsing));
	suite->addTest(new CppUnit::TestCaller<Radiometer_test>(
			   "testWrongRadiometerNameParsing", 
			   &Radiometer_test::testWrongRadiometerNameParsing));
	return suite;
    }
};

//////////////////////////////////////////////////////////////////////

class Bytestream_test : public CppUnit::TestFixture {
private:
    bytestream_t bytestream;

public:
    void testBytestreamFromWords() {
	std::vector<uint16_t> words = { 0x1234, 0x5678, 0x9ABC };

	vector_to_bytestream(words, bytestream);

	CPPUNIT_ASSERT_EQUAL((size_t) bytestream.size(), 6UL);

	CPPUNIT_ASSERT_EQUAL((byte_t) 0x34, bytestream[0]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x12, bytestream[1]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x78, bytestream[2]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x56, bytestream[3]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0xBC, bytestream[4]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x9A, bytestream[5]);
    }

    void testBytestreamFromDoubleWords() {
	std::vector<uint32_t> double_words = { 0x12345678, 0x13579BDF };

	vector_to_bytestream(double_words, bytestream);

	CPPUNIT_ASSERT_EQUAL((size_t) bytestream.size(), 8UL);

	CPPUNIT_ASSERT_EQUAL((byte_t) 0x78, bytestream[0]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x56, bytestream[1]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x34, bytestream[2]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x12, bytestream[3]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0xDF, bytestream[4]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x9B, bytestream[5]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x57, bytestream[6]);
	CPPUNIT_ASSERT_EQUAL((byte_t) 0x13, bytestream[7]);
    }

    static CppUnit::Test * suite() {
	CppUnit::TestSuite * suite = new CppUnit::TestSuite("Bytestream_test");
	suite->addTest(new CppUnit::TestCaller<Bytestream_test>(
			   "testBytestreamFromWords", 
			   &Bytestream_test::testBytestreamFromWords));
	suite->addTest(new CppUnit::TestCaller<Bytestream_test>(
			   "testBytestreamFromDoubleWords", 
			   &Bytestream_test::testBytestreamFromDoubleWords));
	return suite;
    }
};

////////////////////////////////////////////////////////////////////

class Frequency_table_test : public CppUnit::TestFixture {
private:
    bytestream_t bytestream;
    frequency_table_t freq_table;

public:
    void setUp() {
	bytestream.push_back(3);
	bytestream.push_back(3);
	bytestream.push_back(2);
	bytestream.push_back(6);
	bytestream.push_back(4);
	bytestream.push_back(7);
	bytestream.push_back(6);
	bytestream.push_back(3);
	bytestream.push_back(5);
	bytestream.push_back(1);

	build_frequency_table(bytestream, freq_table);
    }

    void testLength() {
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table.size(), 7U);
    }

    void testInclusion() {
	CPPUNIT_ASSERT(freq_table.find(1) != freq_table.end());
	CPPUNIT_ASSERT(freq_table.find(2) != freq_table.end());
	CPPUNIT_ASSERT(freq_table.find(3) != freq_table.end());
	CPPUNIT_ASSERT(freq_table.find(4) != freq_table.end());
	CPPUNIT_ASSERT(freq_table.find(5) != freq_table.end());
	CPPUNIT_ASSERT(freq_table.find(6) != freq_table.end());
	CPPUNIT_ASSERT(freq_table.find(7) != freq_table.end());
    }

    void testFrequencies() {
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table[1], 1U);
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table[2], 1U);
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table[3], 3U);
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table[4], 1U);
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table[5], 1U);
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table[6], 2U);
	CPPUNIT_ASSERT_EQUAL((unsigned) freq_table[7], 1U);
    }

    void testEntropy() {
	/* The result (computed analytically) is 
	 *
	 * log_2(10) / 2 + log_2(5) / 5 + 3 log_2(10/3) / 10
	 */
	CPPUNIT_ASSERT(std::fabs(entropy_from_frequency_table(freq_table) -
				 2.64643934467102) < 1e-6);
    }

    static CppUnit::Test * suite() {
	CppUnit::TestSuite * suite = new CppUnit::TestSuite("Frequency_table_test");
	suite->addTest(new CppUnit::TestCaller<Frequency_table_test>(
			   "testLength", 
			   &Frequency_table_test::testLength));
	suite->addTest(new CppUnit::TestCaller<Frequency_table_test>(
			   "testInclusion",
			   &Frequency_table_test::testInclusion));
	suite->addTest(new CppUnit::TestCaller<Frequency_table_test>(
			   "testFrequencies",
			   &Frequency_table_test::testFrequencies));
	suite->addTest(new CppUnit::TestCaller<Frequency_table_test>(
			   "testEntropy",
			   &Frequency_table_test::testEntropy));
	return suite;
    }
};

////////////////////////////////////////////////////////////////////

class RLE_test : public CppUnit::TestFixture {
public:
    void testRLECompression() {
	std::vector<uint32_t> input_stream;
	Byte_buffer_t output_stream;

	input_stream.push_back(5);
	input_stream.push_back(5);
	input_stream.push_back(5);
	input_stream.push_back(6);
	input_stream.push_back(6);
	input_stream.push_back(4);
	input_stream.push_back(4);
	input_stream.push_back(3);
	input_stream.push_back(2);

	rle_compression(input_stream.data(),
			input_stream.size(),
			output_stream);

	CPPUNIT_ASSERT_EQUAL(40U, (unsigned int) output_stream.buffer.size());

	CPPUNIT_ASSERT_EQUAL((uint32_t) 3, output_stream.read_uint32());
	CPPUNIT_ASSERT_EQUAL((uint32_t) 5, output_stream.read_uint32());
			                
	CPPUNIT_ASSERT_EQUAL((uint32_t) 2, output_stream.read_uint32());
	CPPUNIT_ASSERT_EQUAL((uint32_t) 6, output_stream.read_uint32());
			                
	CPPUNIT_ASSERT_EQUAL((uint32_t) 2, output_stream.read_uint32());
	CPPUNIT_ASSERT_EQUAL((uint32_t) 4, output_stream.read_uint32());
			                
	CPPUNIT_ASSERT_EQUAL((uint32_t) 1, output_stream.read_uint32());
	CPPUNIT_ASSERT_EQUAL((uint32_t) 3, output_stream.read_uint32());
			                
	CPPUNIT_ASSERT_EQUAL((uint32_t) 1, output_stream.read_uint32());
	CPPUNIT_ASSERT_EQUAL((uint32_t) 2, output_stream.read_uint32());	
    }

    static CppUnit::Test * suite() {
	CppUnit::TestSuite * suite = new CppUnit::TestSuite("RLE_test");
	suite->addTest(new CppUnit::TestCaller<RLE_test>(
			   "testRLECompression", 
			   &RLE_test::testRLECompression));
	return suite;
    }
};

////////////////////////////////////////////////////////////////////

class Poly_fit_encoder_test : public CppUnit::TestFixture {
private:
    Vector_of_frames_t vector_of_frames;

public:
    void compare_frames(const Frame_t & frame1, 
			const Frame_t & frame2) {

	CPPUNIT_ASSERT_EQUAL((int) frame1.num_of_elements,
			     (int) frame2.num_of_elements);

	CPPUNIT_ASSERT_EQUAL((int) frame1.parameters.size(),
			     (int) frame2.parameters.size());

	CPPUNIT_ASSERT_MESSAGE("Mismatch in the value of the parameters "
			       "in two Frame_t objects",
			       std::equal(frame1.parameters.begin(),
					  frame1.parameters.end(),
					  frame2.parameters.begin()));

    }

    void setUp() {
	std::vector<double> params1 { 1.0, 2.0 };
	vector_of_frames.push_back(Frame_t(10, params1));

	std::vector<double> params2 { 10.0 };
	vector_of_frames.push_back(Frame_t(6, params2));
    }

    void testFramesToRawBuffer() {
	Byte_buffer_t raw_buffer;
	for(auto & cur_frame : vector_of_frames) {
	    cur_frame.write_to_buffer(raw_buffer);
	}

	CPPUNIT_ASSERT_EQUAL((int) 16,
			     (int) raw_buffer.size());

	// (Fake) elements per frame
	CPPUNIT_ASSERT_EQUAL(10, (int) raw_buffer.read_uint8());
	// Number of parameters
	CPPUNIT_ASSERT_EQUAL(2, (int) raw_buffer.read_uint8());
	// The parameters
	CPPUNIT_ASSERT_EQUAL((float) 1.0, raw_buffer.read_float());
	CPPUNIT_ASSERT_EQUAL((float) 2.0, raw_buffer.read_float());

	// (Fake) elements per frame
	CPPUNIT_ASSERT_EQUAL(6, (int) raw_buffer.read_uint8());
	// Number of parameters
	CPPUNIT_ASSERT_EQUAL(1, (int) raw_buffer.read_uint8());
	// The parameter
	CPPUNIT_ASSERT_EQUAL((float) 10.0, raw_buffer.read_float());
    }

    void testRawBufferToFrames() {
	Byte_buffer_t raw_buffer;
	for(auto & cur_frame : vector_of_frames) {
	    cur_frame.write_to_buffer(raw_buffer);
	}

	Frame_t test_frame;

	test_frame.read_from_buffer(raw_buffer);
	compare_frames(vector_of_frames[0], test_frame);

	test_frame.read_from_buffer(raw_buffer);
	compare_frames(vector_of_frames[1], test_frame);
    }

    void testEncoding() {
	std::vector<double> values { 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 10.0};

	Byte_buffer_t output_buffer;
	size_t num_of_frames = 0;
	size_t num_of_uncompressed_frames = 0;
	// We fit the data against a y = a_1 + a_2 x model, with x = 0..2
	poly_fit_encode(values, 3, 2, 1e6, output_buffer,
			num_of_frames, num_of_uncompressed_frames);

	CPPUNIT_ASSERT_EQUAL((int) 3, (int) num_of_frames);
	CPPUNIT_ASSERT_EQUAL((int) 1, (int) num_of_uncompressed_frames);

	/* The total size is: 2 * (2 + 4 * 2) + (2 + 4) = 26 (the
	 * first two frames have a 2-byte header plus two 4-byte
	 * floating points, then there is the last frame which only
	 * contains one floating point).
	 */
	CPPUNIT_ASSERT_EQUAL((int) 26, (int) output_buffer.size());

	// First frame. The fit is y = 1 + x (a_1 = 1, a_2 = 1)

	// Number of elements (these are: 1, 2, 3)
	CPPUNIT_ASSERT_EQUAL((int) 3, (int) output_buffer.read_uint8());
	// Number of parameters in the fit (a_1 and a_2)
	CPPUNIT_ASSERT_EQUAL((int) 2, (int) output_buffer.read_uint8());
	// a_1
	CPPUNIT_ASSERT_EQUAL((int) 1, (int) output_buffer.read_float());
	// a_2
	CPPUNIT_ASSERT_EQUAL((int) 1, (int) output_buffer.read_float());

	// Second frame. The fit is y = 5 + 2x (a_1 = 5, a_2 = 2)

	// Number of elements (these are: 5, 7, 9)
	CPPUNIT_ASSERT_EQUAL((int) 3, (int) output_buffer.read_uint8());
	// Number of parameters in the fit (a_1 and a_2)
	CPPUNIT_ASSERT_EQUAL((int) 2, (int) output_buffer.read_uint8());
	// a_1
	CPPUNIT_ASSERT_EQUAL((int) 5, (int) output_buffer.read_float());
	// a_2
	CPPUNIT_ASSERT_EQUAL((int) 2, (int) output_buffer.read_float());

	// Last frame

	// Number of elements (only one: 10)
	CPPUNIT_ASSERT_EQUAL((int) 1, (int) output_buffer.read_uint8());
	// Number of parameters in the fit (just one, and it's not a
	// real parameter)
	CPPUNIT_ASSERT_EQUAL((int) 1, (int) output_buffer.read_uint8());
	// This is the last element
	CPPUNIT_ASSERT_EQUAL((int) 10, (int) output_buffer.read_float());
    }

    void testDecoding() {
	std::vector<double> values { 1.0, 2.0, 3.0, 5.0, 7.0, 9.0, 10.0};

	Byte_buffer_t buffer;
	size_t dummy1 = 0, dummy2 = 0;
	// We fit the data against a y = a_1 + a_2 x model, with x = 0..2
	poly_fit_encode(values, 3, 2, 1e6, buffer, dummy1, dummy2);

	std::vector<double> reconstructed;
	poly_fit_decode(values.size(), buffer, reconstructed);

	CPPUNIT_ASSERT_MESSAGE("The reconstructed data stream does not "
			       "match the original",
			       std::equal(reconstructed.begin(),
					  reconstructed.end(),
					  values.begin()));
    }

    static CppUnit::Test * suite() {
	CppUnit::TestSuite * suite = new CppUnit::TestSuite("Poly_fit_encoder_test");
	suite->addTest(new CppUnit::TestCaller<Poly_fit_encoder_test>(
			   "framesToRawBuffer",
			   &Poly_fit_encoder_test::testFramesToRawBuffer));
	suite->addTest(new CppUnit::TestCaller<Poly_fit_encoder_test>(
			   "rawBufferToFrames",
			   &Poly_fit_encoder_test::testRawBufferToFrames));
	suite->addTest(new CppUnit::TestCaller<Poly_fit_encoder_test>(
			   "testEncoding",
			   &Poly_fit_encoder_test::testEncoding));
	suite->addTest(new CppUnit::TestCaller<Poly_fit_encoder_test>(
			   "testDecoding",
			   &Poly_fit_encoder_test::testDecoding));
	return suite;
    }
};

////////////////////////////////////////////////////////////////////

class File_IO_test : public CppUnit::TestFixture {
public:
    void testFileHeaderIO() {
	FILE * f = fopen("./delete_me.bin", "wb");
	
	Squeezer_file_header_t source(SQZ_DETECTOR_POINTINGS);

	source.date_year = 2013;
	source.date_month = 12;
	source.date_day = 25;

	source.time_hour = 12;
	source.time_minute = 2;
	source.time_second = 58;

	source.radiometer.horn = 18;
	source.radiometer.arm = 1;
	source.od = 163;
	source.first_obt = 1.0;
	source.last_obt = 2.0;
	source.first_scet_in_ms = 3.0;
	source.last_scet_in_ms = 4.0;

	source.number_of_chunks = 5;

	source.write_to_file(f);
	fclose(f);

	f = fopen("./delete_me.bin", "rb");

	Squeezer_file_header_t test(SQZ_NO_DATA);
	test.read_from_file(f);

	fclose(f);

#define CHECK_FIELD(name, type)				\
	CPPUNIT_ASSERT_EQUAL((type) source.name, (type) test.name)

	CHECK_FIELD(file_type_mark[0], int);
	CHECK_FIELD(file_type_mark[1], int);
	CHECK_FIELD(file_type_mark[2], int);
	CHECK_FIELD(file_type_mark[3], int);

	CHECK_FIELD(date_year, int);
	CHECK_FIELD(date_month, int);
	CHECK_FIELD(date_day, int);

	CHECK_FIELD(time_hour, int);
	CHECK_FIELD(time_minute, int);
	CHECK_FIELD(time_second, int);

	CHECK_FIELD(radiometer.horn, int);
	CHECK_FIELD(radiometer.arm, int);

	CHECK_FIELD(od, int);

	CHECK_FIELD(first_obt, int);
	CHECK_FIELD(last_obt, int);

	CHECK_FIELD(first_scet_in_ms, int);
	CHECK_FIELD(last_scet_in_ms, int);

	CHECK_FIELD(number_of_chunks, int);

#undef CHECK_FIELD
    }

    void testChunkHeaderIO() {
	FILE * f = fopen("./delete_me.bin", "wb");
	
	Squeezer_chunk_header_t source;

	source.number_of_bytes = 16532;
	source.number_of_samples = 723465;
	source.chunk_type = 15;

	source.compression_error.min_abs_error = 1.0;
	source.compression_error.max_abs_error = 2.0;
	source.compression_error.mean_abs_error = 3.0;
	source.compression_error.mean_error = 4.0;

	source.write_to_file(f);
	fclose(f);

	f = fopen("./delete_me.bin", "rb");

	Squeezer_chunk_header_t test;
	test.read_from_file(f);

	fclose(f);

#define CHECK_FIELD(name, type)				\
	CPPUNIT_ASSERT_EQUAL((type) source.name, (type) test.name)

	CHECK_FIELD(chunk_mark[0], int);
	CHECK_FIELD(chunk_mark[1], int);
	CHECK_FIELD(chunk_mark[2], int);
	CHECK_FIELD(chunk_mark[3], int);

	CHECK_FIELD(number_of_bytes, uint64_t);
	CHECK_FIELD(number_of_samples, uint32_t);
	CHECK_FIELD(chunk_type, uint32_t);

	CHECK_FIELD(compression_error.min_abs_error, double);
	CHECK_FIELD(compression_error.max_abs_error, double);
	CHECK_FIELD(compression_error.mean_abs_error, double);
	CHECK_FIELD(compression_error.mean_error, double);

#undef CHECK_FIELD
    }

    static CppUnit::Test * suite() {
	CppUnit::TestSuite * suite = new CppUnit::TestSuite("File_IO_test");
	suite->addTest(new CppUnit::TestCaller<File_IO_test>(
			   "testFileHeaderIO", 
			   &File_IO_test::testFileHeaderIO));
	suite->addTest(new CppUnit::TestCaller<File_IO_test>(
			   "testChunkHeaderIO", 
			   &File_IO_test::testChunkHeaderIO));
	return suite;
    }
};

//////////////////////////////////////////////////////////////////////

class Byte_buffer_test : public CppUnit::TestFixture {
public:
    void testRead() {
	const uint8_t * raw_buffer = (uint8_t *)
	    "\x01" // uint8_t
	    "\x02\x03" // uint16_t
	    "\x04\x05\x06\x07" // uint32_t
	    "\x08\x09\x0A\x0B\x0C\x0D\x0E\x0F";  // uint64_t
	size_t raw_buffer_length = 15;
	Byte_buffer_t buffer(raw_buffer_length, raw_buffer);

	CPPUNIT_ASSERT_EQUAL((int) raw_buffer_length, 
			     (int) buffer.buffer.size());

	CPPUNIT_ASSERT_EQUAL((int) 0x01, (int) buffer.read_uint8());
	CPPUNIT_ASSERT_EQUAL((int) 0x0203, (int) buffer.read_uint16());
	CPPUNIT_ASSERT_EQUAL((int) 0x04050607, (int) buffer.read_uint32());
	CPPUNIT_ASSERT_EQUAL((int) 0x08090A0B0C0D0E0F, (int) buffer.read_uint64());
    }

    void testWrite() {
	const uint8_t * raw_buffer = (uint8_t *)
	    "\x01" // uint8_t
	    "\x02\x03" // uint16_t
	    "\x04\x05\x06\x07" // uint32_t
	    "\x08\x09\x0A\x0B\x0C\x0D\x0E\x0F";  // uint64_t
	Byte_buffer_t buffer;

	buffer.append_uint8(0x01);
	buffer.append_uint16(0x0203);
	buffer.append_uint32(0x04050607);
	buffer.append_uint64(0x08090A0B0C0D0E0F);

	CPPUNIT_ASSERT_EQUAL(15, (int) buffer.buffer.size());
	CPPUNIT_ASSERT(std::equal(buffer.buffer.begin(), buffer.buffer.end(), 
				  &raw_buffer[0]));
    }

    void testFloatingPoint() {
	Byte_buffer_t buffer;
	buffer.append_float(123.0);
	buffer.append_double(456.0);

	CPPUNIT_ASSERT_EQUAL((float) 123.0, buffer.read_float());
	CPPUNIT_ASSERT_EQUAL((double) 456.0, buffer.read_double());
    }

    void testReadAfterEnd() {
	Byte_buffer_t buffer;
	buffer.append_uint8(0);

	buffer.read_uint8();
	CPPUNIT_ASSERT_THROW(buffer.read_uint8(),
			     std::out_of_range);
    }

    static CppUnit::Test * suite() {
	CppUnit::TestSuite * suite = new CppUnit::TestSuite("Byte_buffer_test");
	suite->addTest(new CppUnit::TestCaller<Byte_buffer_test>(
			   "testRead",
			   &Byte_buffer_test::testRead));
	suite->addTest(new CppUnit::TestCaller<Byte_buffer_test>(
			   "testWrite",
			   &Byte_buffer_test::testWrite));
	suite->addTest(new CppUnit::TestCaller<Byte_buffer_test>(
			   "testFloatingPoint",
			   &Byte_buffer_test::testFloatingPoint));
	suite->addTest(new CppUnit::TestCaller<Byte_buffer_test>(
			   "testReadAfterEnd",
			   &Byte_buffer_test::testReadAfterEnd));

	return suite;
    }
};

////////////////////////////////////////////////////////////////////

int
main(void)
{
    CppUnit::TextUi::TestRunner runner;
    runner.addTest(Radiometer_test::suite());
    runner.addTest(Bytestream_test::suite());
    runner.addTest(Frequency_table_test::suite());
    runner.addTest(RLE_test::suite());
    runner.addTest(Poly_fit_encoder_test::suite());
    runner.addTest(Byte_buffer_test::suite());
    runner.addTest(File_IO_test::suite());
    runner.run();
    return 0;
}
