import os
import sys

ovf_py_dir = os.path.abspath(os.path.join(os.path.dirname( __file__ ), ".."))
sys.path.insert(0, ovf_py_dir)

from ovf import ovf

import numpy as np

import unittest

##########

class TestState(unittest.TestCase):
    def test_nonexistent(self):
        print("----- ovf test nonexistent")
        with ovf.ovf_file("nonexistent.ovf") as ovf_file:
            print("found:      ", ovf_file.found)
            print("is_ovf:     ", ovf_file.is_ovf)
            print("n_segments: ", ovf_file.n_segments)
            segment = ovf.ovf_segment()
            success = ovf_file.read_segment_header(0, segment)
            if success != ovf.OK:
                print("read_segment_header failed: ", ovf_file.get_latest_message())
            self.assertFalse( success == ovf.OK )
        print("----- ovf test nonexistent done")

    def test_write(self):
        print("----- ovf test writing")
        with ovf.ovf_file("testfile_py.ovf") as ovf_file:
            data = np.zeros((2, 2, 1, 3), dtype='d')
            data[0,1,0,:] = [3.0, 2.0, 1.0]
            segment = ovf.ovf_segment(n_cells=[2,2,1], valuedim=3)
            success = ovf_file.write_segment(segment, data)
            if success != ovf.OK:
                print("write_segment failed: ", ovf_file.get_latest_message())
            self.assertTrue( success == ovf.OK )
            data[0,1,0,:] = [4.0, 5.0, 6.0]
            success = ovf_file.append_segment(segment, data)
            if success != ovf.OK:
                print("append_segment failed: ", ovf_file.get_latest_message())
            self.assertTrue( success == ovf.OK )
        print("----- ovf test writing done")

        print("----- ovf test reading")
        with ovf.ovf_file("testfile_py.ovf") as ovf_file:
            print("found:      ", ovf_file.found)
            print("is_ovf:     ", ovf_file.is_ovf)
            print("n_segments: ", ovf_file.n_segments)
            segment = ovf.ovf_segment()
            success = ovf_file.read_segment_header(0, segment)
            if success != ovf.OK:
                print("read_segment_header failed: ", ovf_file.get_latest_message())
            data_shape = (segment.n_cells[0], segment.n_cells[1], segment.n_cells[2], 3)
            print("data shape: ", data_shape)
            self.assertTrue( success == ovf.OK )
            data = np.zeros(data_shape, dtype='f')
            success = ovf_file.read_segment_data(0, segment, data)
            if success != ovf.OK:
                print("read_segment_data failed: ", ovf_file.get_latest_message())
            print("first segment:  ", data[0,1,0,:])
            self.assertTrue( success == ovf.OK )
            success = ovf_file.read_segment_data(1, segment, data)
            if success != ovf.OK:
                print("read_segment_data failed: ", ovf_file.get_latest_message())
            print("second segment: ", data[0,1,0,:])
            self.assertTrue( success == ovf.OK )
        print("----- ovf test reading done")

#########

def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestState))
    return suite

if __name__ == '__main__':
    suite = suite()
    runner = unittest.TextTestRunner()
    success = runner.run(suite).wasSuccessful()
    sys.exit(not success)