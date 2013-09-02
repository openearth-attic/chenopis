#!/usr/bin/env python
import os
import logging
import unittest
import inspect


import openearthtools.modelapi.bmi

logger = logging.getLogger(__name__)
LibSwan = openearthtools.modelapi.bmi.BMI

__file__ = inspect.getfile(inspect.currentframe())
DIRNAME = os.path.dirname(os.path.abspath(__file__))
LIBSWANNAME = os.path.join(DIRNAME, '../src/.libs/libchenopis.dylib')



class SwanBMI(unittest.TestCase):
    """Tests that don't require initialization of SWAN"""
    def setUp(self):
        self.swan = LibSwan(LIBSWANNAME, os.path.join(DIRNAME, '../cases/a11refra/model_io'))
        self.swan.initialize("INPUT")
    def tearDown(self):
        self.swan.finalize()
        del self.swan
    def test_main(self):
        """Test if main function can run"""
        self.swan.main()
        logger.info("main ran")

    def test_initialize(self):
        """Test if initialize function can run"""
        # setup and tear down, run initialize

    def test_update(self):
        """Test if update function can run """
        self.swan.update(2.0)
    def test_n_attributes(self):
        """Test if update function can run """
        n = self.swan.get_n_attributes()
        self.assertEqual(n, 40)
    def test_attribute_name(self):
        """Test if we can get an attribute name"""
        name = self.swan.get_attribute_name(19)
        self.assertEqual(name, "CFL")
    def test_attribute_type(self):
        """Test if we can get an attribute type"""
        type_ = self.swan.get_attribute_type("CFL")
        self.assertEqual(type_, "BMI_INT")
    def test_double_attribute(self):
        """Test if we can get a double type"""
        value = self.swan.get_double_attribute("CFL")
        self.assertEqual(-1.0, value)
    def test_int_attribute(self):
        """Test if we can get an int attribute"""
        value = self.swan.get_int_attribute("MSC")
        self.assertEqual(41, value)
    def test_string_attribute(self):
        """Test if we can get a string attribute"""
        value = self.swan.get_string_attribute("CFL")
        self.assertEqual('some value', value)
    def test_n_vars(self):
        """Test if we can get the number of variables """
        n = self.swan.get_n_vars()
        self.assertEqual(n, 82)
    def test_get_time_step(self):
        """Test if we can get the timestep"""
        # 0 before read
        # 1000... after read
        dt = self.swan.get_time_step()
        self.assertEqual(dt, 10000000000.0)
    def test_var_name(self):
        """Test if we can get a variable name"""
        name = self.swan.get_var_name()
        self.assertEqual(name, "VEL")
    def test_get_var_i(self):
        """Test if we can get a variable number"""
        i = self.swan.get_var_i("VEL")
        self.assertEqual(i, 5)
    # def test_var_type(self):
    #     """Test if we can get the variable type"""
    #     type_ = self.swan.get_var_type("VEL")
    #     self.assertEqual(type_, "BMI_DOUBLE")
    def test_var_unit(self):
        """Test if we can get the variable unit"""
        unit = self.swan.get_var_unit("VEL")
        self.assertEqual(unit, "m/s")
    def test_update(self):
        """Test if we can do a timestep"""
        self.swan.update(0.1)
    def test_3d_float(self):

        msc = self.swan.get_int_attribute('MSC')
        mdc = self.swan.get_int_attribute('MDC')
        mcgrd = self.swan.get_int_attribute('MCGRD')
        print mdc, msc,  mcgrd
        arraytype = ndpointer(dtype='float32',
                              ndim=3,
                              shape=(mdc, msc, mcgrd)[::-1],
                              flags='F')
        x = self.swan.get_3d_float('AC2', arraytype)
        self.assertEqual(0, x[-1,-1,-1])

if __name__=='__main__':
    unittest.main()
