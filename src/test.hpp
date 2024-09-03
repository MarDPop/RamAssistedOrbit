#pragma once

#include <cmath>
#include <exception>
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <functional>

class test_exception : public std::runtime_error {
    std::string msg;
public:
    test_exception(const std::string& arg, const char* file, int line) :
    std::runtime_error(arg) {
        std::ostringstream o;
        o << file << ":" << line << ": " << arg;
        msg = o.str();
    }
    ~test_exception() throw() {}
    const char* what() const throw() {
        return msg.c_str();
    }
};

#define throw_line(arg) throw test_exception(arg, __FILE__, __LINE__);

class Test
{
protected:

    std::vector<std::function<bool()>> _test_cases;

public:

    const std::string name;

    template<typename T>
    static bool assertEqual(T a, T b, const char* file, int line) 
    {
        if(a == b)
        {
            return true;
        }
        throw test_exception("Assertion Failed! Value expected: " + std::to_string(a) 
            + " Value was: " + std::to_string(b), file, line);
        return false;
    }

    template<typename T>
    static bool assertLess(T a, T b, const char* file, int line) 
    {
        if(a > b)
        {
            return true;
        }
        throw test_exception("Assertion Failed! Value should be less than: " + std::to_string(a) + " Value was: " 
            + std::to_string(b), file, line);
        return false;
    }

    template<typename T>
    static bool assertGreater(T a, T b, const char* file, int line) 
    {
        if(a < b)
        {
            return true;
        }
        throw test_exception("Assertion Failed! Value should be greater than: " + std::to_string(a) + " Value was: " 
            + std::to_string(b), file, line);
        return false;
    }

    static bool assertClose(double a, double b, double tol, const char* file, int line)
    {
        if(fabs(a-b) < tol)
        {
            return true;
        }
        throw test_exception("Assertion Failed! Value expected: " + std::to_string(a) + " Value was: " 
            + std::to_string(b) + " Tolerance = "  + std::to_string(tol), file, line);
        return false;
    }

    static bool assertArrayClose(double* a, double* b, unsigned n, double tol, const char* file, int line)
    {
        bool ret = true;
        for(auto i = 0u; i < n; i++) {
            if(fabs(a[i] - b[i]) > tol)
            {
                throw test_exception("Assertion Failed! Value expected: " + std::to_string(a[i]) +" Value was: " 
                    + std::to_string(b[i]) +" Tolerance = " + std::to_string(tol), file, line);
                ret = false;
            }
        }
        return ret;
    }

    static bool assertZero(double a, double tol, const char* file, int line)
    {
        if(fabs(a) < tol)
        {
            return true;
        }
        throw test_exception("Assertion Failed! Expected 0, Value was: " + std::to_string(a) + " Tolerance = " 
            + std::to_string(tol), file, line);
        return false;
    }

    Test(std::string name_) : name(name_) {}

    bool run() 
    {
        bool passed = true;
        for(auto test : _test_cases) {
            try 
            {
                if(!test()) 
                {
                    passed = false;
                }
            }
            catch (const test_exception& e)
            {
                std::cerr << name << " failed: " << e.what() << std::endl;
                passed = false;
            }
        }
        return passed;
    }
};

#define assertEqual(a, b) Test::assertEqual(a, b, __FILE__, __LINE__); 

#define assertArrayClose(a, b, n, tol) Test::assertArrayClose(a, b, n, tol, __FILE__, __LINE__); 

#define assertClose(a, b, tol) Test::assertClose(a, b, tol, __FILE__, __LINE__); 

#define assertZero(a, tol) Test::assertZero(a, tol, __FILE__, __LINE__); 

class TestMath : public virtual Test
{
public:
    TestMath() : Test("TestMath") {
        this->_test_cases.emplace_back(std::bind(&TestMath::testZYRotation, this));
    }

    bool testZYRotation();
};

class TestEarth : public virtual Test
{
public:
    TestEarth() : Test("TestEarth") {
        this->_test_cases.emplace_back(std::bind(&TestEarth::testLLA, this));
    }

    bool testLLA();
};

class TestPhysics : public virtual Test
{
    bool test_area_ratio();

    bool test_aero_quantities();

public:

    TestPhysics() : Test("TestPhysics") {
        this->_test_cases.emplace_back(std::bind(&TestPhysics::test_area_ratio, this));
        this->_test_cases.emplace_back(std::bind(&TestPhysics::test_aero_quantities, this));
    }
};

class TestRamjet : public virtual Test
{
public:
    TestRamjet() : Test("TestRamjet") {
        this->_test_cases.emplace_back(std::bind(&TestRamjet::testRamjet, this));
    }

    bool testRamjet();
};

class Tests
{
public:

    static std::vector<Test> tests;

    static void runAll();
};