#pragma once

#include <algorithm>
#include <stdexcept>

template<typename T>
struct wrapped_array 
{
    T* const data;

    wrapped_array(unsigned size) : data(new T[size]) {}

    ~wrapped_array()
    {
        delete[] data;
    }

    const T& operator[](unsigned idx) const 
    {
        return data[idx];
    }

    T& operator[](unsigned idx) 
    {
        return data[idx];
    }
};

template<typename T>
struct dynamic_array 
{
    T* const data;

    const unsigned size;

    dynamic_array(unsigned size_) : data(new T[size_]), size(size_) {}

    dynamic_array(const dynamic_array& other) : data(new T[other.size]), size(other.size) 
    {
        std::copy_n(other.data, size, this->data);
    }

    dynamic_array(dynamic_array&& other) : data(std::move(other.data)), size(other.size) {}

    dynamic_array& operator=(const dynamic_array& other) // copy assignment
    {
        if(other.size != this->size) 
        {
            throw std::range_error("size of dynamic arrays must be the same to copy");
        }
        std::copy_n(other.data, size, this->data);
    }

    dynamic_array& operator=(dynamic_array&& other) 
    {
        if(other.size != this->size) 
        {
            throw std::range_error("size of dynamic arrays must be the same to copy");
        }
        std::copy_n(other.data, size, this->data);
    }

    ~dynamic_array()
    {
        delete[] data;
    }

    const T* begin() const 
    {
        return data;
    }

    T* begin() 
    {
        return data;
    }

    const T* end() const 
    {
        return data + size;
    }

    T* end() 
    {
        return data + size;
    }

    const T& operator[](unsigned idx) const 
    {
        return data[idx];
    }

    T& operator[](unsigned idx) 
    {
        return data[idx];
    }
};