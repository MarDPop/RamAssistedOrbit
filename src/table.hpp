#pragma once

#include <algorithm>

template<typename T>
class Table_Linear
{
protected:

    static const T* _search_linear(T x, T* _x)
    {
        const T* ptr = _x + _n;
        while(--ptr != _x && *ptr > x) {}
        return ptr;
    }

    static const T* _search_binary(T x)
    {
        return std::lower_bound(_x, _x + _n, x);
    }

    T* const _x;

    T* const _y;

    T* const _delta;

    T* const (*_search)(T x, T* _x);

    unsigned _n;

public:

    constexpr unsigned MAX_LINEAR_SIZE = 32u;

    Table_Linear(const T* x, const T* y, unsigned n) : _x(new T[n]), _y(new T[n]), 
        _delta(new T[n]), _search(n > MAX_LINEAR_SIZE ? &Table_Linear::_search_binary : &Table_Linear::_search_linear), _n(n)
    {
        assert(n > 1);
        memcpy(_x, x, n*sizeof(T));
        memcpy(_y, y, n*sizeof(T));
        for(auto i = 1u; i < n; i++) 
        {
            _delta[i - 1] = (_y[i] - y[i - 1])/(_x[i] - x[i - 1]);
        }
        _delta[_n - 1] = _delta[_n - 2];
    }

    ~Table_Linear() 
    {
        delete[] _x;
        delete[] _y;
        delete[] _delta;
    }

    T get(T x) const
    {
        const T* x_ptr = _search(x, _x);
        const T dx = x - *x_ptr;
        const auto dist = x_ptr - _x;
        return _y[dist] + _delta[dist]*dx;
    }

};