#pragma once

#include <algorithm>
#include <vector>

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


class nested_table
{
    struct node {

        static std::vector<double> diff(std::vector<double> values)
        {
            std::vector<double> dvalues(values.size() - 1);
            for(unsigned i = 1; i < values.size(); i++)
            {
                dvalues[i - 1] = 1.0/(values[i] - values[i - 1]);
            }
            return dvalues;
        }

        const std::vector<double> values;
        const std::vector<double> dvalues;
        const std::vector<node> next;
        const std::vector<double> (node::*get)(const double*) const;

        node(std::vector<double> values_, std::vector<node> next_) : 
            values(values_), 
            dvalues(diff(values_)), 
            next(next_), 
            get((next_.size() < 2 ? &node::get_last : &node::get_intermediate))
        {}

        std::vector<double> get_last(const double* x) const {
            return values;
        }

        std::vector<double> get_intermediate(const double* x) const {
            auto x_level = std::clamp(*x, values.front(), values.back());
            unsigned i = 1;
            while(values[i] < x_level) i++;

            auto f = (x_level - values[i - 1])*dvalues[i - 1];

            x++;
            auto y = next[i - 1].get(x);
            const auto y1 = next[i].get(x);

            for(unsigned j = 0; j < y1.size(); j++) {
                y[j] += f*(y1[j] - y[j]);
            }

            return y;
        }
        
    };

    node _head;

public:

    std::vector<double> get(const std::vector<double>& x) const {
        return _head.get(&x[0]);
    }
};


