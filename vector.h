#ifndef VECTOR_H
#define VECTOR_H

#include <initializer_list>
#include <cmath>
#include <cstdlib>
#include <assert.h>
#include <ostream>

template <size_t size>

class Vector
{
public:
    Vector()
    {
        static_assert(size > 0, "Wrong size\n");
    }

    Vector(std::initializer_list<double> list)
    {
        static_assert(size > 0, "Wrong size\n");

        for (size_t i = 0; i < size && i < list.size(); i++)
            dataArray[i] = *(list.begin() + i);
    }

    Vector& operator+= (const Vector& b) 
    {
        for (size_t i = 0; i < size; i++)
            dataArray[i] += b.dataArray[i]; 
        return (*this);
    }

    Vector& operator-= (const Vector& b) 
    {
        for (size_t i = 0; i < size; i++)
            dataArray[i] -= b.dataArray[i];  
        return (*this);
    }

    Vector& operator*= (double b) 
    {
        for (size_t i = 0; i < size; i++)
            dataArray[i] *= b; 
        return (*this);
    }

    Vector operator+ (const Vector& b) const
    {
        Vector a1 = (*this);
        for (size_t i = 0; i < size; i++)
            a1.dataArray[i] += b.dataArray[i]; 

        return a1;
    }

    Vector operator- (const Vector& b) const
    {
        Vector a1 = (*this);
        for (size_t i = 0; i < size; i++)
            a1.dataArray[i] -= b.dataArray[i]; 

        return a1;
    }

    Vector operator- () const
    {
        Vector a1 = (*this);
        for (size_t i = 0; i < size; i++)
            a1.dataArray[i] = - a1.dataArray[i]; 

        return a1;
    }

    Vector operator* (double b) const
    {
        Vector a1 = (*this);
        for (size_t i = 0; i < size; i++)
            a1.dataArray[i] *= b; 

        return a1;
    }

    friend Vector operator*(double b, const Vector& a)
    {
        Vector a1 = a;
        for (size_t i = 0; i < size; i++)
            a1.dataArray[i] *= b; 

        return a1;
    }

    Vector operator/ (double b) const
    {
        Vector a1 = (*this);
        for (size_t i = 0; i < size; i++)
            a1.dataArray[i] /= b; 

        return a1;
    }

    double& operator[] (size_t ind)
    {
        if (ind >= size)
            throw ("Wrong index");

        return dataArray[ind];
    }

    double manhattanDistance()
    {
        double norm = 0;
        for (size_t i = 0; i < size; i++)
            norm += fabs(dataArray[i]);

        return norm;
    }

    double euclidDistance()
    {
        double norm = 0;
        for (size_t i = 0; i < size; i++)
            norm += dataArray[i] * dataArray[i];

        return sqrt(norm);
    }

    friend std::ostream& operator<< (std::ostream& os, const Vector& v)
    {
        for (size_t i = 0; i < size - 1; i++)
            os << v.dataArray[i] << " ";

        os << v.dataArray[size - 1] << std::endl;

        return os;
    } 

private:
    double dataArray[size];
};

#endif //VECTOR_H
