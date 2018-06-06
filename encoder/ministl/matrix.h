/*
C++ provides many container classes via the Standard Template Library (STL). The most appropriate choices here are std::valarray and std::vector. The std::vector template class is much more frequently utilised because it is more general. std::vector can handle many different types (including pointers and smart pointer objects), whereas std::valarray is designed solely for numerical values and thus the compiler can make certain optimisations.

At first glance std::valarray would seem like a great choice to provide storage for our matrix values. However, in reality it turns out that compiler support for the numerical optimisations that std::valarray is supposed to provide does not really exist. Not only that but the std::valarray has a poorer application programming interface (API) and isn't as flexible as a std::vector. Thus it makes sense to use the std::vector template class for our underlying storage mechanism.

Supposing that our matrix has M
rows and N columns, we could either create a single std::vector of length N×M or create a "vector of vectors". The latter creates a single vector of length M, which takes a std::vector<T> of types as its type. The inner vectors will each be of length N

. We will utilise the "vector of vectors" approach. The primary reason to use such a mechanism is that we gain a good API, helpful in accessing elements of such a vector "for free". We do not need to look up the correct element as we would need to do in the single large vector approach. The declaration for this type of storage mechanism is given by:

std::vector<std::vector<T> >

Where T
is our type placeholder. Note: For nearly all of the quantitative work we carry out, we will use the double precision type for numerical storage.


 * auto x = range(10); // [0, ..., 9]
auto y = range(2, 20); // [2, ..., 19]
auto z = range(10, 2, -2); // [10, 8, 6, 4]
 * 
 * 
 * #include <vector>
#include <stdexcept>

template <typename IntType>
std::vector<IntType> range(IntType start, IntType stop, IntType step)
{
  if (step == IntType(0))
  {
    throw std::invalid_argument("step for range must be non-zero");
  }

  std::vector<IntType> result;
  IntType i = start;
  while ((step > 0) ? (i < stop) : (i > stop))
  {
    result.push_back(i);
    i += step;
  }

  return result;
}

template <typename IntType>
std::vector<IntType> range(IntType start, IntType stop)
{
  return range(start, stop, IntType(1));
}

template <typename IntType>
std::vector<IntType> range(IntType stop)
{
  return range(IntType(0), stop, IntType(1));
}
 * 
 *  */

#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <iostream>

template <typename T> class Matrix {
 private:
  std::vector<std::vector<T> > mat;
  unsigned rows;
  unsigned cols;

 public:
  Matrix(unsigned _rows, unsigned _cols, const T& _initial);
  
  Matrix(std::initializer_list< std::vector<T> > l): mat(l),rows(0),cols(0) {
         //std::cout << "constructed with a " << l.size() << "-element list\n";

         rows = mat.size();
         
        for (const auto & tab : mat) {
            
            if( cols &&  cols != tab.size())
            {
                
                display( );
                throw "all the column should be same size";
              
            }
           cols = tab.size();   
        }
    }
  

  void display( )
  {
      // Prints welcome message...
    std::cout << "Start ..." << std::endl;

    for( const auto & tab : mat )
    {
        for( const auto & x : tab )
        {
            std::cout << x << " ";
        }
        std::cout << std::endl;
    }
    
    std::cout << "end ..." << std::endl;
  }
  
  Matrix(const Matrix<T>& rhs);
  virtual ~Matrix();

  // Operator overloading, for "standard" mathematical matrix operations                                                                                                                                                          
  Matrix<T>& operator=(const Matrix<T>& rhs);

  // Matrix mathematical operations                                                                                                                                                                                               
  Matrix<T> operator+(const Matrix<T>& rhs);
  Matrix<T>& operator+=(const Matrix<T>& rhs);
  Matrix<T> operator-(const Matrix<T>& rhs);
  Matrix<T>& operator-=(const Matrix<T>& rhs);
  Matrix<T> operator*(const Matrix<T>& rhs);
  Matrix<T>& operator*=(const Matrix<T>& rhs);
  Matrix<T> transpose();

  // Matrix/scalar operations                                                                                                                                                                                                     
  Matrix<T> operator+(const T& rhs);
  Matrix<T> operator-(const T& rhs);
  Matrix<T> operator*(const T& rhs);
  Matrix<T> operator/(const T& rhs);

  // Matrix/vector operations                                                                                                                                                                                                     
  std::vector<T> operator*(const std::vector<T>& rhs);
  std::vector<T> diag_vec();

  // Access the individual elements                                                                                                                                                                                               
  T& operator()(const unsigned& row, const unsigned& col);
  const T& operator()(const unsigned& row, const unsigned& col) const;
  
   Matrix<T> slice(int start, int end)
   {
       std::vector<std::vector<T> > mat1(mat.begin() , mat.begin()+1); 
   }
 
  // Matrix<T>  operator[](const unsigned& row)
  // {
   //    return *this;
  // }
   
    
  // Access the row and column sizes                                                                                                                                                                                              
  unsigned get_rows() const;
  unsigned get_cols() const;

};

#include "matrix.cpp"

#endif