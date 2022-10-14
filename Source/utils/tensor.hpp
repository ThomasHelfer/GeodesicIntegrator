#ifndef TENSOR_HPP_
#define TENSOR_HPP_


#include "DimensionDefinitions.hpp"

#define IDX_SPACEDIM 4

/// This class implements a tensor with given rank, element data type, and
/// dimension.  By default the dimension is equal to IDX_SPACEDIM.
template <int rank, class data_t, int size = IDX_SPACEDIM> class tensor
{
    template <int, class, int> friend class tensor;
    typedef typename tensor<rank - 1, data_t, size>::arr_t arr_t[size];
    arr_t arr;

  public:
    tensor() {}

    //    
    //    tensor(std::initializer_list<data_t> list) :
    //        arr (list)
    //    {}

    template <typename... T>  tensor(T... data) : arr{data...} {}

    operator arr_t &() { return arr; }

    operator const arr_t &() const { return arr; }
};

template <class data_t, int size> class tensor<0, data_t, size>
{
    template <int, class, int> friend class tensor;
    typedef data_t arr_t;
    arr_t arr;

  public:
    tensor() {}

    
    tensor(data_t val) : arr(val) {}

    operator arr_t &() { return arr; }

    operator const arr_t &() const { return arr; }
};

#endif /* TENSOR_HPP_ */