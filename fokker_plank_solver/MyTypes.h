#ifndef MYTYPES_H_INCLUDED
#define MYTYPES_H_INCLUDED

class Array3d{
public:
    Array3d(unsigned size1, unsigned size2,unsigned size3);
    Array3d(unsigned size1, unsigned size2,unsigned size3,double * arr_in);
    double& operator() (unsigned i1,unsigned i2, unsigned i3);
    double operator() (unsigned i1,unsigned i2, unsigned i3) const;
    ~Array3d(){delete[] data_;};
    double* getp_data() const{return data_;};
private:
    unsigned size1_, size2_, size3_;
    double* data_;
};

inline
Array3d::Array3d(unsigned size1, unsigned size2,unsigned size3,double * arr_in)
    : size1_ (size1)
    , size2_ (size2)
    , size3_ (size3)
{
    data_ = arr_in;
}

inline
Array3d::Array3d(unsigned size1, unsigned size2,unsigned size3)
    : size1_ (size1)
    , size2_ (size2)
    , size3_ (size3)
{
    data_ = new double[size1*size2*size3];
}

inline
double& Array3d::operator() (unsigned i1,unsigned i2, unsigned i3){
  return data_[i3+size3_*i2+size3_*size2_*i1];
}

inline
double Array3d::operator() (unsigned i1,unsigned i2, unsigned i3) const{
  return data_[i3+size3_*i2+size3_*size2_*i1];
};

class Array1d{
public:
    Array1d(unsigned size1);
    Array1d(unsigned size1,double* arr_in);
    double& operator() (unsigned i1);
    double operator() (unsigned i1) const;
    ~Array1d(){delete[] data_;};
    double* getp_data() const{return data_;};
    long get_size1()const{return  size1_;}
private:
    unsigned size1_;
    double* data_;
};
inline
Array1d::Array1d(unsigned size1, double* arr_in)
    : size1_ (size1)
{
    data_ = arr_in;
}
inline
Array1d::Array1d(unsigned size1)
    : size1_ (size1)
{
    data_ = new double[size1];
}

inline
double& Array1d::operator() (unsigned i){
  return data_[i];
}
inline
double Array1d::operator() (unsigned i) const{
  return data_[i];
};

#endif


