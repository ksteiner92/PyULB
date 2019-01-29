//
// Created by klaus on 06.01.19.
//

#ifndef LBM_UTILS_H
#define LBM_UTILS_H

#include <tuple>

#include "mesh.h"

template<class T>
inline T sqr(T x)
{
   return x * x;
}

struct RandomHash
{
   typedef unsigned long long int Ullong;
   static inline Ullong int64(Ullong u)
   {
      Ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
      v ^= v >> 21;
      v ^= v << 37;
      v ^= v >> 4;
      v *= 4768777513237032717LL;
      v ^= v << 20;
      v ^= v >> 41;
      v ^= v << 5;
      return v;
   }

   static inline double doub(Ullong u)
   {
      return 5.42101086242752217e-20 * int64(u);
   }
};

template<class T>
struct Hash
{
   virtual unsigned long long operator()(const T& key) const noexcept = 0;
};

struct NullHash : public Hash<unsigned long long int>
{
   typedef unsigned long long int Ullong;

   inline Ullong operator()(const Ullong& key) const noexcept override
   {
      return key;
   }
};

struct LineHash : public Hash<std::tuple<int, int>>
{
   typedef unsigned long long int Ullong;

   inline Ullong operator()(const std::tuple<int, int>& key) const noexcept override
   {
      return RandomHash::int64(std::get<0>(key)) - RandomHash::int64(std::get<1>(key));
   }
};

struct TriangleHash : public Hash<std::tuple<int, int, int>>
{
   typedef unsigned long long int Ullong;

   inline Ullong operator()(const std::tuple<int, int, int>& key) const noexcept override
   {
      return (RandomHash::int64(std::get<0>(key))
              ^ RandomHash::int64(std::get<1>(key))
              ^ RandomHash::int64(std::get<2>(key)));
   }
};

#endif //LBM_UTILS_H
