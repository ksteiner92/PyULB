//
// Created by klaus on 06.01.19.
//

#ifndef LBM_UTILS_H
#define LBM_UTILS_H

#include <tuple>
#include <vector>
#include <algorithm>

#include "mesh.h"

typedef unsigned long long int ullong;

template<class T>
inline T sqr(T x)
{
   return x * x;
}

template<class T>
typename std::vector<T>::iterator insertSorted(std::vector<T>& vec, const T& item)
{
   return vec.insert(std::upper_bound(vec.begin(), vec.end(), item), item);
}

template<class ForwardIt, class T, class Compare=std::less<>>
ForwardIt binary_find(ForwardIt first, ForwardIt last, const T& value, Compare comp={})
{
   // Note: BOTH type T and the type after ForwardIt is dereferenced
   // must be implicitly convertible to BOTH Type1 and Type2, used in Compare.
   // This is stricter than lower_bound requirement (see above)

   first = std::lower_bound(first, last, value, comp);
   return first != last && !comp(value, *first) ? first : last;
}

struct RandomHash
{
   static inline ullong int64(std::size_t u)
   {
      ullong v = u * 3935559000370003845LL + 2691343689449507681LL;
      v ^= v >> 21;
      v ^= v << 37;
      v ^= v >> 4;
      v *= 4768777513237032717LL;
      v ^= v << 20;
      v ^= v >> 41;
      v ^= v << 5;
      return v;
   }

   static inline double doub(std::size_t u)
   {
      return 5.42101086242752217e-20 * int64(u);
   }
};

template<class T>
struct Hash
{
   virtual ullong operator()(const T& key) const noexcept = 0;
};

struct NullHash : public Hash<std::size_t>
{
   inline ullong operator()(const std::size_t& key) const noexcept override
   {
      return key;
   }
};

struct LineHash : public Hash<std::tuple<std::size_t, std::size_t>>
{
   inline ullong operator()(const std::tuple<std::size_t, std::size_t>& key) const noexcept override
   {
      return RandomHash::int64(std::get<0>(key)) - RandomHash::int64(std::get<1>(key));
   }
};

struct TriangleHash : public Hash<std::tuple<std::size_t, std::size_t, std::size_t>>
{
   inline ullong operator()(const std::tuple<std::size_t, std::size_t, std::size_t>& key) const noexcept override
   {
      return (RandomHash::int64(std::get<0>(key))
              ^ RandomHash::int64(std::get<1>(key))
              ^ RandomHash::int64(std::get<2>(key)));
   }
};

#endif //LBM_UTILS_H
