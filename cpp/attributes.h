//
// Created by klaus on 19.01.19.
//

#ifndef LBM_ATTRIBUTES_H
#define LBM_ATTRIBUTES_H

#include <iostream>
#include <vector>

class BaseAttribute
{
public:
   BaseAttribute(std::size_t size) : size(size) {}

   std::size_t getSize() const
   {
      return size;
   }

protected:
   std::size_t size;
};

template<class T>
class Attribute : public BaseAttribute
{
public:
   Attribute(std::size_t n) : BaseAttribute(n), values(n) {}

   Attribute(std::size_t n, const T& def) : BaseAttribute(n), values(n, def) {}

   void setValue(std::size_t idx, const T &value)
   {
      if (idx >= values.size()) {
         std::cout << "Index out of bounds" << std::endl;
         std::cout.flush();
         throw std::out_of_range("Index out of range");
      }
      values[idx] = value;
   }

   const T& getValue(std::size_t idx) const
   {
      if (idx >= values.size())
         throw std::out_of_range("Index out of range");
      return values[idx];
   }

   T& getValue(std::size_t idx)
   {
      if (idx >= values.size())
         throw std::out_of_range("Index out of range");
      return values[idx];
   }

private:
   std::vector<T> values;
};


#endif //LBM_ATTRIBUTES_H
