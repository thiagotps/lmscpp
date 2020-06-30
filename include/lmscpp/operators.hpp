#ifndef __OPERATORS_H_
#define __OPERATORS_H_

#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>

namespace SymEngine
{
  namespace OverloadedOperators
  {
    RCP<const Basic>
    operator+(const RCP<const Basic> &l, const RCP<const Basic> &r)
    {
      return add(l,r);
    }

    RCP<const Basic>
    operator-(const RCP<const Basic> &l, const RCP<const Basic> &r)
    {
      return sub(l,r);
    }

    RCP<const Basic>
    operator*(const RCP<const Basic> &l, const RCP<const Basic> &r)
    {
      return mul(l,r);
    }

  }
}

#endif // __OPERATORS_H_
