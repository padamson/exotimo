#ifndef MACROS_H_
#define MACROS_H_

#define USE( i, j, k, l, val ) {
  irs = m*(j-1) + i;
  itu = m*(l-1) + k;
  iru = m*(l-1) + i;
  its = m*(j-1) + k;
  G(irs) = G(irs) + 2.0 * R(itu) * val;
  G(iru) = G(iru) - R(its) * val;
}

#endif
