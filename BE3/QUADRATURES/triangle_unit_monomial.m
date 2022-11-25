function value = triangle_unit_monomial ( ex, ey )

%*****************************************************************************80
%
%% triangle_unit_monomial integrates a monomial over the unit triangle.
%
%  Discussion:
%
%    This routine integrates a monomial of the form
%
%      x^ex y^ey
%
%    where the exponents are nonnegative integers.  Note that
%    if the combination 0^0 is encountered, it should be treated
%    as 1.
%
%    Integral ( over unit triangle ) x^m y^n dx dy = m% * n% / ( m + n + 2 )%
%
%    The integration region is:
%
%      0 <= X
%      0 <= Y
%      X + Y <= 1.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    16 April 2019
%
%  Author:
%
%    John Burkardt
%
%  Parameters:
%
%    Input, integer EX, EY, the exponents.
%
%    Output, real VALUE, the integral of the monomial.
%
  value = 1.0;

  k = ex;

  for i = 1 : ey
    k = k + 1;
    value = value * i / k;
  end

  k = k + 1;
  value = value / k;

  k = k + 1;
  value = value / k;

  return
end
