function v = monomial_value_2d ( n, ex, ey, x, y )

%*****************************************************************************80
%
%% monomial_value_2d evaluates a monomial in x and y.
%
%  Discussion:
%
%    This routine evaluates a monomial of the form
%
%      x^ex * y^ey
%
%    The combination 0.0^0 is encountered is treated as 1.0.
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
%    Input, integer N, the number of evaluation points.
%
%    Input, integer EX, EY, the exponents.
%
%    Input, real X(N), Y(N), the point coordinates.
%
%    Output, real V(N), the monomial values.
%
  v = ones ( n, 1 );

  if ( 0 ~= ex )
    v(1:n) = v(1:n) .* x(1:n) .^ ex;
  end

  if ( 0 ~= ey )
    v(1:n) = v(1:n) .* y(1:n) .^ ey;
  end

  return
end
