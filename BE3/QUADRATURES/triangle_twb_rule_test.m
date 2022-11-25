function triangle_twb_rule_test ( )

%*****************************************************************************80
%
%% triangle_twb_rule_test tests triangle_twb_rule.
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
  addpath ( '../triangle_twb_rule' );

  timestamp ( );
  fprintf ( 1, '\n' );
  fprintf ( 1, 'triangle_twb_rule_test\n' );
  fprintf ( 1, '  MATLAB version\n' );
  fprintf ( 1, '  Test triangle_twb_rule.\n' );

  degree_max = 5;
  triangle_unit_quad_test ( degree_max );
%
%  Terminate.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'triangle_twb_rule_test\n' );
  fprintf ( 1, '  Normal end of execution.\n' );
  fprintf ( 1, '\n' );
  timestamp ( );

  rmpath ( '../triangle_twb_rule' );

  return
end
function triangle_unit_quad_test ( degree_max )

%*****************************************************************************80
%
%% triangle_unit_quad_test tests rules for the unit triangle.
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
%    Input, int DEGREE_MAX, the maximum total degree of the monomials to check.
%
  fprintf ( 1, '\n' );
  fprintf ( 1, 'triangle_unit_quad_test\n' );
  fprintf ( 1, '  Approximate monomial integrals in triangle with TWB rules.\n' );

  degree = 0;
  ex = 0;
  ey = degree;

  while ( true )

    fprintf ( 1, '\n' );
    fprintf ( 1, '  Monomial: x^%d y^%d\n', ex, ey );

    for strength = 1 : 25

      n = twb_rule_n ( strength );

      if ( n < 1 )
        continue;
      end

      w = twb_rule_w ( strength );
      x = twb_rule_x ( strength );
      y = twb_rule_y ( strength );
      v = monomial_value_2d ( n, ex, ey, x, y );
      q = w' * v;
      fprintf ( 1, '  %6d  %6d  %14.6g\n', strength, n, q );

    end

    q = triangle_unit_monomial ( ex, ey );
    fprintf ( 1, '   Exact          %14.6g\n', q );

    if ( ex < degree )
      ex = ex + 1;
      ey = ey - 1;
    elseif ( degree < degree_max )
      degree = degree + 1;
      ex = 0;
      ey = degree;
    else
      break;
    end

  end

  return
end
