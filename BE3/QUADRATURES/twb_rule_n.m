function n = twb_rule_n ( strength )

%*****************************************************************************80
%
%% twb_rule_n returns the order of a TWB rule of given strength.
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
%  Reference:
%
%    Mark Taylor, Beth Wingate, Len Bos, 
%    Several new quadrature formulas for polynomial integration in the triangle, 
%    http://arxiv.org/abs/math/0501496v2,
%    08 February 2007.
%
%  Parameters:
%
%    Input, integer STRENGTH, the desired strength.
%    1 <= STRENGTH.
%
%    Output, integer N, the order of the rule, if it exists.
%    Otherwise, the value -1 is returned.
%
  if ( strength == 1 )
    n = 1;
  elseif ( strength == 2 )
    n = 3;
  elseif ( strength == 4 )
    n = 6;
  elseif ( strength == 5 )
    n = 10;
  elseif ( strength == 7 )
    n = 15;
  elseif ( strength == 9 )
    n = 21;
  elseif ( strength == 11 )
    n = 28;
  elseif ( strength == 13 )
    n = 36;
  elseif ( strength == 14 )
    n = 45;
  elseif ( strength == 16 )
    n = 55;
  elseif ( strength == 18 )
    n = 66;
  elseif ( strength == 20 )
    n = 78;
  elseif ( strength == 21 )
    n = 91;
  elseif ( strength == 23 )
    n = 105;
  elseif ( strength == 25 )
    n = 120;
  else
    n = -1;
  end

  return
end
