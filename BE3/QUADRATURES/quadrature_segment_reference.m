function [Ng1D,xg1D,wg1D] = quadrature_segment_reference(ordre)

Ng1D = ordre+1;
[xg1D,wg1D] = gauss1D(Ng1D);
xg1D = (xg1D+1)/2;
wg1D = wg1D/2;

end 