function [  ] = plot_with_errorbars( x, y, dy )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here



hE     = errorbar(x, y, dy);

set(hE                            , ...
  'LineStyle'       , 'none'      , ...
  'Marker'          , '.'         , ...
  'Color'           , [.3 .3 .3]  );
end

