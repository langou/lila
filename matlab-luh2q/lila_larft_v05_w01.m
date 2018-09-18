   function [ T ] = lila_larft_v05_w01( m, n, i, mt, A, T )
%
      ilo = i;
      ihi = ilo + n - 1;
%
      for ii = 1:n
%
         T( ilo:ilo,ilo:ilo ) = lapack_larft( A( ilo:ilo,ilo:ilo ) );
% 
      end
%
   end
