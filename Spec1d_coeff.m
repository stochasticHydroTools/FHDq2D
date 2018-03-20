function [specGR_record, specGG_record, specRR_record,specGR_record_half, specGG_record_half, specRR_record_half] = ...
      Spec1d_coeff(c, dx, dy, iter, NX, NY, spectrum1d_gap, ...
      specGR_record, specGG_record, specRR_record, specGR_record_half, specGG_record_half, specRR_record_half)
% Computes S(k) statistics as HydroGrid does
spec2d_RR = abs(real2fs2d(c(:,:,2),dx,dy)).^2;
spec2d_GG = abs(real2fs2d(c(:,:,1),dx,dy)).^2;
spec2d_GR = real(real2fs2d(c(:,:,2),dx,dy).*conj(real2fs2d(c(:,:,1),dx,dy)));
specGG_record(:,iter/spectrum1d_gap+1) = specGG_record(:,iter/spectrum1d_gap+1) + spec2d_GG(1,2:NX/2)';
specRR_record(:,iter/spectrum1d_gap+1) = specRR_record(:,iter/spectrum1d_gap+1) + spec2d_RR(1,2:NX/2)';
specGR_record(:,iter/spectrum1d_gap+1) = specGR_record(:,iter/spectrum1d_gap+1) + spec2d_GR(1,2:NX/2)';
specGR_record_half = specGR_record_half + spec2d_GR(1,2:NX/2)';
specGG_record_half = specGG_record_half + spec2d_GG(1,2:NX/2)';
specRR_record_half = specRR_record_half + spec2d_RR(1,2:NX/2)';

end

