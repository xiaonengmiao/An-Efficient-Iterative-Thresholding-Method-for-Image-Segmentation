function [u_hat1,u_hat2,u_hat3,u_hat4,u_hat5,u_hat6]=HeatConv(dt,u1,u2,u3,u4,u5,u6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fast evaluation of convelution using FFT by the basic principle
% g*h=sum(\frac{1}{dx*dy}iFFT(FFT(g).FFT(h)))
% We may do some extension to make the image to be perodic in both x, y
% direction, That depends on the profile of the image and can be modified
% if needed. In this simple code, we just assume the nonzero part are away
% from the boundary of computational domain.
% Input: dt --- artificial time step
%        u  --- characteristic function of different regions.
% Output:u_hat --- diffused value of the charecteristic function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% size of the image
[M,N] = size(u1);
x1 = (0:M-1)*2*pi/M-pi;
x2 = (0:N-1)*2*pi/N-pi;
dx1 = x1(2)-x1(1);
dx2 = x2(2)-x2(1);
xx1 = repmat(x1,N,1);
xx2 = repmat(x2,M,1)';
G_dt = (1/(4*pi*dt)*exp(-(xx1.^2+xx2.^2)/(4*dt)))';
K = fft2(G_dt);
% 
% u_hat1 = real(dx1*dx2*ifftshift(ifft2(fft2(u1).*fft2(G_dt))));
% u_hat2 = real(dx1*dx2*ifftshift(ifft2(fft2(u2).*fft2(G_dt))));
% if nargin ==4 
%    u_hat3 = real(dx1*dx2*ifftshift(ifft2(fft2(u3).*fft2(G_dt))));
% elseif nargin ==5 
%   u_hat4 = real(dx1*dx2*ifftshift(ifft2(fft2(u4).*fft2(G_dt))));
% elseif nargin ==6
%    u_hat5 = real(dx1*dx2*ifftshift(ifft2(fft2(u5).*fft2(G_dt))));
% elseif nargin ==7 
%    u_hat6 = real(dx1*dx2*ifftshift(ifft2(fft2(u6).*fft2(G_dt))));
% end

u_hat1 = real(dx1*dx2*ifftshift(ifft2(fft2(u1).*K)));
u_hat2 = real(dx1*dx2*ifftshift(ifft2(fft2(u2).*K)));
if nargin ==4 
   u_hat3 = real(dx1*dx2*ifftshift(ifft2(fft2(u3).*K)));
elseif nargin ==5 
    u_hat3 = real(dx1*dx2*ifftshift(ifft2(fft2(u3).*K)));
    u_hat4 = real(dx1*dx2*ifftshift(ifft2(fft2(u4).*K)));
elseif nargin ==6
    u_hat3 = real(dx1*dx2*ifftshift(ifft2(fft2(u3).*K)));
    u_hat4 = real(dx1*dx2*ifftshift(ifft2(fft2(u4).*K)));
   u_hat5 = real(dx1*dx2*ifftshift(ifft2(fft2(u5).*K)));
elseif nargin ==7 
    u_hat3 = real(dx1*dx2*ifftshift(ifft2(fft2(u3).*K)));
    u_hat4 = real(dx1*dx2*ifftshift(ifft2(fft2(u4).*K)));
   u_hat5 = real(dx1*dx2*ifftshift(ifft2(fft2(u5).*K)));
   u_hat6 = real(dx1*dx2*ifftshift(ifft2(fft2(u6).*K)));
end


end