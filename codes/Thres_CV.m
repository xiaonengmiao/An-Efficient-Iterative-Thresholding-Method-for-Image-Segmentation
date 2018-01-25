%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   
%  This code implement An Efficient Iterative Thresholding 
%  Method for Image Segmentation
%-----------------------------------------------------------------
%  usage of variables:
%        input: 
%                I: input image, method: choose from provided
%                such as chan, vector, multiphase, multiphasevector
%        output: 
%                segment of I
%
%-----------------------------------------------------------------
%  Created by lihaohan on 08/10/2016.
%  Copyright @ 2016 HKUST. All rights reserved.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%


function [u1,u2]=Thres_CV(I,method)
% -- resize original image
%    s = 512./min(size(I,1),size(I,2));
%    if s<1
%        I = imresize(I,s);
%    end
   M = size(I,1);
   N = size(I,2);
%%
% -- choose method  
switch lower(method)
    case 'chan'
        if size(I,3) == 3
            P = rgb2gray(uint8(I));
            P = double(P);
        elseif size(I,3) == 2
            P = 0.5.*(double(I(:,:,1))+double(I(:,:,2)));
        else
            P = double(I);
        end
        P = P./max(max(abs(P)));
    case 'vector'
        P = double(I);
        for i = 1:3
            P(:,:,i) = P(:,:,i)./max(max(abs(P(:,:,i))));
        end
    case 'multiphase'
        if size(I,3) == 3
            P = rgb2gray(uint8(I));
            P = double(P);
        elseif size(I,3) == 2
            P = 0.5.*(double(I(:,:,1))+double(I(:,:,2)));
        else
            P = double(I);
        end
        P = P./max(max(abs(P)));
    case 'multiphasevector'
        P = double(I);
        for i = 1:3
            P(:,:,i) = P(:,:,i)./max(max(abs(P(:,:,i))));
        end
    otherwise
        error('!invalid method')
end

%%
% -- looping
switch lower(method)


    case 'chan'
        %% Standard ChanVese 
        dt = 0.03;  % time step
        alpha = 0.01; % parameter alpha
        lamda = alpha * sqrt(pi)/sqrt(dt);
        
        % -- initial contour
        u1 = zeros(M,N);
        S = 10;
        u1(S:M-S,S:N-S) = ones(M-2*S+1,N-2*S+1);
        u2 = ones(M,N)-u1;
        
        % -- plot images
        figure();
        subplot(2,2,1); imshow(I); title('Input Image');
        subplot(2,2,2); imshow(P); title('Noised Image');
        subplot(2,2,3); imshow(P); hold; contour(u1, [0.5 0.5], 'r','LineWidth',1); title('initial contour');
        subplot(2,2,4); title('Segmentation');
        
        u1p = - ones(size(u1)); % previous u1
        u2p = - ones(size(u1)); % previous u2
        change = 1;        
        tic;        
        k = 0;
        while change > 1e-12
            change = norm(u1-u1p) + norm(u2-u2p);
            u1p = u1;
            u2p = u2;
            k = k+1;
            [f1,f2] = daterm(P,u1,u2);
            [uh1,uh2] = HeatConv(dt,u1,u2); % convolution with heat kernel
            index1 = f1+2*lamda*uh2;
            index2 = f2+2*lamda*uh1;
            u1 = double(index1<=index2); % thresholding if u1>u2 then u1=1 vise versa
            u2 = 1-u1;

            imshow(P);
            hold on
            contour(u1, [0.5 0.5], 'r','LineWidth',2);
            hold off
            axis off 
            axis square
            title([num2str(k) ' Iterations']);
            drawnow
%             pause(0.5);
        end
        toc;

        figure();
        imshow(P);
        hold on
        contour(u1, [0.5 0.5], 'r','LineWidth',1.3);
        hold off
        title([num2str(k) ' Iterations']);
    
    case 'vector'
        %% color image approach 
        f1 = P(:,:,1); % first color channel
        f2 = P(:,:,2); % second
        f3 = P(:,:,3); % third
        S = 50;
        k=0;
        dt = 0.0002;  % time step
        
        % -- initial contour
        lamda = 0.09;
        u1 = zeros(M,N);
        u1(S:M-S,S:N-S) = ones(M-2*S+1,N-2*S+1);
        u2 = ones(M,N)-u1;
        
        % -- plot images
        figure();
        subplot(2,2,1); imshow(I); title('Input Image');
        subplot(2,2,2); imshow(P); hold; contour(u1, [0.5 0.5], 'k','LineWidth',1); title('initial contour');
        subplot(2,2,3); title('Segmentation');

        u1p = - ones(size(u1)); % previous u1
        u2p = - ones(size(u1)); % previous u2
        change = 1;

        tic;
        while change > 1e-12
        
        change = norm(u1-u1p) + norm(u2-u2p);
        u1p = u1;
        u2p = u2;
        k = k+1;
        [f11,f21] = daterm(f1,u1,u2); % data term of first color channel
        [f12,f22] = daterm(f2,u1,u2); % data term of second color channel
        [f13,f23] = daterm(f3,u1,u2); % data term of third color channel
        [uh1,uh2] = HeatConv(dt,u1,u2); % heat kernel convolution

        index1 = f11+f12+f13-6*lamda*uh1;
        index2 = f21+f22+f23-6*lamda*uh2;
        u1 = double(index1<=index2); % thresholding if u1>u2 then u1=1 vise versa
        u2 = 1-u1;            
        imshow(I);
        hold on
        contour(u1, [0.5 0.5], 'k','LineWidth',2);
        hold off
        title([num2str(k) ' Iterations']);
        drawnow;
%         pause(0.02)
        end
        toc;

        figure()
        imshow(I);
        hold on
        contour(u1, [0.5 0.5], 'k','LineWidth',1.3);
        hold off
        title([num2str(k) ' Iterations']);
        figure(); imshow(u1); title('Global Region-Based Segmentation');
    case 'multiphase'
        %% multiphase approach
        k = 0;
        dt = 0.025;  % time step
        alpha = 0.005;
        lamda = alpha * sqrt(pi)/sqrt(dt);

        % -- initial contour
        u1 = zeros(M,N);
        u1(40:80,100:200) = 1;
        u2 = zeros(M,N);
        u2(120:140,100:200) = 1;
        u3 = 1 - u1 - u2;
        
        % -- plot images
        figure();
        subplot(2,2,1); imshow(I); title('Input Image');
        subplot(2,2,2); imshow(I); hold; contour(u1 * 0 + u2 * 0.5 + u3 * 1, [0.1 0.9], 'r','LineWidth',1); title('initial contour');
        subplot(2,2,3); title('Segmentation');
         
        u1p = - ones(size(u1)); % previous u1
        u2p = - ones(size(u1)); % previous u2
        u3p = - ones(size(u1)); % previous u3
        change = 1;

        tic;
        while change>1e-12
            change = norm(u1-u1p) + norm(u2-u2p) + norm(u3-u3p);
            u1p = u1;
            u2p = u2;
            u3p = u3;
            k = k+1;
            [f1,f2,f3] = daterm(P,u1,u2,u3); % data term
            [uh1,uh2,uh3] = HeatConv(dt,u1,u2,u3); % heat kernel convolution
            index1 = f1-2*lamda*uh1;
            index2 = f2-2*lamda*uh2;
            index3 = f3-2*lamda*uh3;
            % -- thresholding: if ui is the largest value then ui=1 and
            % uj=0 for j!=i
            u1 = double(index1<=index2).*double(index1<=index3);
            u2 = double(index2<index1).*double(index2<=index3).*double(u1==0);
            u3 = 1-u1-u2;
                       
            imshow(I);
            hold on
            contour(u1 * 0 + u2 * 0.5 + u3 * 1, 'r', 'LineWidth', 1);
            hold off
            title([num2str(k) ' Iterations']);
            drawnow

        end
        toc;
        figure(); imshow(u1 * 0 + u2 * 0.5 + u3 * 1); title('Global Region-Based Segmentation');
    case 'multiphasevector'
        %% multiphase color image approach
        f1 = P(:,:,1); % first color channel
        f2 = P(:,:,2); % second
        f3 = P(:,:,3); % third
        k = 0;
        dt = 0.01;  % time step
        alpha = 0.003;
        lamda = alpha * sqrt(pi)/sqrt(dt); % parameter lamda
        
        % -- initial contour
        u1 = zeros(M,N);
        u1(100:200,100:300) = 1;
        u2 = zeros(M,N);
        u2(220:350,100:300) = 1;
        u3 = zeros(M,N);
        u3(100:200,400:450) = 1;
        u4 = 1 - u1 - u2 - u3;
    
        % -- plot images
        figure();
        subplot(2,2,1); imshow(I); title('Input Image');
        subplot(2,2,2); imshow(I); hold; contour(u1 * 0 + u2 * 0.3 + u3 * 0.6 + u4, [0.1 0.9], 'k','LineWidth',1); title('initial contour');
        subplot(2,2,3); title('Segmentation');
        
        u1p = - ones(size(u1)); % previous u1
        u2p = - ones(size(u1)); % previous u2
        u3p = - ones(size(u1)); % previous u3
        u4p = - ones(size(u1));
        change = 1;

        tic;
        while change>1e-12
        
            change = norm(u1-u1p) + norm(u2-u2p) + norm(u3-u3p) + norm(u4-u4p);
            u1p = u1;
            u2p = u2;
            u3p = u3;
            u4p = u4;
            k = k+1;
            [f11,f21,f31,f41] = daterm(f1,u1,u2,u3,u4); % data term of first color channel
            [f12,f22,f32,f42] = daterm(f2,u1,u2,u3,u4); % data term of second color channel
            [f13,f23,f33,f43] = daterm(f3,u1,u2,u3,u4); % data term of third color channel
            [uh1,uh2,uh3,uh4] = HeatConv(dt,u1,u2,u3,u4); % heat kernel convolution
            index1 = f11+f12+f13-6*lamda*uh1;
            index2 = f21+f22+f23-6*lamda*uh2;
            index3 = f31+f32+f33-6*lamda*uh3;
            index4 = f41+f42+f43-6*lamda*uh4;
            
            % -- thresholding: if ui is the largest value then ui=1 and
            % uj=0 for j!=i
            u1 = double(index1<=index2).*double(index1<=index3).*double(index1<=index4);
            u2 = double(index2<index1).*double(index2<=index3).*double(index2<=index4).*double(u1==0);
            u3 = double(index3<index1).*double(index3<index2).*double(index3<=index4).*double(u2==0).*double(u1==0);    
            u4 = 1-u1-u2-u3;
                    
            imshow(I);
            hold on
            contour(u1 * 0 + u2 * 0.3 + u3 * 0.6+u4, 'k', 'LineWidth', 2);
            hold off
            title([num2str(k) ' Iterations']);
            drawnow
         
        end
        toc;
         
        figure()
        imshow(I); hold; contourf(u1 * 0 + u2 * 0.3 + u3 * 0.6 + u4); 
        title('Global Region-Based Segmentation');
    otherwise
        error('wrong');
end
                  

        
        
        
        
        
        
        
        
        