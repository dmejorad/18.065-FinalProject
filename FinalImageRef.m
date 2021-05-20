clear
ref = imread('Middle_orig-page-001.jpg');
target = imread('Right_orig-page-001.jpg');
refSize = size(ref);
targSize = size(target);
ref = double(ref(refSize(1)/2:refSize(1)/2+2000,1000:5000,1));
target = double(target(refSize(1)/2:refSize(1)/2+2000,1000:5000,1));
% target = targetog(200:1000,1000:2000);
% figure(2); imagesc(targetog); colormap gray

% ref = double(imread('Test1Ref.png'));
% target = double(imread('Test1RefTran.png'));
% refSize = size(ref);
% targSize = size(target);
% ref = ref + (2.*rand(refSize)+1);
% target = target + (2.*rand(targSize)+1);
% ref = ref(34:684,100:750);
% target = target(44:694,120:770);


figure(2); imagesc(target); colormap gray
figure(1); imagesc(ref); colormap gray
% figure(3); imagesc(target); colormap gray


[FX,FY] = gradient(target); %Precalculate gradient of target image
coord_ref = imref2d(size(ref));



tx = 750;
ty = 0;
Omega = [tx,ty]';
M = [1 0 Omega(1);...
     0 1 Omega(2);...
     0 0 1 ];

tform = affine2d(M');

[init] = imwarp(target,tform,'OutputView',coord_ref);
% [FXT] = imwarp(FX,tform,'OutputView',coord_ref);
% [FYT] = imwarp(FY,tform,'OutputView',coord_ref);
figure(3); imagesc(init); colormap gray




numofiter = 500; % Define how long the minimization should run
eta = 0.1; % Learning rate 
etaog = eta;
etaadapt = true;
drawall = false
%%%%%%%%%%%%%%%%% Initialize the starting parameters
% tx = 0;
% ty = 0;
Omega = [tx,ty]';
M = [1 0 Omega(1);...
     0 1 Omega(2);...
     0 0 1 ];
tform = affine2d(M');
[init] = imwarp(target,tform,'OutputView',coord_ref);

initdel = ref -init;
error = zeros(numofiter,1);
% clf(11)
figure(11)
h = animatedline;
xlim([0 numofiter+1])
bad = false;
badcum = 0;
b = ProgressBar(numofiter, ...
        'UpdateRate', inf, ...
        'IsActive', true);
stepsize = zeros(numofiter,1);
for ii = 1:numofiter
    b(1, [], []);
    M = [1 0 Omega(1);...
         0 1 Omega(2);...
         0 0 1];
    tform = affine2d(M');
    [trans] = imwarp(target,tform,'OutputView',coord_ref);
    [FXT] = imwarp(FX,tform,'OutputView',coord_ref);
    [FYT] = imwarp(FY,tform,'OutputView',coord_ref);
    
    delta = ref-trans;
    gx = sum(2 .*delta(:).* FXT(:));
    gy = sum(2 .*delta(:).*FYT(:));
    error(ii) = sum(delta(:).^2);
    Grad = [gx gy]';
    Grad = Grad/norm(Grad);
    
    perchange = error(ii)/error(1);
    if etaadapt

        if ii>2 & bad == false
           if error(ii)>error(ii-1)
               Grad = [0,0]';
               eta = eta * (perchange);
               bad = true;
               badcum = badcum +1;
           elseif error(ii)<error(ii-1)
               if badcum > -1
                   Omega = Omega - 1.5*eta.*Grad;
               end

           end
        end
    end
    
    
    Omega = Omega - eta.*Grad;
    Omega(2) = ty;
    stepsize(ii) = eta;
    if bad ==true
        eta = etaog;
        bad = false;
    end
    
%     figure(11)
%     addpoints(h,ii,error(ii)/error(1));
%     drawnow
%     title('Error')
% 
% 
    if drawall == true
        figure(4); imagesc(delta); colormap gray
        title(['Final Delta iter#:',num2str(ii)])
        drawnow 
        figure(11)
        addpoints(h,ii,error(ii)/error(1));
        drawnow
        title('Error')

    end
    
end
b.release();

figure(12)
hold on
plot(1:numofiter,error/max(error)); 
title('Error')
xlim([0 numofiter+1])
% ylim([0,1])
if drawall == false
    figure(4); imagesc(delta); colormap gray
    title('Delta')
end

figure(7)
subplot(1,3,1)
imagesc(initdel); colormap gray
title('Initial Difference')
subplot(1,3,2)
imagesc(delta); colormap gray
title('Final Difference')
subplot(1,3,3)
plot(1:numofiter,error/max(error)); 
title('Error vs Iterations')
ylabel('Norm(Error)')
xlabel('Iteration #')

figure(14); imagesc(ref+init); colormap gray
title('initial')

% M = [1 0 Omega(1)-1000;...
%      0 1 Omega(2)-200;...
%      0 0 1 ];
% 
% tform = affine2d(M');
% [final] = imwarp(target,tform,'OutputView',coord_ref);
figure(13); imagesc(trans+ref); colormap gray
title('Combined Image')
%%
save(sprintf('Newtx%dty%dNumIter%deta%d.mat',tx, ty,numofiter,eta))

