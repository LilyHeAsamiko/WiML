format longG

I = imread('2891.jpg');
I1 = imread('2986.jpg');
I_ = I ;
I0 = zeros(size(I));
imshow(I,[])
[m, n, k] = size(I);

for d = 1: k % according to staining 
    I_(:,:,d) = mean_filter(I(:,:,d),4);    
end
%feature extraction(features can be used for cell detection)
%kse =  Haarlikefeature(I0(:,:,1), winL);
%Maximum level to which to perform the 2-D Haar transform, specified as a positive integer. The default value depends on the length of the input signal, x.
%If both the row and column sizes of x are powers of two, the 2-D Haar transform is obtained down to level log2(min(size(x))).
%If both the row and column sizes of x are even, but at least one is not a power of two, level is equal to floor(log2(min(size(x)/2))).
level = floor(log2(min(size(I_(:,:,1))/2)));
[a,h,v,d] = haart2(I_(:,:,1),level,'integer');

nn = 1;%length(h) choose batch size
ah = size(h{nn},1);
bh = size(h{nn},2);
winL1h = uint32(floor(m/ah));
winL2h = uint32(floor(n/bh));

av = size(v{nn},1);
bv = size(v{nn},2);
winL1v = uint32(floor(m/av));
winL2v = uint32(floor(n/bv));

ad = size(d{nn},1);
bd = size(d{nn},2);
winL1d = uint32(floor(m/ad));
winL2d = uint32(floor(n/bd));
%   [Imaxx, indx] = max(I_);% can be set according to the staining color as [R,G,B]
%    [Imaxy, indy] =
%    max(Imaxx);..................................................................................................-------------------------------------------------------------------------------------------------------------------------------------.......................................................................................................................................................................................................................................................................................................................................................................
%    indy = find(Imaxx == max(Imaxx));
    I11 = I_(:, :, 1);
%    threshold = 0.8*max(max(I11)); % Imaxy*0.05;
    lt =  0.8*(max(max(I11)));%max(max(I_))*0.05;
    
    
    I0(I11> lt) = I11(I11> lt);
    % label;
    Yh = getlb(I11, lt, ah, bh, winL1h, winL2h);
    Yv = getlb(I11, lt, av, bv, winL1v, winL2v);
    Yd = getlb(I11, lt, ad, bd, winL1d, winL2d);

%    Indx = indx(indy);
%    center = [Indx', indy'];
%    [idx, idy] = find(I_ >= threshold);
%    I0(idx, idy, 1)= 1;

% for n1 = 1:length(h)
%      Sh = repmat(0.0000001, [1, n1]);
%      Sh(n1)= features(h{n1});
% end
% for n2 = 1:length(v)
%      Sv = repmat(0.0000001, [1, n2]);
%      Sv(n2)= features(v{n2});
% end
% for n3 = 1:length(d)
%     Sd = repmat(0.0000001, [1, n3]);
%     Sd(n3)= features(d{n3});
% end
%S = [Sh;Sv;Sd];

S= [h{nn}(:)';v{nn}(:)';d{nn}(:)'];
kse = S;
Nn = size(S, 2);
F = kse;
%F = F(F> 0);   
Y = [Yh(1:Nn);Yv(1:Nn);Yd(1:Nn)];

W= repmat(1/3, [3, Nn]);
N = 50;
% sum(mean(F>0==(Y==1),1))
% sum(mean(F>mean(F,1)==(Y==1),1))
% sum(mean(abs(repmat(1,[1,Nn])./(1+exp(-F*Nn)))>0.5==(Y==1),1))
F(abs(repmat(1,[1,Nn])./(1+exp(-F*Nn)))>0.5)= 1;
F(abs(repmat(1,[1,Nn])./(1+exp(-F/Nn)))<0.5)= 0;
%F(F>=0)= 1;
%F(F<0)= 0;
Beta = W;
while N >0 
    change = 0;
    for l = 1:size(F,1)
        w = W(l,:);
        err = w.*abs(F(l,:) - Y(l,:));
        beta = err./(1-err);
        wt = w.*beta;
        wt = wt/sum(wt);
        errt = wt.*abs(F(l,:) - Y(l,:));
        if sum(err>errt)~= 0      
            err(err >errt) = errt(err >errt);
            W(l,err >errt) = w(err >errt);
            Beta(l,err >errt) = beta(err >errt);
            change = l;
        end
    end
    N = N-1;
    if change == 0
        break;
    end
end
row = find(sum(Beta, 2) ~= 0,2);
Beta = Beta(row,:);
alpha = -log(Beta);
fnew = F(row,:); 
Ypred = sum(alpha.*fnew, 1);
Ylb = sum(Beta.*Y(row,:), 1);
YLB = zeros(1, size(Ylb,2));
YLB(Ylb >= 0.5*sum(Beta, 1)) = 1;
YLB(Ylb< 0.5*sum(Beta, 1)) = 0;
YPRED = zeros(1, size(Ypred,2));
YPRED (Ypred >= 0.5*sum(alpha, 1)) = 1;
YPRED (Ypred < 0.5*sum(alpha, 1)) = 0;

acc = sum(sum(YPRED == YLB))/length(Ylb)
train = YLB;
train(YLB == 1) = 255;
train = reshape(train, 515,650);
YPRED(YPRED == 1) = 255;
YPRED = reshape(YPRED, 515,650);


figure,
subplot(2,2,1)
imshow(imresize(I,0.5),[])
title('origin')
subplot(2,2,2)
imshow(imresize(I0,0.5),[])
title('labeled ROI of dentrites')
subplot(2,2,3)
imshow(train,[])
title('labeled dentrites after train')
subplot(2,2,4)
imshow(YPRED,[])
title('predicted dentrites')


[cx, cy] = find(I11> lt);
figure,
subplot(4,1,1)
scatter(cx,cy);
title('ROI points');  
subplot(4,1,2)
Z = linkage([cx, cy],'average','chebychev');
%Z = linkage([cx, cy],'ward','euclidean');
lastTen = Z(end-9:end,:)
plot(lastTen(:,3))
title('ward distance v.s clusters')
subplot(4,1,3)
c = cluster(Z, "maxclust",4, 'Depth',4);
crosstab(c)
cutoff = median([Z(end-3,3) Z(end-2,3) Z(end-1,3)]);
dendrogram(Z,'ColorThreshold',cutoff)

subplot(4,1,4)
c1 = cluster(Z, "maxclust",6, 'Depth',8);
crosstab(c1)
cutoff = median([Z(end-5,3) Z(end-4,3) Z(end-3,3) Z(end-2,3) Z(end-1,3)]);
dendrogram(Z,'ColorThreshold',cutoff)

% T = cluster(Z,'cutoff',3,'Depth',4);
% gscatter(X(:,1),X(:,2),T)

%{
lastTen =

                     10078                     10142          48.8455781823415
                     10123                     10143          54.4532826313052
                     10138                     10144          57.2559115762243
                      8753                     10145          99.8121730860675
                     10134                     10146          105.358383534137
                      9959                     10110            110.8836996337
                     10135                     10148          201.435291270902
                     10140                     10150          334.881857914797
                     10147                     10151           430.87552374308
                     10149                     10152          688.663277943107

ans =

           1
        1373
        2180
        1408
         115


ans =

           1
          24
          91
         837
        1343
        1373
        1408

}%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function S = features(Matrix)
[m, n] = size(Matrix);
S=sum(sum(repmat(1/(m*n),[m,n]).*Matrix));
end

function kse = Haarlikefeature(Image, winL)
[m, n] = size(Image);
I0 = Image;
%edge 1( 10:tlbr01, 11:bltr10, 12:tlbr10, 13:bltr01, 14:h01, 15:h10 ),line 2(20:v010, 21:v101, 22:tlbr101, 23: tlbr010, 24:h101, 25: trbl010, 26:h010), center-surround 3
I0_ = zeros(size(I0));
%idx = [idx - 1;max(idx);max(idx) + 1];
%idy = [idy - 1;max(idy);max(idy) + 1];
for i = 1: length(idx)
    for dx = min(max(idx(i) - 1,1),m):min(max(idx(i) + 1,1),m)
        for dy = min(max(idy(i) - 1,1),n):min(max(idy(i) + 1,1),n)
            for N = -1:1:1
                if I0(min(max(dx+N,1), m),min(max(dy+N,1),n))== 0
                    if I0(dx,min(max(dy+N,1), n))== 1
                        if sum(sum(I0(min(max(dx+2*N,1), m):N:dx,dy:N:min(max(dy+2*N,1),n))==0)) == 8
                            I0_(min(max(dx+2*N,1), m):N:dx,dy:N:min(max(dy+2*N,1),n)) = 3;
                        elseif I0(min(max(dx+2*N,1), m),min(max(dy+N,1), n)) == 1
                            I0_(min(max(dx+2*N,1), m),min(max(dy+N,1), n)) = 24;
                            I0_(min(max(dx+N,1), m),min(max(dy+N,1),n)) = 24;
                            I0_(dx,min(max(dy+N,1), n))= 24;
                        elseif I0(min(max(dx+2*N,1), m),min(max(dy+N,1),n)) == 0 || I0(min(max(dx-N,1), m),min(max(dy+N,1),n)) == 1
                            I0_(min(max(dx+2*N,1), m),min(max(dy+N,1), n)) = 14;
                            I0_(min(max(dx+N,1), m),min(max(dy+N,1),n)) = 14;
                        end
                    elseif I0(min(max(dx+2*N,1), m),min(max(dy+2*N,1),n)) == 0 || I0(min(max(dx-N,1), m),min(max(dy-N,1),n)) == 0
                        I0_(min(max(dx+2*N,1), m),min(max(dy+2*N,1),n)) = 12;
                        I0_(min(max(dx-N,1), m),min(max(dy-N,1),n)) = 12;
                        I0_(dx,dy) = 12;
                    elseif I0(min(max(dx+2*N,1), m),min(max(dy+2*N,1),n)) == 1
                        I0_(min(max(dx+2*N,1), m),min(max(dy+2*N,1),n)) = 22;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n))=22;
                        I0_(min(max(dx+N,1), m),min(max(dy+2*N,1)), n)=22;
                    elseif I0(dx, dy)== 1 && I0(min(max(dx+2*N,1), m), min(max(dy+2*N,1),n))== 1
                        I0_(min(max(dx+2*N,1),m),min(max(dy+2*N,1),n)) = 23;
                        I0_(dx, dy) = 23;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n))= 23;
                    elseif I0(min(max(dx+N,1), m),dy) == 1 && I0(min(max(dx+N,1), m),min(max(dy-N,1), n)) == 0
                        I0_(min(max(dx+N,1), m),dy) = 20;
                        I0_(min(max(dx+N,1), m),min(max(dy-N,1), n)) = 20;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n))= 20;
                    elseif I0(min(max(dx+N,1), m),min(max(dy+2*N,1),n)) == 1 && I0(min(max(dx+N,1), m),min(max(dy+3*N,1), n)) == 0
                        I0_(min(max(dx+N,1), m),min(max(dy+2*N,1),n)) = 20;
                        I0_(min(max(dx+N,1), m),min(max(dy+3*N,1), n)) = 20;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n))= 20;
                    elseif I0(min(max(dx+N,1), m),min(m+ax(dy+2*N,1),n)) == 1 && (I0(min(max(dx+N,1), m),min(max(dy+3*N,1), n)) == 1 || I0(min(max(dx+N,1), m),dy == 0))
                        I0_(min(max(dx+N,1), m),min(max(dy+2*N,1),n)) = 14;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n)) = 14;
                    elseif I0(min(max(dx+2*N,1), m),min(max(dy+N,1),n)) == 1 && I0(min(max(dx+3*N,1), m),min(max(dy+N,1), n)) == 0
                        I0_(min(max(dx+2*N,1), m),min(max(dy+N,1),n)) = 26;
                        I0_(min(max(dx+3*N,1), m),min(max(dy+N,1), n)) = 26;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n))= 26;
                    elseif I0(min(max(dx+2*N,1), m),min(max(dy+N,1),n)) == 1 && (I0(min(max(dx+3*N,1), m),min(max(dy+N,1), n)) == 1 || I0(dx,min(max(dy+N,1), n)) == 0)
                        I0_(min(max(dx+2*N,1), m),min(max(dy+N,1),n)) = 15;
                        I0_(min(max(dx+3*N,1), m),min(max(dy+N,1), n)) = 15;
                    elseif I0(min(max(dx+2*N,1), m),dy) == 1 && (I0(min(max(dx+3*N,1), m), min(max(dy+N,1), n)) == 1 || I0(dx, min(max(dy+2*N,1), n)) == 0)
                        I0_(min(max(dx+2*N,1), m),idy(i)) = 13;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n)) = 13;
                    elseif I0(dx,min(max(dy+2*N,1), n)) == 1 && (I0(min(max(dx+2*N,1), m),dy) == 0 || I0( min(max(dx-N,1), n), min(max(dy+3*N,1), n)) == 1)
                        I0_(min(max(dx+2*N,1), m),dy) = 11;
                        I0_(min(max(dx-N,1), n), min(max(dy+3*N,1), n)) = 11;
                    elseif I0(min(max(dx-N,1), m),min(max(dy+N,1),n)) == 0 && I0(dx,min(max(dy+N,1), n)) == 1
                        I0_(min(max(dx-N,1), m),min(max(dy+N,1),n)) = 26;
                        I0_(min(max(dx+3*N,1), m),min(max(dy+N,1), n)) = 26;
                        I0_(dx,min(max(dy+N,1), n))= 26;
                    elseif I0(min(max(dx+2*N,1), m),dy) == 1 && I0(dx,min(max(dy+2*N,1), n)) == 1
                        I0_(min(max(dx+2*N,1), m),dy) = 25;
                        I0_(dx,min(max(dy+2*N,1), n)) = 25;
                        I0_(min(max(dx+N,1), m),min(max(dy+N,1),n))=25;
                    elseif I0(min(max(dx+N,1), m),dy) == 1 && I0(min(max(dx+N,1), m), min(max(dy+2*N,1), n)) == 1
                        I0_(min(max(dx+N,1), m),dy) = 21;
                        I0_(min(max(dx+N,1),min(max(dy+N,1), n))) = 21;
                        I0_(min(max(dx+N,1), m),min(max(dy+2*N,1),n)) = 21;
                    elseif I0(min(max(dx+2*N,1), m),min(max(dy+2*N,1), n)) == 1 && (I0(dx,dy) == 0 || I0( min(max(dx+2*N,1),m), min(max(dy+2*N,1),n)) == 0) 
                        I0_(min(max(dx+2*N,1), m),min(max(dy+2*N,1), n)) = 10;
                        I0_(min(max(dx+N,1),m),min(max(dy+N,1), n)) = 10;
                    end
                end
            end
        end
    end
end
kse = zeros(uint32(m*n/winL^2), uint32(m/winL),uint32(n/winL), 2);
II = zeros(size(I0_));
II_ = zeros(size(I0_));
[xxi, yyi] = find(I0_ > 0);
        for i = max(min(xxi)-winL+1,1):winL:min(max(xxi)+winL-1,m)
            for j = max(min(yyi)-winL+1,1):winL:min(max(yyi)+winL-1,n)
                w0 = winL^2/(m*n);
                ui = 0;
                vi = 0;
                xi = i;%(i-(min(yyi)-winL+1))/winL+1;
                yi = j;%(j-(min(yyi)-winL+1))/winL+1;
                if I0(xi, yi) == 0
                    wi = -w0;  
                    f = 1; 
                    S = 0;
                end
                while  I0_(xi + ui, yi + vi) >= 0
%                     switch I0_(xi + ui, yi + vi)
%                         case 10
%                             f = 1;
%                         case 11
%                             f = 2;
%                         case 12
%                             f = 3;
%                         case 13
%                             f = 4;
%                         case 14
%                             f = 5;
%                         case 15
%                             f = 6;
%                         case 20
%                             f = 7;
%                         case 21
%                             f = 8;
%                         case 22
%                             f = 9;
%                         case 23
%                             f = 10;
%                         case 24
%                             f = 11;
%                         case 25
%                             f = 12;
%                         case 26
%                             f = 13;
%                         case 3
%                             f = 14;
%                     end
                    if vi < winL & yi + vi < n 
                        if I0_(xi + ui, yi + vi)== 0
                            while I0_(xi + ui, yi + vi + 1)== 0
                                vi = vi + 1;
                                if vi == winL & yi + vi < n
                                    if I0_(xi + ui, yi + vi -1)== 0
                                        break;
                                    else
                                        II_(xi + ui, yi + vi -1) = sum(II_(xi + ui, yi: yi+ vi -1));
                                        S = S + wi*(II_(xi + ui, yi + vi -1)+II_(xi + ui, yi));
                                        kse(f, xi, yi + vi -1, :) = [S, vi-1];
                                        yi = yi + vi; 
                                        f = f + 1;
                                        S = 0;
                                    end
                                end
                            end
                        end
                    end
                    II(xi + ui, yi + vi) = I0_(xi + ui, yi + vi);
                    if vi< winL & yi + vi < n
                        while I0_(xi + ui, yi + vi + 1) == I0_(xi+  ui, yi + vi)
                            while vi< winL & yi + vi < n 
                                if I0(xi+ ui, yi + vi +1) == I0(xi+ ui, yi+ vi)
                                    II_(xi + ui, yi + vi) = I0(xi+ ui, yi+ vi);
                                    vi = vi +1;
                                end
                            end
                            if I0(xi+ ui, yi + vi +1) ~= I0(xi+ ui, yi+ vi) & I0_(xi + ui, yi + vi +1) == I0_(xi+  ui, yi + vi)
                                II_(xi + ui, yi + vi) = I0(xi+ ui, yi+ vi);
                                II_(xi + ui, yi + vi) = sum(II(xi + ui, yi:yi + vi));
                                II(xi + ui, yi + vi) = I0(xi + ui, yi + vi);
                                S = S + wi*(II_(xi + ui, yi + vi)+II_(xi + ui, yi));
                                kse(f,  xi, yi + vi , :) = [S, vi];
                                yi = yi + vi;
                                winL = winL - vi;
                                vi = 0;
                                f = f+1;
                                wi = -wi;
                                S = 0;
                                while vi< winL & yi + vi + 1< n 
                                    if I0(xi+ ui, yi + vi +1) == I0(xi+ ui, yi+ vi)                               
                                         II_(xi + ui, yi + vi) = I0(xi+ ui, yi+ vi);
                                         vi = vi +1
                                    end
                                end
                                if I0(xi + ui, yi + vi +1) ~= I0(xi+  ui, yi + vi)
                                    II_(xi + ui, yi + vi) = I0(xi+ ui, yi+ vi);
                                    II_(xi + ui, yi + vi) = sum(II(xi + ui, yi:yi + vi));
                                    II(xi + ui, yi + vi) = I0(xi + ui, yi + vi);
                                    S = S + wi*(II_(xi + ui, yi + vi)+II_(xi + ui, yi));
                                    kse(f,  xi, yi + vi , :) = [S, vi];
                                    yi = yi + vi;
                                    winL = winL - vi;
                                    vi = 0;
                                    f = f+1;
                                    wi = -wi;
                                    S = 0;
                                end
                                while I0_(xi + ui, yi + vi)== 0
                                    vi = vi + 1;
                                    if vi == winL| yi + vi == n
                                       break;
                                    end
                                end
                            end  
                        end
                    end
                    if vi< winL & yi + vi + 1 < n
                        if I0_(xi + ui, yi + vi + 1) ~= I0_(xi+  ui, yi + vi) 
                            II(xi+ ui, yi + vi) = I0(xi+ ui, yi + vi);
                            II_(xi+ ui, yi + vi) = sum(II(xi+ ui, yi:yi + vi));
                            S = S + wi*(II_(xi, yi + vi)+II_(xi, yi));
                            kse(f, xi, yi + vi, :) = [S, vi];
                            yi = yi + vi;
                            f = f + 1;
                            vi = 0;
                        end
                        if I0_(xi + ui, yi + vi + 1)== 0
                            while vi< winL & yi + vi + 1 < n
                                if I0_(xi + ui, yi + vi + 1)== 0
                                    vi = vi + 1;
                                    if vi == winL| yi + vi == n
                                       break;
                                    end
                                end
                            end
                            if vi < winL & yi + vi < n
                                yi = yi + vi;
                                winL = winL - vi;
                                vi = 0;  
                            end
                        end
                    end
                    if vi == winL | yi + vi == n 
                        if I0_(xi + ui, yi + vi)== 0
                            yi = j; 
                            vi = 0;
                        else
                            II(xi+ ui, yi + vi) = I0(xi+ ui, yi + vi);
                            II_(xi+ ui, yi + vi) = sum(II(xi+ ui, yi:yi + vi));
                            S = S + wi*(II_(xi+ ui, yi + vi)+II_(xi+ ui, yi));
                            kse(f, xi+ ui, yi + vi, :) = [S, vi];
                            yi = yi + vi;
                            f = f + 1;
                            vi = 0;
                        end
                    end
                    while ui< winL & xi + ui < m 
                        if I0_(xi + ui +1, yi: yi + vi) == I0_(xi + ui, yi: yi + vi)
%                            xi = xi + ui;
                            II(xi + ui +1, yi: yi + vi) = I0_(xi + ui, yi: yi + vi); 
                            ui = ui + 1;
                            if I0(xi + ui, yi: yi + vi) ~= I0(xi + ui -1, yi: yi + vi)
                                II_(xi + ui -1, yi + vi) = sum(II(xi: xi + ui -1, yi + vi));
                                II_(xi, yi + vi) = sum(II(xi, yi: yi + vi));
                                S = S + wi*(II_(xi, yi) +  II_(xi + ui -1, yi + vi) - II_(xi, yi + vi) - II_(xi + ui -1, yi));
                                kse(f, xi+ ui -1, yi + vi, 1) = kse(f, xi+ ui -1, yi + vi, 1) +S;
                                f = f + 1;
                                wi = -wi;
                                II(xi + ui , yi + vi + 1) = I0(xi + ui, yi + vi + 1);
                                xi = xi+ ui;
                                ui = 0;
                                S = 0;
                            elseif  ui== winL | xi + ui== m
                                II_(xi + ui, yi) = sum(II(xi: xi + ui, yi + vi));
                                kse(f, xi+ ui -1, yi + vi, 1) = kse(f, xi+ ui -1, yi + vi, 1) - wi*(II_(xi, yi) +  II_(xi + ui -1, yi + vi) - II_(xi, yi + vi) - II_(xi + ui -1, yi));
                                if xi + ui < m 
                                    xi = xi+ ui;
                                    II(xi + ui, yi + vi) = I0(xi: xi + ui, yi + vi);
                                    f = f+1;
                                else
                                    f = 1;
                                end                                    
                            elseif  I0_(xi + ui, yi: yi + vi) ~= I0_(xi + ui -1, yi: yi + vi)
                                II_(xi + ui -1, yi + vi) = sum(II(xi: xi + ui -1, yi + vi));
                                II_(xi, yi + vi) = sum(II(xi, yi: yi + vi));
                                II_ (xi + ui, yi + vi) = I0(xi + ui, yi + vi);
                                S = S - wi*(II_(xi, yi) +  II_(xi + ui -1, yi + vi) - II_(xi, yi + vi) - II_(xi + ui -1, yi));
                                kse(f, xi+ ui -1, yi + vi, 1) = kse(f, xi+ ui -1, yi + vi, 1) - wi*(II_(xi, yi) +  II_(xi + ui -1, yi + vi) - II_(xi, yi + vi) - II_(xi + ui -1, yi));
                                S = 0;
                                f = f + 1;
                                xi = xi + ui;
                                yi = yi + vi;
                                ui = 0;    
                                vi = 0;    
                            end
                        end
                    end
                end                 
            end
         end
end




%I1(:,:,3)= I1(:,:,2);
figure,
subplot(121)
imshow(I, [])
subplot(122)
imshow(I1,[])