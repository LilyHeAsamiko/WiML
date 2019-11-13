
  %
   Ipred= imresize(Ipred, 0.1);
  I11 = imresize(I11, 0.1);
  Ipred= sparse(double(Ipred(:)));
  I11 = sparse(double(I11(:)));
 [h, p] = ttest(Ipred , I11)
 
  %discrete
  
 pcd = I11/max(I11);
 yc = ones(1,length(pcd))./(1+exp(-pcd));
 yc1 = ones(1,length(pcd))./(1+exp(-pcd-log(1)));
 
 pred = Ipred/max(Ipred);
 ypred = ones(1,length(pred))./(1+exp(-pred));
 ypred1 = ones(1,length(pred))./(1+exp(-pred-log(1)));

 [h, p] = ttest(yc, ypred)% h = 1, significantly different

%F0
 fa_ac = (yc+pcd).*(ypred+pred);
 %F01
 fa_ac1 = (yc1+pcd).*(yc1-pcd*log(1));
 %F10
 f1a1_ac = (ypred-pred*log(1)).*(ypred+pred);
 %F1
 f1a1_ac1 = (yc-pcd*log(1)).*(ypred-pred*log(1));


 
 % fk = (gam(1) - pad(1:4).*log(1.1)).*(gam(2) - pacd(1:4).*log(1.1));
 F0 = [mean(fa_ac(:)), mean(fa_ac1(:))];
 F1 = [mean(f1a1_ac(1)), mean(f1a1_ac1(1))];
 

 nH = sum([ttest(F0(1),F1(1))==0 ,ttest(F0(2),F1(1))==0 ,ttest(F0(1),F1(2))==0 ,ttest(F0(2),F1(2))==0]);
 phi =  3/4+0.5*nH/4
 psi = 0.5+0.5*(3/4-1/4)
 
 phi_ = mean(f1a1_ac1)-mean(fa_ac)+0.5
 psi_ = max(F0(:)-F1(:))
 
 
 
 

 
  %
 
 a_80 = [-4, -3, -2, -1, 2, 2, 3, 3, 6, 7 ,7, 10, 14, 14, 17, 18, 24, 26, 29, 30, 31, 55, 72];
 a_80c = [30, 33, 36, 47, 47, 48, 48, 49, 58, 60, 64, 64, 67, 72, 75, 76, 77, 79, 84, 86, 86, 88, 88, 89, 96];
 
 theta = 1;
 K =2;
 %discrete
 pad = a_80/max(a_80);
 ya = ones(1,length(pad))./(1+exp(-pad));
 ya1 = ones(1,length(pad))./(1+exp(-pad-log(theta)));
  
 pacd = a_80c/max(a_80c);
 yac = ones(1,length(pacd))./(1+exp(-pacd));
 yac1 = ones(1,length(pacd))./(1+exp(-pacd-log(theta)));   

[h, p] = ttest(ya, yac(1:23))% h = 1, significantly different
 
 gam=[-3.101, 0.9425];
 %F0
 fa_ac = (ya+pad).*(yac(1:23)+pacd(1:23));
 %F01
 fa_ac1 = (ya+pad).*(yac1(1:23)-pacd(1:23)*log(1));
 %F10
 f1a1_ac = (ya1-pad*log(1)).*(yac(1:23)+pacd(1:23));
 %F1
 f1a1_ac1 = (ya1-pad*log(1)).*(yac1(1:23)-pacd(1:23)*log(1));

 
 % fk = (gam(1) - pad(1:4).*log(1.1)).*(gam(2) - pacd(1:4).*log(1.1));
 F0 = [mean(fa_ac(:)), mean(fa_ac1(:))];
 F1 = [mean(f1a1_ac(1)), mean(f1a1_ac1(1))];
 

 nH = sum([ttest(F0(1),F1(1))==0 ,ttest(F0(2),F1(1))==0 ,ttest(F0(1),F1(2))==0 ,ttest(F0(2),F1(2))==0]);
 phi =  3/4+0.5*nH/4;
 psi = 0.5+0.5*(3/4-1/4);
 
 phi_ = mean(fa_ac)-mean(f1a1_ac1)+0.5
 psi_ = max(F0(:)-F1(:))