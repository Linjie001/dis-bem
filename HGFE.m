function [Hi,Gi,VFi]=HGFE(GaussXW,Node,LEC,Mat,AG,ElNo,np,Node_dis)
% ENCo,xp,GaussXW,GaussPt
 global Element d1 d2
    	En1=Element(ElNo,13);  En2=Element(ElNo,14);   
       xp=Node_dis(np,2:3)';    
	   Den=Mat(1);    NU=Mat(3);    SG=Mat(4);
       Hi=zeros(2,4);	Gi=zeros(2,4);	VFi=zeros(2,1);
       c1=1/(d1+d2);
       DLE=Node(2,2:3)-Node(1,2:3);
       LEE=sum(DLE.*DLE,2).^0.5;
       DLE=DLE./[LEE LEE];
       NEE=DLE*[0 -1;1 0];
       Cor=[0.013320243;0.079750427;0.19787102;0.35415398;0.52945857;0.70181452;0.84937932;0.95332645];
       Wi=[0.16441660;0.23752560;0.22684198;0.17575408;0.11292402;0.057872212;0.020979074;0.0036864071];
       %Cor=[0.0090425944 0.053971054 0.13531134 0.24705169 0.38021171 0.52379159 0.66577472 0.79419019 0.89816102 0.96884798];
       %Wi=[0.12095474 0.18636310 0.19566066 0.17357723 0.13569597 0.093647084 0.055787938 0.027159893 0.0095151992 0.0016381586];
  for iG=1:size(GaussXW,2)
        x=GaussXW(1,iG);
        xl=([1-x x+1]/2*Node(1:2,2:3)).'; 
        r=xl-xp;
        rl=sqrt(r.'*r);            
        drd1=r(1)/rl;
        drd2=r(2)/rl;
        drdn=[drd1 drd2]*NEE';
        NShape=[c1*(d2-x) 0 c1*(d1+x) 0;0 c1*(d2-x) 0 c1*(d1+x)];
        %求H矩阵
        if   En1~=np && En2~=np
             p11=drdn*((1-2*NU)+2*drd1^2);
             p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
             p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
             p22=drdn*((1-2*NU)+2*drd2^2); 
             P=-1/4/pi/(1-NU)/rl*LEE/2*[p11 p12;p21 p22];
              FPN=P*NShape;  
             Hi=Hi+FPN*GaussXW(2,iG);
              %求G矩阵
            u11=(4*NU-3)*log(rl)+drd1^2-(7-8*NU)/2;
            u12=drd1*drd2;
            u21=u12;
            u22=(4*NU-3)*log(rl)+drd2^2-(7-8*NU)/2;
            U=1/8/pi/SG/(1-NU)*LEE/2*[u11 u12;u21 u22];
            Gi=Gi+U*NShape*GaussXW(2,iG); 
        elseif En1==np || En2==np
               if En1==np %单元首节点为源点
                     for nreg=1:2
                         if  nreg==1
                                xsi=-(d1+1)/2+ x*(1-d1)/2;
                                 dxdxb= (1-d1)/2;
                         else
                                xsi= (1-d1)/2+x*(1+d1)/2;  
                               dxdxb=(1+d1)/2;
                         end
                           xl=([1-xsi xsi+1]/2*Node(1:2,2:3)).'; 
                           r=xl-xp;
                           rl=sqrt(r.'*r);            
                           drd1=r(1)/rl;
                           drd2=r(2)/rl;
                           drdn=[drd1 drd2]*NEE';
                           NShape=[c1*(d2-xsi) 0 c1*(d1+xsi) 0;0 c1*(d2-xsi) 0 c1*(d1+xsi)];
                           %求H矩阵，形函数为零
                           p11=drdn*((1-2*NU)+2*drd1^2);
                           p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
                           p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
                           p22=drdn*((1-2*NU)+2*drd2^2); 
                           P=-1/4/pi/(1-NU)/rl*LEE/2*[p11 p12;p21 p22];
                           FPN=P*NShape;              
                             Hi(:,3:4)=Hi(:,3:4)+FPN(:,3:4)*GaussXW(2,iG)*dxdxb;
                             %求G矩阵，非奇异部分
                           u11=drd1^2-(7-8*NU)/2;
                           u12=drd1*drd2;
                           u21=u12;
                           u22=drd2^2-(7-8*NU)/2;
                           U=1/8/pi/SG/(1-NU)*LEE/2*[u11 u12;u21 u22];
                         %  NShape=[c1*(d2-xsi) 0 c1*(d1+xsi) 0;0 c1*(d2-xsi) 0 c1*(d1+xsi)];
                           FUN=U*NShape;
                          Gi=Gi+FUN*GaussXW(2,iG)*dxdxb; 
                             %求G矩阵，奇异部分，形函数为零
                           u1=(4*NU-3)*log(rl);
                           U1=1/8/pi/SG/(1-NU)*LEE/2*[u1 0;0 u1];
                           FUN1=U1*NShape;
                           Gi(:,3:4)=Gi(:,3:4)+FUN1(:,3:4)*GaussXW(2,iG)*dxdxb;
                     end
             elseif En2==np
                    for nreg=1:2
                        if nreg==1
                            xsi= (d2-1)/2+ x*(1+d2)/2;  
                             dxdxb=(1+d2)/2;
                        else
                             xsi= (d2+1)/2+x*(1-d2)/2;  
                               dxdxb=(1-d2)/2;  
                        end
                              xl=([1-xsi xsi+1]/2*Node(1:2,2:3)).'; 
                              r=xl-xp;
                               rl=sqrt(r.'*r);            
                              drd1=r(1)/rl;
                             drd2=r(2)/rl;
                             drdn=[drd1 drd2]*NEE';
                             NShape=[c1*(d2-xsi) 0 c1*(d1+xsi) 0;0 c1*(d2-xsi) 0 c1*(d1+xsi)];
                            %求H矩阵，形函数为零
                              p11=drdn*((1-2*NU)+2*drd1^2);
                              p12=drdn*(2*drd1*drd2)-(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
                              p21=drdn*(2*drd1*drd2)+(1-2*NU)*(drd1*NEE(2)-drd2*NEE(1));
                              p22=drdn*((1-2*NU)+2*drd2^2); 
                              P=-1/4/pi/(1-NU)/rl*LEE/2*[p11 p12;p21 p22];
                             FPN=P*NShape;
                            Hi(:,1:2)=Hi(:,1:2)+FPN(:,1:2)*GaussXW(2,iG)*dxdxb;
                           %求G矩阵，非奇异部分
                            u11=drd1^2-(7-8*NU)/2;
                           u12=drd1*drd2;
                           u21=u12;
                           u22=drd2^2-(7-8*NU)/2;
                           U=1/8/pi/SG/(1-NU)*LEE/2*[u11 u12;u21 u22];
                           NShape=[c1*(d2-xsi) 0 c1*(d1+xsi) 0;0 c1*(d2-xsi) 0 c1*(d1+xsi)];
                           FUN=U*NShape;
                          Gi=Gi+FUN*GaussXW(2,iG)*dxdxb; 
                          %求G矩阵，奇异部分，形函数为零
                             u2=(4*NU-3)*log(rl);
                             U2=1/8/pi/SG/(1-NU)*LEE/2*[u2 0;0 u2];
                              FUN1=U2*NShape;
                              Gi(:,1:2)=Gi(:,1:2)+FUN1(:,1:2)*GaussXW(2,iG)*dxdxb;
                    end
               
               end
            

            %U奇异部分，形函数非零，分离出的非奇异部分
                 if En1==np
                    for  nreg=1:2
                         if  nreg==1
                              xsi=-d1-(1-d1)*(1-x)/2;
                             dxdxb= -1+d1;
                             dxbdxp=-0.5;
                         else
                                                   %奇异部分的非奇异部分
                               xsi= -d1+(1+d1)*(1+x)/2;
                               dxdxb= 1+d1;
                               dxbdxp=0.5;
                         end
                          ra=LEE/2*abs(dxdxb);
                           u1=1/8/pi/SG/(1-NU)*LEE/2*(4*NU-3)*c1*(d2-xsi)*log(ra)*dxbdxp*dxdxb*GaussXW(2,iG);
                           Gi(1,1)=Gi(1,1)+u1;Gi(2,2)=Gi(2,2)+u1;
                    end
                elseif En2==np
                     for  nreg=1:2
                         if  nreg==1
                            xsi=d2-(1-x)/2*(1+d2);
                            dxdxb= -(1+d2);
                             dxbdxp=-0.5;
                         else
                                 xsi= d2+(1-d2)*(1+x)/2;
                                dxdxb= (1-d2);
                               dxbdxp=0.5;
                         end
                              ra=LEE/2*abs(dxdxb);
                               u1=1/8/pi/SG/(1-NU)*LEE/2*(4*NU-3)*c1*(d1+xsi)*log(ra)*dxbdxp*dxdxb*GaussXW(2,iG);
                               Gi(1,3)=Gi(1,3)+u1;Gi(2,4)=Gi(2,4)+u1;
                     end
                 end
                 %求G矩阵，奇异部分的特殊积分
          if  iG==1 && (En1==np || En2==np)
               for ig=1:8
                   x=Cor(ig);
                  if En1==np
                       for nreg=1:2
                          if  nreg==1
                                 xsi=-d1-(1-d1)*x;
                                 dxdxb= (1-d1);
                          else
                             xsi= -d1+(1+d1)*x;
                             dxdxb= 1+d1;
                          end
                                 u1=1/8/pi/SG/(1-NU)*LEE/2*(3-4*NU)*c1*(d2-xsi)*dxdxb*Wi(ig);
                               Gi(1,1)=Gi(1,1)+u1;Gi(2,2)=Gi(2,2)+u1;
                               a(nreg)=u1;
                       end 
                  else
                      for nreg=1:2
                          if   nreg==1
                              xsi=d2-x*(1+d2);
                              dxdxb= (1+d2);
                         else
                             xsi= d2+(1-d2)*x;
                             dxdxb= (1-d2);
                         end
                            u1=1/8/pi/SG/(1-NU)*LEE/2*(3-4*NU)*c1*(d1+xsi)*dxdxb*Wi(ig);
                            
                           Gi(1,3)=Gi(1,3)+u1;Gi(2,4)=Gi(2,4)+u1;
                      end
                  end  
               end
          end
        end
  end
                 
	 %%计算体积力Bi
	    VForce=rl*(2*log(rl)+1)*(Den*AG*([drd1 drd2]*NEE')-...
                             1/2/(1-NU)*Den*AG'*[drd1 drd2]'*NEE');
     	VFi=VFi-VForce*GaussXW(2,iG)*LEE/2/8/pi/SG ;                            
  

