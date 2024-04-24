
% mac=0.2;
% M=0.023;
% Ru=8.314;
% D=0.0213;


x_cont=max(find(T_surf>Ttr));




    %%%%%%%%%%%%%%%%%%  求声速传热极限  %%%%%%
    p_vapor=exp(11.9463-12633.73./T_vapor-0.4672*log(T_vapor))*1e6;
    rho_vapor=exp(-57.566+0.18157*T_vapor-2.2885e-4*power(T_vapor,2)+1.5614e-7*power(T_vapor,3)-5.5058e-11*power(T_vapor,4)+7.8615e-15*power(T_vapor,5))*0.001;
    hfg_vapor=(393.37*(1-T_vapor/2503.7)+4398.6*(1-T_vapor/2503.7)^0.29302)*1000;

    Q_sonic=rho_vapor*3.1415926*rMin*rMin*hfg_vapor*sqrt(1.333*Ru/2/2.333/M*T_vapor);     %% 声速传热极限 W
    %%%%%%%%%%%%%%%%%%%%


    if Q_evap_sonic==0
        Rg=Q_evap/delta_T_vapor;
    else
        Rg=Q_evap_sonic/delta_T_vapor;
    end

    Q_lv(1,1:x_cont)=-q_lv(1:x_cont)*2*3.1415926*rMin*dx;
    T_surf_new(1,1:x_cont)=T_surf(1:x_cont)-Q_lv(1:x_cont)/Rg;
    p_surf_new=exp(11.9463-12633.73./T_surf_new-0.4672*log(T_surf_new))*1e6;   
    
    Q_evap_sonic=0;

%%%%%%%%%%%%%%%%  根据已知界面温度 和蒸汽能量守恒关系，求解恒定蒸汽温度
iterMax = 10000;
for iter = 1: iterMax
     if mflag(1,n_points_x)==0
        delta_T_vapor=T_vapor-Ttr;
     else
        delta_T_vapor=T_surf_new(1)-T_surf_new(n_points_x);
     end
     q_evap=0;
     q_cond=0;
     p_vapor=exp(11.9463-12633.73./T_vapor-0.4672*log(T_vapor))*1e6;
     hfg_vapor=(393.37*(1-T_vapor/2503.7)+4398.6*(1-T_vapor/2503.7)^0.29302)*1000;
     q_lv(1:x_cont) =alpha_evapcond*(p_surf_new(1:x_cont)./sqrt(T_surf_new(1:x_cont))-p_vapor./sqrt(T_vapor))*hfg_vapor;

%  q_lv(1:x_cont) =alpha_evapcond*(p_surf(1:x_cont)./sqrt(T_surf(1:x_cont))-p_vapor./sqrt(T_vapor))*hfg_vapor;
        for i=1:x_cont
            if q_lv(i)>0
            q_evap=q_evap+q_lv(i);  %%%  实际传热量  W
            else
            q_cond=q_cond+q_lv(i);    
            end
        end
     Q_evap=q_evap*2*3.1415926*rMin*dx;
        if q_evap>abs(q_cond)
            T_vapor=T_vapor+0.0001;
        else
            T_vapor=T_vapor-0.0001;
        end
    
%         if abs(q_cond+q_evap)<300
%             break;
%         end
   
end



%%%%%%%%%%%%%%%%  判断计算传热量是否超过声速传热极限
    if (Q_evap>Q_sonic)&&(Q_evap<Q_in)

    %%%%  迭代求解声速传热极限时的恒定蒸汽温度
    for iter_sonic=1:10000
            
            p_vapor=exp(11.9463-12633.73./T_vapor-0.4672*log(T_vapor))*1e6;
            hfg_vapor=(393.37*(1-T_vapor/2503.7)+4398.6*(1-T_vapor/2503.7)^0.29302)*1000;
            Q_lv=0;
            q_lv(1:x_cont) =alpha_evapcond*(p_surf_new(1:x_cont)./sqrt(T_surf_new(1:x_cont))-p_vapor./sqrt(T_vapor))*hfg_vapor;
            for i=1:x_cont
                if q_lv(i)>0
                    Q_lv=Q_lv+q_lv(i)*2*3.1415926*dx*rMin;
                end
            end
    
            if Q_lv<Q_sonic
                T_vapor=T_vapor-0.01;
            else
                T_vapor=T_vapor+0.01;
            end

            if abs(Q_lv-Q_sonic)<3
                break;
            end

    end

    
    Q_evap_sonic=Q_lv;
    Q_cond_sonic=0;
   
    for i=1:x_cont
        if q_lv(i)<0
         Q_cond_sonic=Q_cond_sonic+q_lv(i)*2*3.1415926*rMin*dx;      %%%  实际传热量  
        end
    end

    Q_evap_sonic/Q_cond_sonic
    for i=1:x_cont  
        if q_lv(i)<0
            if Q_cond_sonic~=0
                 q_lv(i)=-q_lv(i)*Q_evap_sonic/Q_cond_sonic;  %%%  实际传热量  W
            end
        end
    end
  
    end
