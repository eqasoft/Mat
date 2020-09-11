function Kd = klCompose(f,nElem,Material,Lx,klorder,KL_xi,lambda)


if length(f)>1

    Kd = sfem_fun_beam3(f,nElem,Material{1},Lx);

    for iKL = 1:klorder
     Kd_KL = sfem_fun_beam3(f,nElem,Material{iKL+1},Lx);
     for iKd = 1:length(Kd)
      Kd{iKd} = Kd{iKd} + sqrt(lambda(iKL))*Kd_KL{iKd}*KL_xi(iKL);
     end
    end

else
    
    Kd = sfem_fun_beam3(f,nElem,Material{1},Lx);
    for iKL = 1:klorder
      Kd_KL = sfem_fun_beam3(f,nElem,Material{iKL+1},Lx);
      Kd = Kd + sqrt(lambda(iKL))*Kd_KL*KL_xi(iKL);
    end

end