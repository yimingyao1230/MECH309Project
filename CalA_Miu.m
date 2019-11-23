function [A, miu] = CalA_Miu(phi, i, j, Uinf, Minf, gamma, dx, dy)
    %A judgement
 
    u_ = (phi(j,i-1) - phi(j,i+1))/(2*dx) ;
    v_ = (phi(j-1,i) - phi(j-1,i))/(2*dy) ;
    u = (u_ + Uinf);
    U = sqrt (u^2 + v_^2);
    m = U / c;
    if m > 1 % supersonic locally
        A (j,i) = (1 - Minf)^ 2 - (gamma + 1) * Minf^2 / Uinf * (phi (j,i) - phi (j,i-2) ) / ( 2* dx);
    else % subsonic flow
        A (j,i) = (1 - Minf)^ 2 - (gamma + 1) * Minf^2 / Uinf * (phi (j,i+1) - phi (j,i-1) ) / ( 2* dx);
    end

    if A(j,i) > 0
        miu(j,i) = 0;
    else
        miu(j,i)= 1;
    end
end 