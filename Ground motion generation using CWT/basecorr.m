function [vel,despl,cxg,cvel,cdespl] = basecorr(t,xg,CT)

% perform baseline correction

% luis.montejo@upr.edu

imax = 80 ;    % maximum number of iterations
tol  = 0.001; % tolerance (percent of the max)

dt    = t(2)-t(1); np = length(xg); cxg = zeros(size(xg));
vel    = dt*cumtrapz(xg);
despl    = dt*cumtrapz(vel); 

L  = round(CT/(dt))-1;
M  = np-L;

for q = 1:imax

  dU = 0;
  ap = 0;
  an = 0;

  dV = 0;
  vp = 0;
  vn = 0;
  for i = 1:(np-1)
      dU = dU + (t(np)-t(i+1)) * xg(i+1) * dt;
  end

  for i = 1:L+1
      aux = ((L-(i-1))/L)*(t(np)-t(i))*xg(i)*dt;
      if aux >= 0
          ap = ap + aux;
      else
          an = an + aux;
      end
  end

  alfap = -dU/(2*ap);
  alfan = -dU/(2*an);

  for i =2:np
      if i<=L+1
          if xg(i)>0
              cxg(i) = (1 + alfap*(L-(i-1))/L) * xg(i);
          else
              cxg(i) = (1 + alfan*(L-(i-1))/L) * xg(i);
          end
      end
      if i>L+1
          cxg(i) = xg(i);
      end
  end

  xg = cxg;

  %==========================================================

  for i = 1:(np-1)
      dV = dV + xg(i+1) * dt;
  end


  for i = M:np
      auxv = ((i - M)/(np-M))*xg(i)*dt;
      if auxv >= 0
          vp = vp + auxv;
      else
          vn = vn + auxv;
      end
  end

  valfap = -dV/(2*vp);
  valfan = -dV/(2*vn);

  for i =2:np
      if i>=M
          if xg(i)>0
              cxg(i) = (1 + valfap*((i - M)/(np-M))) * xg(i);
          else
              cxg(i) = (1 + valfan*((i - M)/(np-M))) * xg(i);
          end
      end
      if i<M
          cxg(i) = xg(i);
      end
  end
  
  xg     = cxg;
  cvel   = dt*cumtrapz(xg);
  cdespl = dt*cumtrapz(cvel); 
  errv(q) = abs(cvel(length(cvel))/ max(abs(cvel)));
  errd(q) = abs(cdespl(length(cdespl)) / max(abs(cdespl)));

  if errv(q) <= tol & errd(q) <= tol, break, 
  else 
    xg = cxg;
  end
  
end