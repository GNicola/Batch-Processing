function [Tb,Tf] = DIVTIME(Time,n); 
 k = 0;
 m = 0;

 for ii=1 : n
   c = Time(ii);
   if (c < 0)
     k = k + 1;
     T(k,1) = c;
   elseif (c > 0)
     m = m + 1;
     Tf(m,1) = c;
   else 
     k = k + 1;
     m = m + 1;
     Tf(m,1) = c;
     T(k,1) = c;
   end
 end

 for ii=1 : k	
   Tb(ii,1) = T(k+1-ii,1);
 end

 k
 m

end
