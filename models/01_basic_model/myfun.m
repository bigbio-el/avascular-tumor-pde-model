
function [F]= myfun(x,DX,Dc,alpha,gamma,k1,k2,k3,c0,N)
   

         
   for i = 1 : (N)       
         
     % Forward Finite Difference Approximation
    
      if i == 1
          
        P1 = x(i,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
        P1_l = x(i + 1,1)/(x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1));
        PQR = x(i,1) + x(i + N,1) + x(i + 3*N,1);
        PQR_l = x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1);
        Q1 = x(i + N,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
        Q1_l = x(i + N + 1,1)/(x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1));
        S1 = x(i + 3*N,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
        S1_l = x(i + 3*N + 1,1)/(x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1));
        PQNS_s = x(i,1) + x(i + N,1) + x(i + 2*N,1) + x(i + 3*N,1);
        
        
        dp(i,1) = ((P1_l - P1)/DX) * ((PQR_l - PQR)/DX) + P1 * ((PQR_l - 2 * PQR + PQR_l)/DX^2) + G(x(i + 4*N,1)) * x(i,1) * (1 - PQNS_s) - Z(x(i + 4*N,1)) * x(i,1);
        
        dq(i,1) = ((Q1_l - Q1)/DX) * ((PQR_l - PQR)/DX) + Q1 * ((PQR_l - 2 * PQR + PQR_l)/DX^2) + Z(x(i + 4*N,1)) * x(i,1) - H(x(i + 4*N,1)) * x(i + N,1);  
        
        dn(i,1) = H(x(i + 4*N,1)) * x(i + N,1);
        
        ds(i,1) = ((S1_l - S1)/DX) * ((PQR_l - PQR)/DX) + S1 * ((PQR_l - 2 * PQR + PQR_l)/DX^2) + G(x(i + 4*N,1)) * x(i + 3*N,1) * (gamma - PQNS_s);
        
        dc(i,1) = Dc *(x(i + 4*N + 1,1)-2 * x(i + 4*N,1) + x(i + 4*N + 1,1))/(DX^2) + k1 * c0 * (1 - alpha * PQNS_s) - k1 * x(i + 4*N,1) - k2 * x(i,1) * x(i +4*N,1) - k3 * x(i + 3*N,1) * x(i + 4*N,1);
                
      end;
      
     % Central Finite Difference Approximation
      
      if ((i > 1) && (i < N))
          
         P1 = x(i,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
         P1_l = x(i + 1,1)/(x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1));
         P1_r = x(i - 1,1)/(x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1));       
         PQR = x(i,1) + x(i + N,1) + x(i + 3*N,1);
         PQR_l = x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1);
         PQR_r = x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1);
         Q1 = x(i + N,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
         Q1_l = x(i + N + 1,1)/(x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1));
         Q1_r = x(i + N - 1,1)/(x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1));   
         S1 = x(i + 3*N,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
         S1_l = x(i + 3*N + 1,1)/(x(i + 1,1) + x(i + N + 1,1) + x(i + 3*N + 1,1));
         S1_r = x(i + 3*N - 1,1)/(x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1));   
         PQNS_s = x(i,1) + x(i + N,1) + x(i + 2*N,1) + x(i + 3*N,1);
        
        
         dp(i,1) = ((P1_l - P1_r)/(DX*2)) * ((PQR_l - PQR_r)/(DX*2)) + P1 * ((PQR_l - 2 * PQR + PQR_r)/DX^2) + G(x(i + 4*N,1)) * x(i,1) * (1 - PQNS_s) - Z(x(i + 4*N,1)) * x(i,1);
        
         dq(i,1) = ((Q1_l - Q1_r)/(DX*2)) * ((PQR_l - PQR_r)/(DX*2)) + Q1 * ((PQR_l - 2 * PQR + PQR_r)/DX^2) + Z(x(i + 4*N,1)) * x(i,1) - H(x(i + 4*N,1)) * x(i + N,1);  
        
         dn(i,1) = H(x(i + 4*N,1)) * x(i + N,1);
        
         ds(i,1) = ((S1_l - S1_r)/(DX*2)) * ((PQR_l - PQR_r)/(DX*2)) + S1 * ((PQR_l - 2 * PQR + PQR_r)/DX^2) + G(x(i + 4*N,1)) * x(i + 3*N,1) * (gamma - PQNS_s);
        
         dc(i,1) = Dc *(x(i + 4*N + 1,1)-2 * x(i + 4*N,1) + x(i + 4*N - 1,1))/(DX^2) + k1 * c0 * (1 - alpha * PQNS_s) - k1 * x(i + 4*N,1) - k2 * x(i,1) * x(i +4*N,1) - k3 * x(i + 3*N,1) * x(i + 4*N,1);
         
      end;     
      
     % Backward Finite Difference Approximation
    
      if i == N
          
         P1 = x(i,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
         P1_r = x(i - 1,1)/(x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1));
         PQR = x(i,1) + x(i + N,1) + x(i + 3*N,1);
         PQR_r = x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1);
         Q1 = x(i + N,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
         Q1_r = x(i + N - 1,1)/(x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1));
         S1 = x(i + 3*N,1)/(x(i,1) + x(i + N,1) + x(i + 3*N,1));
         S1_r = x(i + 3*N - 1,1)/(x(i - 1,1) + x(i + N - 1,1) + x(i + 3*N - 1,1));
         PQNS_s = x(i,1) + x(i + N,1) + x(i + 2*N,1) + x(i + 3*N,1);
        
        
         dp(i,1) = ((P1 - P1_r)/DX) * ((PQR - PQR_r)/DX) + P1 * ((PQR_r - 2 * PQR + PQR_r)/DX^2) + G(x(i + 4*N,1)) * x(i,1) * (1 - PQNS_s) - Z(x(i + 4*N,1)) * x(i,1);
        
         dq(i,1) = ((Q1 - Q1_r)/DX) * ((PQR - PQR_r)/DX) + Q1 * ((PQR_r - 2 * PQR + PQR_r)/DX^2) + Z(x(i + 4*N,1)) * x(i,1) - H(x(i + 4*N,1)) * x(i + N,1);  
        
         dn(i,1) = H(x(i + 4*N,1)) * x(i + N,1);
        
         ds(i,1) = ((S1 - S1_r)/DX) * ((PQR - PQR_r)/DX) + S1 * ((PQR_r - 2 * PQR + PQR_r)/DX^2) + G(x(i + 4*N,1)) * x(i + 3*N,1) * (gamma - PQNS_s);
        
         dc(i,1) = Dc *(x(i + 4*N - 1,1)-2 * x(i + 4*N,1) + x(i + 4*N - 1,1))/(DX^2) + k1 * c0 * (1 - alpha * PQNS_s) - k1 * x(i + 4*N,1) - k2 * x(i,1) * x(i + 4*N,1) - k3 * x(i + 3*N,1) * x(i + 4*N,1);
                  
      end;
    
      % Create the Vectors
                 
   end;      
   
   F = [dp;dq;dn;ds;dc];