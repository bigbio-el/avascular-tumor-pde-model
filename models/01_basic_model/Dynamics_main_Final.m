

% Sherratt, J., Chaplain, M., 2001. A new mathematical model for avascular
% tumour growth. J. Math. Biol.43(2), 291-312



% dp - Cell density of proliferating cells
            
% dq - Cell density of quiescent cells

% dn - Cell density of necrotic cells

% ds - Cell density of surrounding tissue

% dc - Nutrient/growth factor concentration
            
            
            
%___________________________________________________________________________________%
% BODY OF THE PROGRAM
%___________________________________________________________________________________%

                                 clc
                                 clear
%___________________________________________________________________________________%
% Initial Parameters:
%___________________________________________________________________________________%

t1 = 0;                          % First time point
t2 = 52.5;                       % Last time point
h = 0.0005;                      % Time step
N = 150;                         % Number of nodes

Dc = 10;                         % nutrient diffusion coefficient
alpha = 0.8;                     % 
gamma = 0.4;                     % factor, reflecting that surrounding epithelial cells, motile themselves, will divide at a rate that depends on nutrient level, 
                                 % but that is intrinsically lower than that of proliferating tumour cells  
                                 
k1 = 8;                          % correlation coefficient 1
k2 = 1;                          % correlation coefficient 2
k3 = 1;                          % correlation coefficient 3
c0 = 1;                          % nutrient concentration in the absence of a tumour cell population

X0 = 0;                          % Left boundary
L0 = 265;                        % Right boundary
DX = (L0 - X0)/(N - 1);          % Step difference
XX = X0 : DX : L0;               % Space (x) range  
error = 0.1;                     % Error
step = 0.5;                      % In loop-working step

%___________________________________________________________________________________%
% Initial Input:
%___________________________________________________________________________________%

count = 0;

qq = 0;

for i = 1 : N
    
count = count + 1; 
    
    p(count,1) = 0.125 * exp(-0.1*qq);
    
    s(count,1) = gamma * (1 - 0.01*exp(-0.1*qq));
    
    qq = qq + DX;
    
end;   

x = [p;zeros(N,1);zeros(N,1);s;ones(N,1)];

%___________________________________________________________________________________%                                 
% Using Explicit Euler Method 
%___________________________________________________________________________________%

counter = 0;

for T = t1 : h : (t2 - h)
        
counter = counter + 1;

x = x + (myfun(x,DX,Dc,alpha,gamma,k1,k2,k3,c0,N) * h);

ProliferatingCells = x(1 : N,1);                                % Calculated values: ProliferatingCells

QuiescentCells = x(1 + N : (2 * N),1);                          % Calculated values: QuiescentCells

NecroticCells = x(1 + 2 * N : (3 * N),1);                       % Calculated values: NecroticCells

SurroundingCells = x(1 + 3 * N : (4 * N),1);                    % Calculated values: SurroundingCells

Nutrients = x(1 + 4 * N : (5 * N),1);                           % Calculated values: Nutrients

%Store values

% at the time point t = 7.5
if (counter == (7.5/h))
    
    P75 = ProliferatingCells;
    Q75 = QuiescentCells;
    Ne75 = NecroticCells;
    S75 = SurroundingCells;
    Nu75 = Nutrients;
    
end

% at the time point t = 15
if (counter == (15/h))
    
    P150 = ProliferatingCells;
    Q150 = QuiescentCells;
    Ne150 = NecroticCells;
    S150 = SurroundingCells;
    Nu150 = Nutrients;
    
end

% at the time point t = 22.5
if (counter == (22.5/h))
    
    P225 = ProliferatingCells;
    Q225 = QuiescentCells;
    Ne225 = NecroticCells;
    S225 = SurroundingCells;
    Nu225 = Nutrients;
    
end

% at the time point t = 30
if (counter == (30/h))
    
    P300 = ProliferatingCells;
    Q300 = QuiescentCells;
    Ne300 = NecroticCells;
    S300 = SurroundingCells;
    Nu300 = Nutrients;
    
end

% at the time point t = 37.5
if (counter == (37.5/h))
    
    P375 = ProliferatingCells;
    Q375 = QuiescentCells;
    Ne375 = NecroticCells;
    S375 = SurroundingCells;
    Nu375 = Nutrients;
    
end

% at the time point t = 45
if (counter == (45/h))
    
    P450 = ProliferatingCells;
    Q450 = QuiescentCells;
    Ne450 = NecroticCells;
    S450 = SurroundingCells;
    Nu450 = Nutrients;
    
end

% at the time point t = 52.5
if (counter == (52.5/h))
    
    P525 = ProliferatingCells;
    Q525 = QuiescentCells;
    Ne525 = NecroticCells;
    S525 = SurroundingCells;
    Nu525 = Nutrients;
    
end

end;

figure ; plot(XX,P75,XX,P150,XX,P225,XX,P300,XX,P375,XX,P450,XX,P525) ; grid on ; title('Enlarged model of an avascular tumour growth') ; xlabel('Space x, cm') ; ylabel('Proliferating Cells, p');

figure ; plot(XX,Q75,XX,Q150,XX,Q225,XX,Q300,XX,Q375,XX,Q450,XX,Q525) ; grid on ; title('Enlarged model of an avascular tumour growth') ; xlabel('Space x, cm') ; ylabel('Quiescent Cells, q');

figure ; plot(XX,Ne75,XX,Ne150,XX,Ne225,XX,Ne300,XX,Ne375,XX,Ne450,XX,Ne525) ; grid on ; title('Enlarged model of an avascular tumour growth') ; xlabel('Space x, cm') ; ylabel('Necrotic Cells, n');

figure ; plot(XX,S75,XX,S150,XX,S225,XX,S300,XX,S375,XX,S450,XX,S525) ; grid on ; title('Enlarged model of an avascular tumour growth') ; xlabel('Space x, cm') ; ylabel('Surrounding Tissue Cells, s');

figure ; plot(XX,Nu75,XX,Nu150,XX,Nu225,XX,Nu300,XX,Nu375,XX,Nu450,XX,Nu525) ; grid on ; title('Enlarged model of an avascular tumour growth') ; xlabel('Space x, cm') ; ylabel('Nutrients, c');



