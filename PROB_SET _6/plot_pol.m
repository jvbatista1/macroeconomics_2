function PF = plot_pol(r)
% Function to Plot Household Policy Functions (for highest and lowest
% employment states) given an interest rate r
% PF is the policy function matrix used for subsequent simulation

% Redefine the Global Variables from the Master File
global beta mu theta delta l_grid Pi Na Nl a_grid b % rtp lnl_grid Ni T Tburn 

%% Calculate the Implied Wage for a Given Interest Rate
% =======================================================================
% Because I only solve for this once, I don't bother defining a function
wage = (1-theta)*(theta/(r+delta)^(theta/(1-theta)));

%% Define the Debt Limit
% =======================================================================
% Use the b defined before the function is called
phi = -b;

%% Form the Household Asset Holding Grid
% =======================================================================
a_max   = 20;       % Maximum Capital Holdings (needed for computer to solve)
a_min   = phi;      % Minumum Capital Holdings (mathematically defined in problem)

% TRICK: Define the wealth grid logarithmically to put more grid points at
% areas of the state space that are most non-linear (i.e., around the
% borrowing constraint)
a_grid  = (exp(linspace(log(a_min+1-a_min),log(a_max+1-a_min),Na))-1+a_min)'; % Grid for wealth (endogenous state)

% If, instead, you want a vanilla equally-spaced grid
% a_grid = linspace(a_min,a_max,Na);

%% Calculate the Household Decision Rules
% =======================================================================
% I use Value Function Iteration (VFI) to do this.
% For a complete exposition of VFI, see my notes and code for Ph.D. 21,
% Problem Set 3 on the Faculty Website.

% I carry out the VFI in two stages:
% (1) I build a 3-dimensional grid for contemporaneous utility, where
% dimension 1 refers to today's asset holdings, dimension 2 to
% tomorrow's asset holdings and dimension 3 to today's labour endowment
% shock.
% (2) Using the output from (1), I run a vanilla value function
% iteration.

% This approach to VFI is advantageous as complex inequality
% constraints can be easily applied in stage (1) by setting
% contemporaneous utility to arbitrarily large negative values at grid
% points that violate constraints.

%% Build a 3-Dimensional Contemporaneous Utility Grid
% ===================================================================
Ut = zeros(Na,Na,Nl);   % Dimension 1: Today's assets (a);
                        % Dimension 2: Tomorrow's assets (a');
                        % Dimension 3: Today's labour endowment (l);

% Household Utility
u = @(c) (c^(1-mu)-1)/(1-mu);

for kk = 1:Nl            % Loop Over Labour Supply Endowment Today
    for ii = 1:Na        % Loop Over Savings Today
        for jj = Na:-1:1 % Loop Over Savings Tomorrow
            l  = l_grid(kk);     % Labour Supply Today
            a  = a_grid(ii);     % Savings Today
            ap = a_grid(jj);     % Savings Tomorrow
            % Solve for Consumption at Each Point
            c = wage*l + (1+r)*a - ap;

            % Apply Constraints
            if (ap < phi) || (c < 0)
                % If Tomorrow's Savings Violate Debt Limit or Today's
                % Consumption is Negative, set Utility to be arbitrarily
                % negative
                Ut(ii,jj,kk) = -99999999999999;
            else
                % If constraints satisfied, calculate contemporaneous
                % utility
                Ut(ii,jj,kk) = u(c);
            end
        end
    end
end

%% Vanilla Value Function Iteration
% ===================================================================
% Initial Guess of Value Function
V0 = kron(l_grid,ones(Na,1));      % l_grid is Nl x 1; ones(Na,1) is Na x 1
                                   % V0 is Na x Nl
EVF = zeros(max(size(V0)),min(size(V0)));
V1 = zeros(max(size(V0)),min(size(V0)));
PF = zeros(max(size(V0)),min(size(V0))); % policy function

% VFI Tolerance and Convergence Criteria
vfi_tol  = 0.0001;
vfi_err  = 2;
vfi_iter = 0;

while vfi_err > vfi_tol

    % Calculate the Guess of the Expected Value Function
    for kk = 1:Nl
        % Loop Over Today's Labour Endowment to Solve for Expected Value
        % Function
        EVF(:,kk) = V1*Pi(kk,:)';
        % To be used as the next guess for the value function
    end
    
    for kk = 1:Nl   % Loop Over Today's Labour Endowment Shock
        % Here, I use Kronecker products to save a loop over asset holdings 
        % for today and tomorrow
        %
        % See that, for kk=1, Ut(:,:,1) is a Na x Na matrix
        %
        Objgrid = Ut(:,:,kk) + beta*kron(ones(Na,1),EVF(:,kk)'); % Kron gives Na x Na
        for ii = 1:Na
            % Loop Over Today's Capital Stock to Find Max of Value Function
            [V1(ii,kk),PF(ii,kk)] = max(Objgrid(ii,:));
        end
    end

    vfi_iter = vfi_iter + 1;
    vfi_err  = norm(V1(:) - V0(:));
    iter100 = mod(vfi_iter,50);
    if iter100 == 0
        display(['Number of VF Iterations ',num2str(vfi_iter)]);
    end
    V0 = V1; % update
end

%% Build Household Policy Functions for Saving and Consumption
% ===================================================================
% Policy Function for Assets
AF  = nan(size(PF));
% Policy Function for Consumption
CF  = nan(size(PF));

for kk = 1:Nl
    for  ii = 1:Na
        a  = a_grid(ii);        % Savings Brought into Today
        ap = a_grid(PF(ii,kk)); % Savings Made Today for Tomorrow
        AF(ii,kk) = ap;
        l  = l_grid(kk);
        CF(ii,kk) = wage*l + (1+r)*a - ap; % Consumption possibilities
    end
end

% Plot Policy Function for ap against a
fontsize = 14;
fig1 = figure('units','normalized');
set(fig1,'Color','white','numbertitle','off','name','Policy Function - Savings')
plot(a_grid,AF(:,1),'k-.','LineWidth',1); hold on
plot(a_grid,AF(:,end),'k','LineWidth',2); hold on
plot(a_grid,a_grid,'k:','LineWidth',1); hold off
legend({'Low Employment State','High Employment State','45-degree Line'},...
    'FontSize',fontsize,'Location','SouthEast','Orientation','Vertical','Interpreter','latex')
title('Household Savings Policy Function','FontSize',fontsize,'Interpreter','latex')
xlabel('$a_{t}$','FontSize',fontsize,'Interpreter','latex')
ylabel('$a_{t+1}$','FontSize',fontsize,'Interpreter','latex')
axis('tight')
grid on, grid minor
