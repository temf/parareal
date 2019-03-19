%% Parareal solution for an induction machine model im_3kW
%
% The algorithm is documented in the paper
% "A New Parareal Algorithm for Problems with Discontinuous Sources"
% by
% I. Kulchytska-Ruchka, S. Schöps, Technische Universitaet Darmstadt
% I. Niyonzima, Université Grenoble Alpes
% M. J. Gander, Université de Genève
% see https://arxiv.org/abs/1803.05503 
%
% The model is originally developed by J. Gyselinck, R.V. Sabariego
% see http://onelab.info

% The source code is distributed under the terms of the GNU General
% Public License (GPL) (version 2 or later).

% path to 'getdp', hardcode if necessary
% e.g. getdp_path = '/usr/bin/getdp/getdp';
getdp_path = getdp_findexe();

%% A. Set-up
% A.1. Fine model
ode_f    = strcat(pwd,'/model/model_f/im_3kW.pro');
solver_f = 'odegetdp';
opt_f    = struct('Verbose',0,'TimeStep',1e-6,...
                  'Executable',getdp_path);
% A.2. Coarse model
ode_c    = strcat(pwd,'/model/model_c/im_3kW.pro');
solver_c = 'odegetdp';
opt_c    = struct('Verbose',0,'TimeStep',1e-3,...
                  'Executable',getdp_path);
% A.3. Parameters for Parareal
opt_para   = struct('NProc',20,'NIter',20,'Verbose','on',...
                  'OutputSel',1,'OutputFcn',[],...
                  'AbsTol',1.5e-5,'RelTol',1.5e-5);
% A.4. Parameters for GetDP
opt_getdp  = {'Flag_AnalysisType',1,...
              'Flag_NL',1,...
              'Flag_ImposedSpeed',1,...
              'Nb_max_iter',60,...
              'relaxation_factor', 0.5,...
              'stop_criterion', 1e-6,...
              'NbTrelax', 2,...
              'Flag_SaveAllSteps', 1};
% A.5. Interval and initial value
trange   = [0 0.02];
numdof   = getdp_numdof (ode_c,trange,opt_c,...
                         'Flag_AnalysisType',1,...
                         'Flag_NL',1,...
                         'Flag_ImposedSpeed',1);
init     = zeros(numdof,1);

%% B. Parareal solution
sol_para = odeparareal (ode_c,solver_c,opt_c,...
                        ode_f,solver_f,opt_f,...
                        trange,init,opt_para,opt_getdp{:});

%% C. Sequential solution
sol_seq  = odegetdp(ode_f,trange,init,opt_f,opt_getdp{:});

%% D. Visualization
figure; 
plot(sol_seq.x,sol_seq.y(1,:),...
     sol_para.x,sol_para.y(1,:)); % plot 1st component of solutions
title('Solution');
xticks(trange(1):opt_c.TimeStep:trange(end))
xtickangle(45)
grid
axis tight
legend(['sequential ' num2str(opt_f.TimeStep,'%1.2e')],'parareal')
xlabel('time / s')
ylabel('magnetic vector potential / Wb')

figure; 
semilogy(1:length(sol_para.err)-1,sol_para.err(2:end),'*-');
title('Convergence');
xlabel('iteration')
ylabel('error')