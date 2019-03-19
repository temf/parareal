%Plot reference results

% The source code is distributed under the terms of the GNU General
% Public License (GPL) (version 2 or later).

% read solution from file
ref_sol  = load ('solution.csv');
ref_err  = load ('error.csv');
TimeStep = ref_sol(2,1)-ref_sol(1,1);

figure; 
plot(ref_sol(:,1),ref_sol(:,2),...
     ref_sol(:,1),ref_sol(:,3)); % plot 1st component of solutions
title('Solution');
grid
axis tight
legend(['sequential ' num2str(TimeStep,'%1.2e')],'parareal')
xlabel('time / s')
ylabel('magnetic vector potential / Wb')

figure; 
semilogy(ref_err(:,1),ref_err(:,2),'*-');
title('Convergence');
xlabel('iteration')
ylabel('error')