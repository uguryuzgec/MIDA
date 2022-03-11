clear all
close all
clc
% mex cec14_func.cpp -DWINDOWS
% 17-22 hybrid func. 
% 1-3 unimodal func.
% 4-16 multimodal func.
% 23-28 composition func.
func_num=25; 
runs=1;
D=2; % dimension
Xmin=-100;
Xmax=100;
pop_size=10*D;
iter_max=1000;
fhd=str2func('cec14_func');
empty_solution.cost=[];
empty_solution.position=[];
empty_solution.t=[];
solution=repmat(empty_solution,func_num,runs);
 lb=Xmin;
 ub=Xmax;
 X_suru=lb+(ub-lb).*rand(pop_size,D);

[gbest,gbestval,FES,t]= DA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,X_suru,func_num);

solution(func_num,1).position = gbest;
solution(func_num,1).cost = gbestval-func_num*100;
solution(func_num,1).t = t;
fprintf('---------------------------------------------------------------\n');
fprintf('Optimization with DA \n');
fprintf('Func no: %d -> %d. run : best error = %1.2e\n',func_num,runs,solution(func_num,1).cost);

plot_function
drawnow
% 1. Convergence Analysis
fonk_numara = func_num;
FESindex = [0.01, 0.02, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]*FES;

% % % figure(2)
% % % for j=1:runs
% % % 	semilogy(FESindex,solution(fonk_numara,j).t,'-s','color',rand(1,3));
% % % 	hold on
% % % end
% % % xlabel('Function Evaluations');
% % % ylabel('Error Value');
% % % str = sprintf('Convergence Analysis of FN%d',fonk_numara);
% % % title(str);
% saveas(figure (1),'analiz\fonk.jpeg')

% testing your algorithm for all functions with xxx runs 
% % % for i=1:func_num
% % %     func_num=i;
% % %     for j=1:runs
% % %          [gbest,gbestval,FES,t]= PSO_func(fhd,D,pop_size,iter_max,Xmin,Xmax,func_num);
% % %          solution(i,j).position = gbest;
% % % 		 solution(i,j).cost = gbestval-func_num*100;
% % % 		 solution(i,j).t = t;
% % % 		 fprintf('Func no: %d -> %d. run : best error = %1.2e\n',i,j,solution(i,j).cost);
% % %     end
% % %     fprintf('--------------------------------------------\n')
% % % 	name = strcat('PSO_',num2str(func_num),'_',num2str(D),'.txt'); 
% % %     write([solution(i,:).t],name); 
% % % end

if runs >1
	plot_statistic(solution,func_num,runs);
end

figure (2)
semilogy(FESindex,solution(fonk_numara,1).t,'-db');
xlabel('Function Evaluations');
ylabel('Error Value');
str = sprintf('Convergence Analysis of FN%d',fonk_numara);
title(str);

% Optimization with MIDA 

[i_gbest,i_gbestval,i_FES,i_t]= ImpDA_func(fhd,D,pop_size,iter_max,Xmin,Xmax,X_suru,func_num);

solution(func_num,1).position = i_gbest;
solution(func_num,1).cost = i_gbestval-func_num*100;
solution(func_num,1).t = i_t;
fprintf('---------------------------------------------------------------\n');
fprintf('Optimization with ImpDA\n ');
fprintf('Func no: %d -> %d. run : best error = %1.2e\n',func_num,runs,solution(func_num,1).cost);

% 1. Convergence Analysis
figure (2)
hold on
semilogy(FESindex,solution(fonk_numara,1).t,'-sr');
legend('DA','MIDA')

% % % figure 
% % % for j=1:runs
% % % 	semilogy(i_FESindex,solution(fonk_numara,j).t,'-s','color',rand(1,3));
% % % 	hold on
% % % end
% % % xlabel('Function Evaluations');
% % % ylabel('Error Value');
% % % str = sprintf('Convergence Analysis of FN%d',fonk_numara);
% % % title(str);

if runs >1
	plot_statistic(solution,func_num,runs);
end
