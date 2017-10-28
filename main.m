clc
clear all
close all

%psri = 'yes';
psri = 'no';

% Start time
tic


disp(['Main (PSRI=' psri ')']);

% out = 'yes';
out = 'no';

% plot = 'yes';
plot = 'no';

% Mu Definition
mu = 1;

% Lambda Definition
lambda = 1;

% h_max Array Definition
h = [0.4, 0.3, 0.2, 0.15, 0.1, 0.075];
n = length(h);

% Thickness Parameter Array Definition
t = [0.00001, 0.0001, 0.001, 0.01];
nt = length(t);

% Mesh Building
% meshbuilder(h);

% Quadrature Degree
fdq = 'degree=4';

% Error Arrays
errL2 = zeros(n,nt);
errH1 = zeros(n,nt);

% Computation
for j=1:nt
    
    disp(['t = ', num2str(t(j))]);
    
    pre_comp_sym (mu, lambda, t(j));
    
    for i=1:n
    
      %%% Cleaning WorkSpace
	  clear xv yv vertices edges endpoints boundary boundedges;
    
      %%% Old Mesh Building
      % [xv,yv,vertices,edges,endpoints,boundary,boundedges] = meshgen(h(i));
    
      %%% Meshes Loading
      meshname = ['mesh' num2str(h(i)) '.mat'];
      load(meshname);
      disp(['    Mesh: ', num2str(i), '   - h_max = ', num2str(h(i)), ' loaded:']);
      nver = length(xv);
      nedge = length(endpoints);
      nele = length(edges);    
    
    
      % Error Computation
      [errL2(i,j), errH1(i,j)] = RM (...
                                     xv, yv, ...
                                     vertices, ...
                                     edges, ...
                                     endpoints, ...
                                     boundary, ...
                                     boundedges, ...
                                     fdq, ...
                                     mu, ...
                                     lambda, ...
                                     t(j), ...
                                     psri, ...
                                     h(i), ...
                                     out);
                                       
      disp(['    - errL2 = ', num2str(errL2(i,j))]);
      disp(['    - errH1 = ', num2str(errH1(i,j))]);
                                   
    end % end for i
    
end % end for j

clc

wtime = toc;
fprintf ( 1, '  MAIN took %f seconds to run.\n', wtime );

disp(['h = ', num2str(h)]);
disp(' ');

% Plot Loop
 for j=1:nt
     
  %disp(['t = ', num2str(t(j))]); 
  %disp(['     errL2 = ', num2str(errL2(:,j)')]);
  %disp(['     errH1 = ', num2str(errH1(:,j)')]); 
    
% Figure 1, L2 Error
  figure(1);
  hold on;
 
    subplot(2, 2, j);
    loglog (h, errL2(:,j), '-*r', h, h.^3,'-b');
    grid on;
    % legend ('ErrL2', 'h^3', 'location', 'northeastoutside');
    title (['ErrL2 (t = ', num2str(t(j)), ')']);
  
  if (strcmp(psri,'yes'))
      saveas (1, 'psri_convL2.png');
  else
      saveas (1, 'fem_convL2.png');
  end
  
%  Figure 2, H1 Error
  figure(2);
 
  hold on;
  
      subplot(2, 2, j), loglog (h, errH1(:,j), '-*r', h, h.^2,'-b');
      grid on;
      % legend ('ErrH1', 'h^2', 'location', 'northeastoutside');
      title (['ErrH1 (t = ', num2str(t(j)), ')']);

  if (strcmp(psri,'yes'))
      saveas (2, 'psri_convH1.png');
  else
      saveas (2, 'fem_convH1.png');
  end
  
     
  end

  
% LATEX ERROR TABLE 
 latex = 'yes';
 
if (strcmp(latex,'yes'))
    
    if (strcmp(psri,'yes'))
    
        latexReport('psri_latex_report.tex', errL2, errH1, h, nt, n, t);
        
    else
        
        latexReport('fem_latex_report.tex', errL2, errH1, h, nt, n, t);
        
    end
    
end % end if

