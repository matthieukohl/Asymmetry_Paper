function [par] = PARS(parstr, options);
% PARS  Parameters for plots and analyses.

  %%%%%%%%%%%%%           parse options         %%%%%%%%%%%%%%%%%%%%
  error(nargchk(1, 2, nargin))     
  if nargin ==1 | isempty(options)
    fopts      = [];
  else
    fopts      = fieldnames(options);
  end

  if strmatch('slide', fopts)
    slide = options.slide;
  else 
    slide = 1;
  end


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  
  % flag for file dump of plots, and directory for file dumps
  plot_file           = 0;
  plot_dir            = '/home/pog/www/'; 
    
  % plot parameters
  if (slide)
   % size of figures (in inches)
    figure_width      = 16/2.54;
    figure_height     = figure_width / 1.666;

    % line widths
    line_width        = 1.2;
    thin_line_width   = .4;
    thick_line_width  = 1.6;
    axes_line_width   = 0.8;

    % line color
    pos_line_color    = 'r';
    neg_line_color    = 'b';
    line_color        = 'b';

    % tick length
    tick_length       = [0.01 0.01];


    % line styles
    line_style        = '-';
    thin_line_cstyle  = 'm-';
    thick_line_cstyle = 'r-';

    % marker size
    marker_size       = 10;

    % font sizes
    font_size         = 14;
    clabel_size       = 12;

  else
    % size of figures (in inches)

    % Science 
    %figure_width      = 2.25; % one column
    %figure_width      = 4.75; % two column
    %figure_height     = 2.25;
    
    % Princeton UP
    %figure_width      = 4.5;  % PUP page
    %figure_width      = 2.5;
    %figure_width      = 3.5;
    %figure_height     = figure_width / 1.666;
    
    % AMS journals
    figure_width      = 3.125;
    figure_height     = figure_width / 1.666;
    %figure_height     = figure_width / 1.618 % tapio wv review

    % annual review
    %figure_width      = 3.38;
    %figure_height     = figure_width/1.6;


    % supercriticality figure in JAS
    %figure_width      = 3;
    %figure_height     = figure_width;
    
    % line colors and styles	
    pos_line_color    = 'k';
    neg_line_color    = 'k';
    line_color        = 'k';
    line_style        = '-';
    thin_line_cstyle  = 'k:';
    thick_line_cstyle = 'k-';


    % line widths
    line_width        = .5;  % pt 
    thin_line_width   = .5; % PUP
    thick_line_width  = 1.0;
    axes_line_width   = 1.0;

    % marker size
    marker_size       = 5;

    % font sizes
    font_size         = 8;
    %clabel_size       = 6; % PUP figures
    clabel_size       = 7;

    % tick length
    tick_length       = [0.015 0.015];

    % used for hydro paper
    %axes_line_width=0.7
    %tick_length=[0.025 0.025]

  end

  
  % physical parameters
  radius_earth 	      = 6.37122e6;    % in meters
  gravity 	      = 9.80665;      % in m / s^2

  % choose the following for consistency with the serial analysis code
  gas_constant	      = 287.04;       % J / kg / K 
  Rv                  = gas_constant/0.622;

  kappa               = 2.0/7.0;
  cp                  = gas_constant/kappa;  % specific heat at constant pressure
  cp_v                = 1870;         % specific heat water vapor [J/kg/K]

  cl                  = 4190;         % heat capacity liquid water [J/kg/K]
  gas_constant_v      = 461.50;       % gas constant water vapor [J/kg/K]
  latent_heat_v       = 2.5e6;        % latent heat of evaporation [J/kg]



  omega               = 7.292e-5;     % angular velocity earth [rad / s]
                                      % related to sidereal day (not 24hrs)
  mean_sfc_press      = 1e5;          % mean surface pressure
  secs_per_day        = 86400;        % length of day in seconds
  l_cond              = 2.5e6;        % latent heat cond [J/kg] at zero deg C
  deg                 = pi/180;
  stefan              = 5.67e-8; 


  %sfc_cdf_isolines   = [.5 .5];      % for plotting of median only
  sfc_cdf_isolines    = [.1 .5 .9];  % isolines of scf temp distr to be plotted
  
  eval(['par = ', parstr, ';']);
  
  
