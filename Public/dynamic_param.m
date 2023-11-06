function [nt, taut, T0]=dynamic_param()
global g_nt g_taut
% set dynamic parameters here


if isempty(g_nt) || isempty(g_taut)
    warning('taut or nt does not a value, and they will be initialised as 10.')
    g_nt=10;
    g_taut=10;
end

nt=g_nt;
taut=g_taut;

T0=100;
end