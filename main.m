function main(nt, taut)

if nargin<1
    nt=10; % severity of change
end

if nargin<2
    taut=10; % frequency of change
end

clc
close all
format compact

warning('off', 'all');

path =fileparts(which('main.m'));
addpath(genpath(path));

% create Results folder to save data (including HV, IGD,SP,Runtime)
filename=append('Results_nt', num2str(nt), '_taut',num2str(taut));
path=append(path, filesep, filename);
check_create_folder(path);

set_global(nt, taut); % set global variables. [caution: do not use it in parallel computing]

runmax=1; % max number of independent runs

%list of algorithms available to test: ["VARE", "VARES", "SGEA", "TrRMMEDA", "PPS","MOEADSVR"]; 
%algorithms=["VAREB", "SGEA", "PPS", "TrRMMEDA", "MOEADSVR"]; % algorithms to experiment
algorithms=["VARE"];%, "VARE_MOEAD"];
prob_ids=212:215; % problems to solve: DF(101-114), SDP(201-215, excluding 212/213: as SDP12/13 not supported), FDA(301-305), DMOP(401-403), F(501-510) 

% initialise varibles as []
[igd, hv, sp, runtime]=deal([]);

for i =1:length(algorithms)
    % if i==5 || i==6, runmax=5; end  % reduce number of runs for slow algorithms
    algo_name=algorithms(i);

    % create a folder for the chosen algorithm to save results
    result_path=append(path, filesep, algo_name);
    check_create_folder(result_path);
    %addpath(result_path);

    % run the chosen algorithm
    [ArchiveResult, IGDResult, SPResult, HVResult, TimeResult]=run_algorithm(algo_name,prob_ids,runmax);

    % save individual data; append - combine strings
    parsave(append(result_path, filesep,'IGD.mat'), IGDResult);
    parsave(append(result_path, filesep,'SP.mat'), SPResult);
    parsave(append(result_path, filesep,'HV.mat'), HVResult);
    parsave(append(result_path, filesep,'Time.mat'), TimeResult);
    parsave(append(result_path, filesep,'Pops.mat'), ArchiveResult);

    % concatenate results accross algorithms
    sep_col=ones(length(prob_ids),1);
    igd=[igd,[IGDResult, sep_col]];
    hv=[hv, [HVResult, sep_col]]; 
    sp=[sp, [SPResult, sep_col]]; 
    runtime=[runtime, [TimeResult,sep_col]];
end

% save each metric data in a separate file in order with timestamp
alg_order=strjoin(algorithms,'_');
date_str=datestr(now,'mm-dd-yyyy HH-MM');

writematrix(igd,append(filename,filesep,sprintf('igd_%s_%s.csv', alg_order, date_str)));
writematrix(hv,append(filename,filesep,sprintf('hv_%s_%s.csv', alg_order, date_str)));
writematrix(sp,append(filename,filesep,sprintf('sp_%s_%s.csv', alg_order, date_str)));
writematrix(runtime,append(filename,filesep,sprintf('time_%s_%s.csv', alg_order, date_str)));

end

function [ArchiveResult, IGDResult, SPResult, HVResult, TimeResult]=run_algorithm(algo_name,prob_ids, max_runs)
fh=str2func(algo_name);

count=1;
for probID=prob_ids
    % set population size according to num_objs
    prob=problem(probID);
    M=prob.objDim;
    pop_size =100;
    if M==3, pop_size=105; end % popsize of 105 is used in tri-objective case

    [Archives,Time, IGDs, HVs, SPs]=fh(probID,pop_size, max_runs); % run each algorithm

    IGDResult(count,:)=stats(IGDs);
    SPResult(count,:)=stats(SPs);
    HVResult(count,:)=stats(HVs);
    TimeResult(count,:)=[min(Time) max(Time) median(Time) mean(Time) std(Time)];
    ArchiveResult{count}=Archives;
    count =count+1;
end
end

function out=stats(in)
avg=mean(in,2); % MIGD for each run
out=[min(avg) max(avg) median(avg) mean(avg) std(avg)];
end

% saving results in parallel environments
function parsave(fname, x)
save(fname, 'x')
end

function check_create_folder(folderName)
if ~isfolder(folderName)
    mkdir(folderName)
end
end

function set_global(in_nt, in_taut)
global g_nt g_taut
g_nt=in_nt;
g_taut=in_taut;
end