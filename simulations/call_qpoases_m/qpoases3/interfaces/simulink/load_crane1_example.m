%%
%%	This file is part of qpOASES.
%%
%%	qpOASES -- An Implementation of the Online Active Set Strategy.
%%	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
%%	Christian Kirches et al. All rights reserved.
%%
%%	qpOASES is free software; you can redistribute it and/or
%%	modify it under the terms of the GNU Lesser General Public
%%	License as published by the Free Software Foundation; either
%%	version 2.1 of the License, or (at your option) any later version.
%%
%%	qpOASES is distributed in the hope that it will be useful,
%%	but WITHOUT ANY WARRANTY; without even the implied warranty of
%%	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%%	See the GNU Lesser General Public License for more details.
%%
%%	You should have received a copy of the GNU Lesser General Public
%%	License along with qpOASES; if not, write to the Free Software
%%	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
%%



%%
%%	Filename:  interfaces/simulink/load_example_QProblem.m
%%	Author:    Aude Perrin, Hans Joachim Ferreau
%%	Version:   3.0embedded
%%	Date:      2007-2014
%%


clear all;

load 'benchmarkCRANE1.mat';

[nC,nV] = size(A);
[nV,nP] = size(g);
    
H_data = H;
g_data = g;
A_data = A;
lb_data = lb;
ub_data = ub;
lbA_data = lbA;
ubA_data = ubA;

clear H g A lb ub lbA ubA;


%% setup QP data
simulationTime = ( 0:(nP-1) )' * 0.1;

H.time = simulationTime;
H.signals.values = repmat( H_data(:)',nP,1 );
H.signals.dimensions = length( H_data(:)' );

g.time = simulationTime;
g.signals.values = g_data';
g.signals.dimensions = size(g_data,1);

A.time = simulationTime;
A_tmp = A_data';
A.signals.values = repmat( A_tmp(:)',nP,1 );
A.signals.dimensions = length( A_data(:)' );

lb.time = simulationTime;
lb.signals.values = lb_data';
lb.signals.dimensions = size(lb_data,1);

ub.time = simulationTime;
ub.signals.values = ub_data';
ub.signals.dimensions = size(ub_data,1);

lbA.time = simulationTime;
lbA.signals.values = lbA_data';
lbA.signals.dimensions = size(lbA_data,1);

ubA.time = simulationTime;
ubA.signals.values = ubA_data';
ubA.signals.dimensions = size(ubA_data,1);


%% open corresponding simulink example
open( 'example_QProblem.mdl' );



%%
%%	end of file
%%
