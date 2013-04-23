function out = release(Co, Avd, Bvd, Cvd, Dvd, D01, Kr, Kf, dens, hp, ht, tf, dt)
%
% model.m
%
% Model exported on Apr 15 2013, 22:17 by COMSOL 4.3.0.151.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/users/che/bt/ankeshs/Academics/che496/Transient');

model.name('polycoat3.mph');

model.param.set('Re_cyl', '1');
model.param.set('R_cyl', '.01');
model.param.set('a_ratio', num2str(ht));
model.param.set('rho_fluid', '1080');
model.param.set('mu_fluid', '0.0089');
model.param.set('b_conc', '1');
model.param.set('D_drug_fluid', '3e-8');
model.param.set('t_poly', num2str(hp));
model.param.set('D_drug_polyI', '1e-15');
model.param.set('K_part', num2str(Kf));
model.param.set('D_drug_polyO', '1e-15*0.5');
model.param.set('conc_init', num2str(Co));
model.param.set('time_final', num2str(tf));
model.param.set('time_step', num2str(dt));
model.param.set('K_r', num2str(Kr));
model.param.set('M_drug', '0.853906');
model.param.set('rho_poly', num2str(dens));
model.param.set('A_vd', num2str(Avd));
model.param.set('B_vd', num2str(Bvd));
model.param.set('C_vd', num2str(Cvd));
model.param.set('D_vd', num2str(Dvd));
model.param.set('D_pex', num2str(D01));
model.param.set('sStep', '20');
model.param.set('lStep', '1000');
model.param.set('sEnd', '1000');

model.modelNode.create('mod1');

model.geom.create('geom1', 2);
model.geom('geom1').axisymmetric(true);
model.geom('geom1').lengthUnit('mm');
model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature.create('c1', 'Circle');
model.geom('geom1').feature.create('sq1', 'Square');
model.geom('geom1').feature.create('int1', 'Intersection');
model.geom('geom1').feature.create('c2', 'Circle');
model.geom('geom1').feature('r1').set('pos', {'0' '-3*R_cyl'});
model.geom('geom1').feature('r1').set('size', {'R_cyl' '5*R_cyl'});
model.geom('geom1').feature('c1').set('pos', {'0' '0'});
model.geom('geom1').feature('c1').set('r', 'a_ratio*R_cyl');
model.geom('geom1').feature('sq1').set('pos', {'0.2' '0'});
model.geom('geom1').feature('sq1').set('base', 'center');
model.geom('geom1').feature('sq1').set('size', '0.4');
model.geom('geom1').feature('int1').selection('input').set({'c1' 'sq1'});
model.geom('geom1').feature('c2').set('rot', '270');
model.geom('geom1').feature('c2').set('r', '(1-t_poly)*a_ratio*R_cyl');
model.geom('geom1').feature('c2').set('angle', '180');
model.geom('geom1').run;

model.view.create('view2', 3);
model.view.create('view3', 3);
model.view.create('view4', 3);
model.view.create('view5', 3);

model.physics.create('spf', 'LaminarFlow', 'geom1');
model.physics('spf').selection.set([1]);
model.physics('spf').feature.create('inl1', 'Inlet', 1);
model.physics('spf').feature('inl1').selection.set([2]);
model.physics('spf').feature.create('out1', 'Outlet', 1);
model.physics('spf').feature('out1').selection.set([8]);
model.physics.create('chds', 'DilutedSpecies', 'geom1');
model.physics('chds').selection.set([1]);
model.physics('chds').feature.create('in1', 'Inflow', 1);
model.physics('chds').feature('in1').selection.set([2]);
model.physics('chds').feature.create('out1', 'Outflow', 1);
model.physics('chds').feature('out1').selection.set([8]);
model.physics('chds').feature.create('conc1', 'Concentration', 1);
model.physics('chds').feature('conc1').selection.set([10 13]);
model.physics.create('chds2', 'DilutedSpecies', 'geom1');
model.physics('chds2').selection.set([2]);
model.physics('chds2').feature.create('conc2', 'Concentration', 1);
model.physics('chds2').feature('conc2').selection.set([11 12]);
model.physics('chds2').feature.create('conc3', 'Concentration', 1);
model.physics('chds2').feature('conc3').selection.set([10 13]);
model.physics.create('chds3', 'DilutedSpecies', 'geom1');
model.physics('chds3').selection.set([3]);
model.physics('chds3').feature.create('conc1', 'Concentration', 1);
model.physics('chds3').feature('conc1').selection.set([11 12]);
model.physics('chds3').feature.create('open1', 'OpenBoundary', 1);
model.physics('chds3').feature('open1').selection.set([11 12]);

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('size1', 'Size');
model.mesh('mesh1').feature('size1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('size1').selection.set([1]);
model.mesh('mesh1').feature.create('size2', 'Size');
model.mesh('mesh1').feature('size2').selection.geom('geom1', 1);
model.mesh('mesh1').feature('size2').selection.set([9 10 13]);
model.mesh('mesh1').feature.create('ftri1', 'FreeTri');
model.mesh('mesh1').feature('ftri1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('ftri1').selection.set([1]);
model.mesh('mesh1').feature.create('bl1', 'BndLayer');
model.mesh('mesh1').feature('bl1').selection.geom('geom1', 2);
model.mesh('mesh1').feature('bl1').selection.set([1]);
model.mesh('mesh1').feature('bl1').feature.create('blp1', 'BndLayerProp');
model.mesh('mesh1').feature('bl1').feature('blp1').selection.set([9 10 13]);
model.mesh('mesh1').feature.create('ftri2', 'FreeTri');

model.view('view1').axis.set('xmin', '-0.042106471955776215');
model.view('view1').axis.set('xmax', '0.05210646986961365');
model.view('view1').axis.set('ymin', '-0.04628787189722061');
model.view('view1').axis.set('ymax', '0.03659117966890335');

model.physics('spf').prop('CompressibilityProperty').set('Compressibility', 'Incompressible');
model.physics('spf').prop('PseudoTimeProperty').set('locCFL', '1.3^min(niterCMP-1,9)+if(niterCMP>25,9*1.3^min(niterCMP-25,9),0)+if(niterCMP>50,90*1.3^min(niterCMP-50,9),0)');
model.physics('spf').feature('fp1').set('rho_mat', 'userdef');
model.physics('spf').feature('fp1').set('rho', 'rho_fluid');
model.physics('spf').feature('fp1').set('mu_mat', 'userdef');
model.physics('spf').feature('fp1').set('mu', 'mu_fluid');
model.physics('spf').feature('fp1').set('minput_velocity_src', 'root.mod1.u');
model.physics('spf').feature('inl1').set('U0in', 'Re_cyl*spf.mu/(spf.rho*2*R_cyl*0.001)');
model.physics('spf').feature('out1').set('BoundaryCondition', 'Pressure');
model.physics('chds').feature('cdm1').set('u_src', 'root.mod1.u');
model.physics('chds').feature('cdm1').set('D_0', {'D_drug_fluid'; '0'; '0'; '0'; 'D_drug_fluid'; '0'; '0'; '0'; 'D_drug_fluid'});
model.physics('chds').feature('cdm1').set('minput_concentration_src', 'root.mod1.c');
model.physics('chds').feature('conc1').set('c0', 'c2*K_part');
model.physics('chds').feature('conc1').set('species', '1');
model.physics('chds2').prop('EquationForm').set('form', 'Transient');
model.physics('chds2').prop('Convection').set('Convection', '0');
model.physics('chds2').feature('cdm1').set('D_0', {'1e-4*D_pex*exp(-(A_vd+B_vd*(c2*M_drug/rho_poly))/(C_vd+D_vd*(c2*M_drug/rho_poly)))'; '0'; '0'; '0'; '1e-4*D_pex*exp(-(A_vd+B_vd*(c2*M_drug/rho_poly))/(C_vd+D_vd*(c2*M_drug/rho_poly)))'; '0'; '0'; '0'; '1e-4*D_pex*exp(-(A_vd+B_vd*(c2*M_drug/rho_poly))/(C_vd+D_vd*(c2*M_drug/rho_poly)))'});
model.physics('chds2').feature('cdm1').set('minput_concentration_src', 'root.mod1.c2');
model.physics('chds2').feature('conc2').set('c0', 'c3*K_r');
model.physics('chds2').feature('conc2').set('species', '1');
model.physics('chds2').feature('conc3').set('c0', 'c/K_part');
model.physics('chds2').feature('conc3').set('species', '1');
model.physics('chds3').prop('EquationForm').set('form', 'Transient');
model.physics('chds3').prop('Convection').set('Convection', '0');
model.physics('chds3').feature('cdm1').set('u_src', 'root.mod1.u');
model.physics('chds3').feature('cdm1').set('D_0', {'D_drug_polyI'; '0'; '0'; '0'; 'D_drug_polyI'; '0'; '0'; '0'; 'D_drug_polyI'});
model.physics('chds3').feature('init1').set('c3', 'conc_init');
model.physics('chds3').feature('conc1').set('c0', 'c2/K_r');
model.physics('chds3').feature('conc1').set('species', '1');
model.physics('chds3').feature('open1').selection.active(false);
model.physics('chds3').feature('open1').set('c0', 'c2');
model.physics('chds3').feature('open1').active(false);

model.mesh('mesh1').feature('size').set('hauto', 7);
model.mesh('mesh1').feature('size1').set('table', 'cfd');
model.mesh('mesh1').feature('size1').set('hauto', 8);
model.mesh('mesh1').feature('size2').set('table', 'cfd');
model.mesh('mesh1').feature('size2').set('hauto', 6);
model.mesh('mesh1').feature('bl1').feature('blp1').set('blnlayers', '2');
model.mesh('mesh1').feature('bl1').feature('blp1').set('blhminfact', '5');
model.mesh('mesh1').run;

model.frame('material1').sorder(1);

model.study.create('std2');
model.study('std2').feature.create('stat', 'Stationary');
model.study('std2').feature.create('time', 'Transient');

model.sol.create('sol7');
model.sol('sol7').study('std2');
model.sol('sol7').attach('std2');
model.sol('sol7').feature.create('st1', 'StudyStep');
model.sol('sol7').feature.create('v1', 'Variables');
model.sol('sol7').feature.create('s1', 'Stationary');
model.sol('sol7').feature('s1').feature.create('fc1', 'FullyCoupled');
model.sol('sol7').feature('s1').feature.create('d1', 'Direct');
model.sol('sol7').feature('s1').feature.remove('fcDef');
model.sol('sol7').feature.create('su1', 'StoreSolution');
model.sol('sol7').feature.create('st2', 'StudyStep');
model.sol('sol7').feature.create('v2', 'Variables');
model.sol('sol7').feature.create('t1', 'Time');
model.sol('sol7').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol7').feature('t1').feature.create('d1', 'Direct');
model.sol('sol7').feature('t1').feature.remove('fcDef');

model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset.create('cpl1', 'CutPlane');
model.result.dataset.remove('dset2');
model.result.numerical.create('int1', 'IntSurface');
model.result.numerical('int1').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').feature.create('surf1', 'Surface');
model.result.create('pg2', 'PlotGroup2D');
model.result('pg2').feature.create('con1', 'Contour');
model.result.create('pg3', 'PlotGroup3D');
model.result('pg3').feature.create('surf1', 'Surface');
model.result.create('pg4', 'PlotGroup2D');
model.result('pg4').feature.create('surf1', 'Surface');
model.result.create('pg5', 'PlotGroup3D');
model.result('pg5').feature.create('surf1', 'Surface');
model.result.create('pg6', 'PlotGroup2D');
model.result('pg6').feature.create('surf1', 'Surface');
model.result.create('pg7', 'PlotGroup3D');
model.result('pg7').feature.create('surf1', 'Surface');
model.result.create('pg8', 'PlotGroup2D');
model.result('pg8').feature.create('surf1', 'Surface');
model.result.create('pg9', 'PlotGroup3D');
model.result('pg9').feature.create('surf1', 'Surface');

model.study('std2').feature('stat').set('activate', {'spf' 'on' 'chds' 'off' 'chds2' 'off' 'chds3' 'off'});
model.study('std2').feature('time').set('tlist', 'range(0,time_step,time_final)');
model.study('std2').feature('time').set('activate', {'spf' 'off' 'chds' 'on' 'chds2' 'on' 'chds3' 'on'});

model.sol('sol7').attach('std2');
model.sol('sol7').feature('st1').name('Compile Equations: Stationary');
model.sol('sol7').feature('st1').set('studystep', 'stat');
model.sol('sol7').feature('v1').set('control', 'stat');
model.sol('sol7').feature('v1').feature('mod1_c3').set('solvefor', false);
model.sol('sol7').feature('v1').feature('mod1_c2').set('solvefor', false);
model.sol('sol7').feature('v1').feature('mod1_c').set('solvefor', false);
model.sol('sol7').feature('s1').set('control', 'stat');
model.sol('sol7').feature('s1').feature('fc1').set('initstep', '0.01');
model.sol('sol7').feature('s1').feature('fc1').set('minstep', '1.0E-6');
model.sol('sol7').feature('s1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol7').feature('st2').name('Compile Equations: Time Dependent (2)');
model.sol('sol7').feature('st2').set('studystep', 'time');
model.sol('sol7').feature('v2').set('control', 'time');
model.sol('sol7').feature('v2').set('initmethod', 'sol');
model.sol('sol7').feature('v2').set('initsol', 'sol7');
model.sol('sol7').feature('v2').set('notsolmethod', 'sol');
model.sol('sol7').feature('v2').set('notsol', 'sol7');
model.sol('sol7').feature('v2').set('notsoluse', 'sol8');
model.sol('sol7').feature('v2').feature('mod1_u').set('solvefor', false);
model.sol('sol7').feature('v2').feature('mod1_p').set('solvefor', false);
model.sol('sol7').feature('t1').set('control', 'time');
model.sol('sol7').feature('t1').set('tlist', 'range(0,time_step,time_final)');
model.sol('sol7').feature('t1').set('maxorder', '2');
model.sol('sol7').feature('t1').feature('fc1').set('maxiter', '5');
model.sol('sol7').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol7').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol7').runAll;

model.result.dataset('rev1').set('startangle', '-90');
model.result.dataset('rev1').set('revangle', '225');
model.result.dataset('cpl1').set('quickplane', 'xy');
model.result.dataset('cpl1').set('quickz', '0.02');
model.result.numerical('int1').set('data', 'cpl1');
model.result.numerical('int1').set('looplevelinput', {'interp'});
model.result.numerical('int1').set('showinterp', {'on' 'off' 'off'});
model.result.numerical('int1').set('interp', {'range(0,sStep,sEnd) range(sEnd,lStep,time_final)'});
model.result.numerical('int1').set('expr', 'chds.tfluxMag_c');
model.result.numerical('int1').set('unit', 'mol/s');
model.result.numerical('int1').set('descr', 'Total flux magnitude');
model.result('pg1').name('Velocity (spf)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg2').name('Pressure (spf)');
model.result('pg2').set('frametype', 'spatial');
model.result('pg2').feature('con1').set('expr', 'p');
model.result('pg2').feature('con1').set('unit', 'Pa');
model.result('pg2').feature('con1').set('descr', 'Pressure');
model.result('pg2').feature('con1').set('number', '40');
model.result('pg3').name('Velocity, 3D (spf)');
model.result('pg3').set('frametype', 'spatial');
model.result('pg4').name('Concentration (chds)');
model.result('pg4').feature('surf1').set('expr', 'c');
model.result('pg4').feature('surf1').set('unit', 'mol/m^3');
model.result('pg4').feature('surf1').set('descr', 'Concentration');
model.result('pg5').name('Concentration, 3D (chds)');
model.result('pg5').feature('surf1').set('expr', 'c');
model.result('pg5').feature('surf1').set('unit', 'mol/m^3');
model.result('pg5').feature('surf1').set('descr', 'Concentration');
model.result('pg6').name('Concentration (chds2)');
model.result('pg6').feature('surf1').set('expr', 'c2');
model.result('pg6').feature('surf1').set('unit', 'mol/m^3');
model.result('pg6').feature('surf1').set('descr', 'Concentration');
model.result('pg7').name('Concentration, 3D (chds2)');
model.result('pg7').feature('surf1').set('expr', 'c2');
model.result('pg7').feature('surf1').set('unit', 'mol/m^3');
model.result('pg7').feature('surf1').set('descr', 'Concentration');
model.result('pg8').name('Concentration (chds3)');
model.result('pg8').feature('surf1').set('expr', 'c3');
model.result('pg8').feature('surf1').set('unit', 'mol/m^3');
model.result('pg8').feature('surf1').set('descr', 'Concentration');
model.result('pg9').name('Concentration, 3D (chds3)');
model.result('pg9').feature('surf1').set('expr', 'c3');
model.result('pg9').feature('surf1').set('unit', 'mol/m^3');
model.result('pg9').feature('surf1').set('descr', 'Concentration');

out = model;
