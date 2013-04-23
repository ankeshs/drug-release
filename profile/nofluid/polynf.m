function out = ploynf(Di, Do, Kr, c0, ar, cd, tf, dt, st)
%
% polynf.m
%
% Model exported on Mar 13 2013, 21:22 by COMSOL 4.3.0.184.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/users/che/bt/ankeshs/Academics/che496/Transient');

model.name('polnf1.mph');

model.param.set('Re_cyl', '1');
model.param.set('R_cyl', '.01');
%model.param.set('a_ratio', '0.4');
model.param.set('a_ratio', num2str(ar));

model.param.set('rho_fluid', '1080');
model.param.set('mu_fluid', '0.0089');
model.param.set('b_conc', '1');
model.param.set('D_drug_fluid', '3e-8');
%model.param.set('t_poly', '0.2');
model.param.set('t_poly', num2str(cd));
%model.param.set('D_drug_polyI', '1e-14');
model.param.set('D_drug_polyI', num2str(Di));

model.param.set('K_part', '0.1');
%model.param.set('D_drug_polyO', '1e-15*0.5');
model.param.set('D_drug_polyO', num2str(Do));

%model.param.set('conc_init', '100');
model.param.set('conc_init', num2str(c0));

%model.param.set('time_final', '6000');
model.param.set('time_final', num2str(tf));
%model.param.set('time_step', '10');
model.param.set('time_step', num2str(dt));
%model.param.set('itm_step', '100');
model.param.set('itm_step', num2str(st));
%model.param.set('K_r', '0.4');
model.param.set('K_r', num2str(Kr));

model.modelNode.create('mod1');

model.geom.create('geom1', 2);
model.geom('geom1').axisymmetric(true);
model.geom('geom1').lengthUnit('mm');
model.geom('geom1').feature.create('r1', 'Rectangle');
model.geom('geom1').feature.create('c1', 'Circle');
model.geom('geom1').feature.create('sq1', 'Square');
model.geom('geom1').feature.create('int1', 'Intersection');
model.geom('geom1').feature.create('c2', 'Circle');
model.geom('geom1').feature('r1').active(false);
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

model.physics.create('chds2', 'DilutedSpecies', 'geom1');
model.physics('chds2').identifier('chds2');
model.physics('chds2').field('concentration').field('c2');
model.physics('chds2').field('concentration').component({'c2'});
model.physics('chds2').selection.set([1]);
model.physics('chds2').feature.create('conc2', 'Concentration', 1);
model.physics('chds2').feature('conc2').selection.set([6 7]);
model.physics('chds2').feature.create('conc3', 'Concentration', 1);
model.physics('chds2').feature('conc3').selection.set([5 8]);
model.physics('chds2').feature.create('fl1', 'Fluxes', 1);
model.physics('chds2').feature('fl1').selection.set([6 7]);
model.physics.create('chds3', 'DilutedSpecies', 'geom1');
model.physics('chds3').identifier('chds3');
model.physics('chds3').field('concentration').field('c3');
model.physics('chds3').field('concentration').component({'c3'});
model.physics('chds3').selection.set([2]);
model.physics('chds3').feature.create('conc1', 'Concentration', 1);
model.physics('chds3').feature('conc1').selection.set([6 7]);

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftri1', 'FreeTri');

model.result.table.create('tbl1', 'Table');
model.result.table.create('tbl2', 'Table');
model.result.table.create('tbl3', 'Table');

model.view('view1').axis.set('xmin', '-0.0026127658784389496');
model.view('view1').axis.set('xmax', '0.006612766534090042');
model.view('view1').axis.set('ymin', '-0.004399999976158142');
model.view('view1').axis.set('ymax', '0.004399999976158142');

model.physics('chds2').prop('EquationForm').set('form', 'Transient');
model.physics('chds2').prop('Convection').set('Convection', '0');
model.physics('chds2').feature('cdm1').set('D_0', {'D_drug_polyO'; '0'; '0'; '0'; 'D_drug_polyO'; '0'; '0'; '0'; 'D_drug_polyO'});
model.physics('chds2').feature('cdm1').set('minput_concentration_src', 'root.mod1.c2');
model.physics('chds2').feature('conc2').set('c0', 'c3*K_r');
model.physics('chds2').feature('conc2').set('species', '1');
model.physics('chds2').feature('conc3').set('species', '1');
model.physics('chds2').feature('fl1').selection.active(false);
model.physics('chds2').feature('fl1').set('species', '1');
model.physics('chds2').feature('fl1').set('N0', 'N3');
model.physics('chds2').feature('fl1').active(false);
model.physics('chds3').prop('EquationForm').set('form', 'Transient');
model.physics('chds3').prop('Convection').set('Convection', '0');
model.physics('chds3').feature('cdm1').set('D_0', {'D_drug_polyI'; '0'; '0'; '0'; 'D_drug_polyI'; '0'; '0'; '0'; 'D_drug_polyI'});
model.physics('chds3').feature('init1').set('c3', 'conc_init');
model.physics('chds3').feature('conc1').set('c0', 'c2/K_r');
model.physics('chds3').feature('conc1').set('species', '1');

model.mesh('mesh1').feature('size').set('hauto', 2);
model.mesh('mesh1').run;

model.frame('material1').sorder(1);

model.result.table('tbl1').comments('Surface Integration 1 (abs(chds.tfluxz_c))');
model.result.table('tbl2').comments('Surface Integration 1 (c2chds2.tfluxMag_c2)');
model.result.table('tbl3').comments('Surface Integration 1 (chds2.tfluxMag_c2)');

model.study.create('std3');
model.study('std3').feature.create('time', 'Transient');

model.sol.create('sol1');
model.sol('sol1').study('std3');
model.sol('sol1').attach('std3');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature.create('t1', 'Time');
model.sol('sol1').feature('t1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('t1').feature.create('d1', 'Direct');
model.sol('sol1').feature('t1').feature.remove('fcDef');

model.result.dataset.create('rev1', 'Revolve2D');
model.result.dataset.create('edg1', 'Edge2D');
model.result.dataset('edg1').selection.set([5 8]);
model.result.dataset.create('rev2', 'Revolve2D');
model.result.dataset('rev2').set('data', 'edg1');
model.result.numerical.create('int1', 'IntSurface');
model.result.numerical('int1').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup2D');
model.result('pg1').feature.create('surf1', 'Surface');
model.result.create('pg2', 'PlotGroup3D');
model.result('pg2').feature.create('surf1', 'Surface');
model.result.create('pg3', 'PlotGroup2D');
model.result('pg3').feature.create('surf1', 'Surface');
model.result.create('pg4', 'PlotGroup3D');
model.result('pg4').feature.create('surf1', 'Surface');
model.result.create('pg5', 'PlotGroup1D');
model.result('pg5').set('probetag', 'none');
model.result('pg5').feature.create('tblp1', 'Table');
model.result('pg5').feature('tblp1').set('probetag', 'none');

model.study('std3').feature('time').set('tlist', 'range(0,time_step,time_final)');

model.sol('sol1').attach('std3');
model.sol('sol1').feature('st1').name('Compile Equations: Time Dependent');
model.sol('sol1').feature('st1').set('studystep', 'time');
model.sol('sol1').feature('v1').set('control', 'time');
model.sol('sol1').feature('t1').set('control', 'time');
model.sol('sol1').feature('t1').set('tlist', 'range(0,time_step,time_final)');
model.sol('sol1').feature('t1').set('maxorder', '2');
model.sol('sol1').feature('t1').feature('fc1').set('maxiter', '5');
model.sol('sol1').feature('t1').feature('fc1').set('jtech', 'once');
model.sol('sol1').feature('t1').feature('d1').set('linsolver', 'pardiso');
model.sol('sol1').runAll;

model.result.dataset('rev1').set('startangle', '-90');
model.result.dataset('rev1').set('revangle', '225');
model.result.numerical('int1').set('data', 'rev2');
model.result.numerical('int1').set('looplevelinput', {'interp'});
model.result.numerical('int1').set('interp', {'range(0,itm_step,time_final)'});
model.result.numerical('int1').set('table', 'tbl3');
model.result.numerical('int1').set('expr', 'chds2.tfluxMag_c2');
model.result.numerical('int1').set('unit', 'mol/s');
model.result.numerical('int1').set('descr', 'Total flux magnitude');
model.result.numerical('int1').setResult;
model.result('pg1').name('Concentration (chds2)');
model.result('pg1').setIndex('looplevel', '501', 0);
model.result('pg2').name('Concentration, 3D (chds2)');
model.result('pg2').setIndex('looplevel', '1', 0);
model.result('pg3').name('Concentration (chds3)');
model.result('pg3').setIndex('looplevel', '306', 0);
model.result('pg3').feature('surf1').set('expr', 'c3');
model.result('pg4').name('Concentration, 3D (chds3)');
model.result('pg4').setIndex('looplevel', '1', 0);
model.result('pg4').feature('surf1').set('expr', 'c3');
model.result('pg5').set('data', 'none');
model.result('pg5').set('xlabel', 'Time');
model.result('pg5').set('ylabel', 'Total flux magnitude (mol/s)');
model.result('pg5').set('xlabelactive', false);
model.result('pg5').set('ylabelactive', false);
model.result('pg5').feature('tblp1').set('table', 'tbl3');

out = model;
