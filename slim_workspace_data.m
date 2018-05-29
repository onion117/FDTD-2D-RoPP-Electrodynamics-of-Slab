close all
clear all
clc

addpath EM_functions
addpath material_data
addpath output_data

sim_casev='negative_matched';


models_fig=[ 'AB ';'MN ';'Chu';'EL ';'AMP'];

models_fig = cellstr(models_fig);

for model_j=1:length(models_fig)

    model=char(models_fig(model_j));
    load(strcat(model,'_Workspace_data_',sim_casev)); 

    clear A3 A4 A4 A5 Bz Bz_at_x Bz_at_y C1 C2 Dx Dx_n_prev
    clear Dy Dy_n_prev Ex Ex_n_prev Ey Ey_n_prev G_x G_x_n_prev
    clear G_y G_y_n_prev Hz_at_z Hz_at_y Hz_n_prev Jex Jey Jmz
    clear Jmz_n_prev_fr Jmzd Jmzd_n_prev Jmzs Jxd Jxs_n_prev
    clear Jyd Jyd_n_prev Jys Jys_n_prev M1 M2 M3 Mz P1 P1Y P2 P2Y
    clear Px Pxd Pxs Py Pyd Pys Tx Ty X3 X4 X5 Y1 Y2 Y3 g_mech_x g_mech_y
    clear t1 t2 t3 t4
    
    model_type=char(models(sim_j));
save(strcat(pwd,'\output_data\',model_type,'_Workspace_data','_',sim_case));
   

    
        
end
disp(' Compression done')

