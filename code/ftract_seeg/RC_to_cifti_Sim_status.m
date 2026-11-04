SimPos = 188; %
time_frame_n = 20;
n = size(I_node_max(1:360,:),1);
for j = 0:0
    %c = cifti_read('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
    c = cifti_read('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
% 
%         rho_amp_sig  = rho_amp;
% rho_prob_sig = rho_prob;
% rho_amp_sig(p_amp > 0.05)   = NaN;
% rho_prob_sig(p_prob > 0.05) = NaN;
%        I_node_max_sim = rho_prob_sig;
        %I_node_max_sim = I_node_max
      %I_node_max_sim = Prob_sorted(SimPos,:);
      %I_node_max_sim = Amp_sorted(SimPos,:);
      I_node_max_sim =(I_node_max(1:360,SimPos));
      I_node_max_sim =log10(I_node_max(1:360,SimPos)) + abs(min(log10(I_node_max(1:360,SimPos))));
      dist_vec  = D_cor(SimPos,1:360);   
      near_mask = dist_vec < 20; 
      near_mask(SimPos) = true;   
      I_node_max_sim(near_mask) = 0;

        clabel0 = c.cdata;
        clabel1  = zeros(size(clabel0,1),1);
        for i = 1:size(clabel0,1)
            if clabel0(i) == 0
                continue
            else
                %clabel1(i) = I_node_max_sim((clabel0(i)));
                clabel1(i) = I_node_max_sim((clabel0(i)));

            end
        end

        c.cdata = clabel1;
        filename_var = strcat('Test_RC_current_log',num2str(SimPos),'t',num2str(j),'.dscalar.nii');
        ciftisavereset(c,filename_var);

end

% SimPos = 265; %259 for dlpfc 188 for M1
% for SimPos = 255:275
%     time_frame_n = 20;
%     n = size(I_node_max(1:360,:),1);
%     for j = 0:0
%         %c = cifti_read('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
%         c = cifti_read('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
% 
% 
%         Stim_site_Node = zeros(360,1);
%         Stim_site_Node(SimPos) = 1;
%         I_node_max_sim = Stim_site_Node;
%         clabel0 = c.cdata;
%         clabel1  = zeros(size(clabel0,1),1);
%         for i = 1:size(clabel0,1)
%             if clabel0(i) == 0
%                 continue
%             else
%                 %clabel1(i) = I_node_max_sim((clabel0(i)));
%                 clabel1(i) = I_node_max_sim((clabel0(i)));
% 
%             end
%         end
% 
%         c.cdata = clabel1;
%         filename_var = strcat('Test_Stim_loc',num2str(SimPos),'t',num2str(j),'.dscalar.nii');
%         ciftisavereset(c,filename_var);
% 
%     end
% end
% 
SimPos = 188; %
time_frame_n = 20;
n = size(I_node_max(1:360,:),1);
for j = 0:20
    %c = cifti_read('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR_Tian_Subcortex_S1.dlabel.nii');
    c = cifti_read('Q1-Q6_RelatedValidation210.CorticalAreas_dil_Final_Final_Areas_Group_Colors.32k_fs_LR.dlabel.nii');
% 
%         rho_amp_sig  = rho_amp;
% rho_prob_sig = rho_prob;
% rho_amp_sig(p_amp > 0.05)   = NaN;
% rho_prob_sig(p_prob > 0.05) = NaN;
if j == 0 
    time_frame_j = 1;
else
    time_frame_j = 500./20 .*j;
end
       I_node_t_sim = -I_node_t(SimPos,:,time_frame_j);
       I_node_t_sim(I_node_t_sim <= 0) = NaN;
       I_node_max_sim =(I_node_t_sim(1:360));
       I_node_max_sim =log10(I_node_max_sim) + abs(min(log10(I_node_max_sim)));
        %I_node_max_sim = I_node_max
      %I_node_max_sim = Prob_sorted(SimPos,:);
      %I_node_max_sim = Amp_sorted(SimPos,:);
      %I_node_max_sim =log10(I_node_max(1:360,SimPos)) + abs(min(log10(I_node_max(1:360,SimPos))));
      % dist_vec  = D_cor(SimPos,1:360);   
      % near_mask = dist_vec < 20; 
      % near_mask(SimPos) = true;   
      % I_node_max_sim(near_mask) = 0;

        clabel0 = c.cdata;
        clabel1  = zeros(size(clabel0,1),1);
        for i = 1:size(clabel0,1)
            if clabel0(i) == 0
                continue
            else
                %clabel1(i) = I_node_max_sim((clabel0(i)));
                clabel1(i) = I_node_max_sim((clabel0(i)));

            end
        end

        c.cdata = clabel1;
        filename_var = strcat('Test_RC_dynamic_current_pos_log',num2str(SimPos),'t',num2str(j),'.dscalar.nii');
        ciftisavereset(c,filename_var);

end