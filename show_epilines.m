T_gt_cross=cross_matrix(T_gt);
E_gt=T_gt_cross*R_gt;
F_gt=invK'*E_gt*invK;
T_final_cross=cross_matrix(T_final);
E_final=T_final_cross*R_final;
F_final=invK'*E_final*invK;
T_algo1_cross=cross_matrix(T_algo1);
E_algo1=T_algo1_cross*R_algo1;
F_algo1=invK'*E_algo1*invK;
T_fivept_cross=cross_matrix(fiveptT);
E_fivept=T_fivept_cross*fiveptR;
F_fivept=invK'*E_fivept*invK;
fig=vgg_gui_F(image2,image1,F_final);